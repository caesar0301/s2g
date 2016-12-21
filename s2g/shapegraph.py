#!/usr/env/bin python
# encoding: utf-8
import fiona
import logging
import networkx as nx
import numpy as np
import progressbar
import sys
from itertools import product
from shapely.geometry import shape, Point, LineString, box, Polygon

from s2g.bonus import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

__all__ = ['ShapeGraph']


class EdgeSegment(object):
    def __init__(self, edge_key, distance, road_segment, line_index=None, line_id=None, cuts_key=None):
        self.edge_key = edge_key
        self.distance = distance
        self.road_segment = road_segment
        self.line_index = line_index
        self.line_id = line_id
        self.cuts_key = cuts_key

    def __str__(self):
        return '{0} [{1:.3f}, line_index {2}, line_id {3}, cuts {4}]'\
            .format(self.edge_key, self.distance, self.line_id, self.cuts_key)

class ShapeGraph(object):
    def __init__(self, geoms=None, shapefile=None, to_graph=True, resolution=1.0,
                 coord_type='lonlat', point_buffer=10e-5, line_id=None):
        """
        Declare and/or construct a graph from a shapefile
        :param geoms: input1, a list of shapely::LineString objects or (lon, lat) coordinates sequence
        :param shapefile: input2, path to shapefile; supported formats are consistent with fiona::open()
        :param to_graph: boolean, indicates if constructing the graph during initialization
        :param resolution: value in kilometers, the maximum road segment distance to become a graph edge
        :param coord_type: string, unused, indicates the type of coordinates
        :param point_buffer: value, the point tolerance (in Euclidean coordinates) to
            perform line linking or point projection
        :param line_id: string, the property field to represent the identify of a road, especially when
            multiple lines consist of a whole road
        """
        assert not (geoms is None and shapefile is None)
        assert not (geoms is not None and shapefile is not None)

        self.line_id_key = line_id
        self.line_ids = {}  # dict, line indexes to line_ids

        if shapefile is not None:
            self.geoms = []
            with fiona.open(shapefile) as source:
                for r in source:
                    s = shape(r['geometry'])
                    if s.geom_type == 'LineString':
                        if line_id is not None:
                            self.line_ids[len(self.geoms)] = r['properties'][line_id]
                        self.geoms.append(s)
                    else:
                        logging.warning('Misc geometry type encountered and omit: {0}'.format(s.geom_type))
        else:
            if line_id is not None:
                logging.warning('Plain coordinates sequence detected, "line_id" is ignored')
            self.geoms = [LineString(line) for line in geoms]

        # parameters
        self.resolution = resolution
        self.coord_type = coord_type
        self.point_buffer = point_buffer
        # components
        self.connected = False
        self.connectivity = []     # list, pairs of line indices
        self.major_components = [] # list of list, line indices

        # Global edge info. Note: node_* are redundant.
        # DO NOT use these to iterate over graph.
        self.node_ids = {}   # dict, node (lon, lat) to ids
        self.node_xy = {}    # dict, node ids to (lon, lat)
        self._edges = {}     # dict, key as self._edge_key, value as (distance,
        # road segment). For edge bridge, the original road segment is recorded.

        # line cuts info for the largest major component
        self.line_cuts = {}  # dict, line index to line cuts. Each cut in *cuts*
        # is a point index of line, including both ends
        self.nodes_counter = 0
        self.graph = nx.Graph()  # generated graph data
        if to_graph:
            self.gen_major_components()
            self.to_networkx()

    def _edge_key(self, p1, p2):
        return tuple(sorted([self.node_ids[p1], self.node_ids[p2]]))

    def _edge_key_by_nodes(self, n1, n2):
        return tuple(sorted([n1, n2]))

    def _register_edge(self, edge, dist, raw_segment, line_index=None, line_id=None, cuts=None):
        p1, p2 = edge
        if isinstance(p1, Point) and isinstance(p2, Point):
            p1 = list(p1.coords)[0]
            p2 = list(p2.coords)[0]
        elif isinstance(p1, tuple) and isinstance(p2, tuple):
            pass
        else:
            raise TypeError('The edge ends should be shapely::Point or 2-tuple')

        if p1 == p2:
            return

        if p1 not in self.node_ids:
            self.node_ids[p1] = self.nodes_counter
            self.nodes_counter += 1
            self.node_xy[self.node_ids[p1]] = p1

        if p2 not in self.node_ids:
            self.node_ids[p2] = self.nodes_counter
            self.nodes_counter += 1
            self.node_xy[self.node_ids[p2]] = p2

        n1 = self.node_ids[p1]
        n2 = self.node_ids[p2]
        if n1 <= n2:
            edge = (n1, n2)
            edge_cuts = cuts if cuts else None
        else:
            edge = (n2, n1)
            edge_cuts = (cuts[1], cuts[0]) if cuts else None

        if edge not in self._edges:
            es = EdgeSegment(edge, dist, raw_segment, line_index=line_index,
                             line_id=line_id, cuts_key=edge_cuts)
            self._edges[edge] = es

    def _remove_edge(self, p1, p2):
        if p1 not in self.node_ids or p2 not in self.node_ids:
            return None
        id1 = self.node_ids[p1]
        id2 = self.node_ids[p2]
        # TODO: remove nodes that do not belong to any edges
        # Now nodes are redundant in the graph.
        # del self.node_ids[p1]
        # del self.node_ids[p2]
        # del self.node_xy[id1]
        # del self.node_xy[id2]
        edge = self._edge_key_by_nodes(id1, id2)
        del self._edges[edge]

    def _update_cut(self, line_index, cut):
        if line_index not in self.line_cuts:
            self.line_cuts[line_index] = set()
        self.line_cuts[line_index].add(cut)

    def _validate_intersection(self, line_index, point):
        line = self.geoms[line_index]
        coords = list(line.coords)
        buffered_point = Point(point).buffer(self.point_buffer)
        touched = False
        if line.intersects(buffered_point):
            touched = True
            cut = point_projects_to_line(point, line)
            pp = coords[cut]
            if pp != point:
                d = great_circle_dist(pp, point)
                self._register_edge((pp, point), d, [pp, point])
            if 0 < cut < len(coords)-1:
                self._update_cut(line_index, cut)
        return touched

    def _validate_pairwise_connectivity(self, ainx, binx):
        a1, a2 = self.geoms[ainx].coords[0], self.geoms[ainx].coords[-1]
        b1, b2 = self.geoms[binx].coords[0], self.geoms[binx].coords[-1]
        return self._validate_intersection(ainx, b1) or \
            self._validate_intersection(ainx, b2) or \
            self._validate_intersection(binx, a1) or \
            self._validate_intersection(binx, a2)

    def road_segment_for_edge(self, edge):
        """
        Get original road segment in point sequence which is
        bridged by graph edges.
        :param edge: a tuple of two node ids
        :return: the point sequence of raw road segment
        """
        edge = self._edge_key_by_nodes(edge[0], edge[1])
        assert edge in self._edges
        return self._edges[edge].road_segment

    def gen_major_components(self):
        """
        Merge individual lines into groups of touching lines.
        :return: a tuple of (connected, connectivity, major_components)
        """
        lines = self.geoms
        L = len(lines)
        graph = nx.Graph()
        logging.info("Validating pair-wise line connections of raw shapefiles (total {0} lines)".format(L))
        neighbors = [(i, j) for i, j in product(range(0, L), range(0, L)) if j > i]
        with progressbar.ProgressBar(max_value=len(neighbors)) as bar:
            for k in range(0, len(neighbors)):
                bar.update(k+1)
                i, j = neighbors[k]
                if self._validate_pairwise_connectivity(i, j):
                    graph.add_edge(i, j)

        # Validate lines connectivity:
        # Pairs of lines which are connected (touched). These line pairs
        # are not guaranteed to be in the same major component.
        self.connectivity = [i for i in graph.edges() if i[0] != i[1]]

        # major connected components
        cc = nx.algorithms.components.connected_components(graph)
        self.major_components = sorted([i for i in cc], key=len, reverse=True)

        # predict if all lines are strongly connected
        self.connected = True if len(self.major_components) == 1 else False

        # statistics
        print("Major components statistics:")
        print("\tTotal components: {0}".format(len(self.major_components)))
        size = [len(c) for c in self.major_components]
        print("\tComponent size: max {0}, median {1}, min {2}, average {3}"\
              .format(np.max(size), np.median(size), np.min(size), np.average(size)))
        print("\tTop comp. sizes: {0}".format(' '.join([str(i) for i in size[0:10]])))

        return self.connected, self.connectivity, self.major_components

    def largest_component(self):
        """
        Get the largest connected component.
        :return: a list of lines consists of the largest component
        """
        return list(self.major_components[0])

    def to_networkx(self):
        """Convert the major component to graph of NetworkX.
        """
        if len(self.major_components) is None:
            self.gen_major_components()

        major = list(self.major_components[0])
        logging.info('Processing the largest component with {0} lines'.format(len(major)))

        # do line cutting
        logging.info('Cutting lines with specific resolution = {0} km'.format(self.resolution))
        with progressbar.ProgressBar(max_value=len(major)) as bar:
            for i in range(len(major)):
                bar.update(i+1)
                line_index = major[i]
                line = self.geoms[line_index]
                line_id = self.line_ids.get(line_index)
                coords = list(line.coords)
                intersects = self.line_cuts[line_index] if line_index in self.line_cuts else set()
                cuts, dist = cut_line(line, self.resolution, intersects)

                # record cut-line info
                [self._update_cut(line_index, c) for c in cuts]

                # record edges info
                for j in range(1, len(cuts)):
                    scut = cuts[j - 1]
                    ecut = cuts[j]
                    self._register_edge((coords[scut], coords[ecut]),
                                        dist[j], coords[scut:ecut+1],
                                        line_index, line_id, (scut, ecut))

        # assemble graph
        for edge, segment in self._edges.items():
            self.graph.add_edge(edge[0], edge[1], weight=segment.distance)

        # validate connectivity of generated graph
        cc = nx.algorithms.components.connected_components(self.graph)
        mc = sorted([i for i in cc], key=len, reverse=True)
        assert len(mc) == 1

        logging.info('Graph created with {0} nodes, {1} edges'.format(
            len(self.graph.edges()), len(self.graph.nodes())))

        return self.graph

    def point_projects_to_node(self, point, distance_tolerance=0.01):
        """
        Return the nearest (in sense of great circle distance) graph node to given point.
        :param point: a tuple of (lon, lat) or shapely Point
        :param distance_tolerance: approximated buffer range in kilometers
        :return: (nearest_node_id, nearest_node_lonlat) or None
        """
        p = Point(point)
        nearest = None
        min_dist = -1
        for nid in self.graph.nodes():
            dist = great_circle_dist(p.coords[0], self.node_xy[nid])
            if dist > distance_tolerance:
                continue
            if min_dist < 0:
                min_dist = dist
                nearest = nid
            else:
                if dist < min_dist:
                    min_dist = dist
                    nearest = nid
        if nearest is None:
            return None
        else:
            return nearest, self.node_xy[nearest]

    def point_projects_to_edges(self, point, distance_tolerance=0.01):
        """
        Project a point to graph edges considering specific distance tolerance.
        Note the tolerance is measured by great circle distance to point per se.
        :param point: a shapely Point instance or (lon, lat) tuple
        :param distance_tolerance: tolerance of distance in km
        :return: a list of projected edges, reversely sorted by offsets.
        """
        point_buffer = distance_to_buffer(distance_tolerance)
        p_buf = Point(point).buffer(point_buffer)
        projected_edges = []
        projected_segments = []
        major = self.largest_component()
        for i in range(0, len(major)):
            line_index = major[i]
            line = self.geoms[line_index]
            if line.intersects(p_buf):
                cuts = sorted(self.line_cuts[line_index])
                for j in range(1, len(cuts)):
                    sinx = cuts[j-1]
                    einx = cuts[j]
                    segment = line.coords[sinx:einx+1]
                    ls = LineString(segment)
                    if ls.intersects(p_buf):
                        edge = self._edge_key(segment[0], segment[-1])
                        offset = ls.distance(Point(point)) # no buffer
                        projected_edges.append((edge, offset))
                        projected_segments.append(segment)
        result = sorted(projected_edges, key=lambda x: x[1], reverse=True)
        edges = list(set([i[0] for i in result]))
        return edges, projected_segments

    def subgraph_within_box(self, bounding_box):
        """
        Extract a subgraph bouneded by a box.
        :param bounding_box: the bounding coordinates in
            (minx, miny, maxx, maxy) or a Polygon instance
        :return: a subgraph of nx.Graph
        """
        if isinstance(bounding_box, Polygon):
            bbox = bounding_box
        else:
            bbox = box(bounding_box[0], bounding_box[1],
                       bounding_box[2], bounding_box[3])
        nbunch = set()
        for edge in self.graph.edges():
            s, e = edge
            if bbox.intersects(LineString([self.node_xy[s], self.node_xy[e]])):
                nbunch.add(s)
                nbunch.add(e)
        return self.graph.subgraph(nbunch)


if __name__ == '__main__':
    pass