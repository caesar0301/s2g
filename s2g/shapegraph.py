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


class ShapeGraph(object):
    def __init__(self, geoms=None, shapefile=None, to_graph=True,
                 resolution=1.0, coord_type='lonlat', point_buffer=10e-5):
        # raw lines
        assert not (geoms is None and shapefile is None)
        assert not (geoms is not None and shapefile is not None)
        if shapefile is not None:
            self.geoms = []
            with fiona.open(shapefile) as source:
                for r in source:
                    s = shape(r['geometry'])
                    if s.geom_type == 'LineString':
                        self.geoms.append(s)
                    else:
                        logging.warning('Misc geometry type encountered and omit: {0}'.format(s.geom_type))
        else:
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
        self._edges = {}     # dict, key as self._edge_key
        # value as (distance, road segment).
        # For edge bridge, the original road segment is recorded.

        # line cuts info for the largest major component
        self.line_cuts = {}  # dict, line index to line (cuts, distances)
        # each cut in *cuts* is a point index of line, including both ends
        # each value in *distances* is the great circle distance of related segment
        self.nodes_counter = 0
        self.graph = nx.Graph()  # generated graph data
        if to_graph:
            self.gen_major_components()
            self.to_networkx()

    def _edge_key(self, p1, p2):
        return tuple(sorted([self.node_ids[p1], self.node_ids[p2]]))

    def _edge_key_by_nodes(self, n1, n2):
        return tuple(sorted([n1, n2]))

    def _register_edge(self, p1, p2, dist, raw_segment):
        assert isinstance(p1, Point) or len(p1) == 2
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
        edge = self._edge_key(p1, p2)
        if edge not in self._edges:
            self._edges[edge] = (dist, raw_segment)

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
        assert edge in self._edges
        del self._edges[edge]

    def _update_cut(self, lid, cut):
        if lid not in self.line_cuts:
            self.line_cuts[lid] = set()
        self.line_cuts[lid].add(cut)

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
                self._register_edge(pp, point,
                                    great_circle_dist(pp, point),
                                    [pp, point])
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

    def road_segment_from_edge(self, edge):
        """
        Get original road segment in point sequence which is
        bridged by graph edges.
        :param edge: a tuple of two node ids
        :return: the point sequence of raw road segment
        """
        edge = self._edge_key_by_nodes(edge[0], edge[1])
        assert edge in self._edges
        return self._edges[edge][1]

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
                lid = major[i]
                line = self.geoms[lid]
                coords = list(line.coords)
                fixed_cuts = self.line_cuts[lid] if lid in self.line_cuts else set()
                cuts, dist = cut_line(line, self.resolution, fixed_cuts)

                # record cut-line info
                [self._update_cut(lid, c) for c in cuts]

                # record edges info
                for j in range(1, len(cuts)):
                    sinx = cuts[j - 1]
                    einx = cuts[j]
                    self._register_edge(coords[sinx],
                                        coords[einx],
                                        dist[j],
                                        coords[sinx:einx+1])

        # assemble graph
        for edge, dist in self._edges.items():
            self.graph.add_edge(edge[0], edge[1], weight=dist)

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
        return [i[0] for i in result], projected_segments

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