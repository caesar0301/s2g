#!/usr/env/bin python
# encoding: utf-8
import logging
import numpy as np
import sys
from itertools import product

import fiona
import networkx as nx
import progressbar
from shapely.geometry import shape, Point, LineString, box, Polygon

from s2g.bonus import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

__all__ = ['ShapeGraph', 'EdgeInfo']


class EdgeInfo(object):
    def __init__(self, edge, distance, line_segment, line_index=None, cuts=None):
        self.edge = edge
        self.distance = distance
        self.line_segment = line_segment
        self.line_index = line_index
        self.cuts = cuts

    def __str__(self):
        return 'EdgeInfo: {0} [{1:.3f}, line_index {2}, cuts {3}]' \
            .format(self.edge, self.distance, self.line_index, self.cuts)


class LineInfo(object):
    def __init__(self, index, coords, props, is_major, cuts):
        self.index = index
        self.coords = coords
        self.props = props
        self.is_major = is_major
        self.cuts = cuts

    def __str__(self):
        return 'LineInfo: {0} [props {1}, is_major {2}, cuts {3}, coords {4}]'.format(
            self.index, self.props, self.is_major, self.cuts, self.coords
        )


class ShapeGraph(object):
    def __init__(self, geoms=None, shapefile=None, to_graph=True, resolution=1.0,
                 coord_type='lonlat', point_buffer=10e-5, properties=list(),
                 geom_count=None):
        """
        Declare and/or construct a graph from a shapefile
        :param geoms: input1, a list of shapely::LineString objects or (lon, lat) coordinates sequence
        :param shapefile: input2, path to shapefile; supported formats are consistent with fiona::open()
        :param to_graph: boolean, indicates if constructing the graph during initialization
        :param resolution: value in kilometers, the maximum road segment distance to become a graph edge
        :param coord_type: string, unused, indicates the type of coordinates
        :param point_buffer: value, the point tolerance (in Euclidean coordinates) to
            perform line linking or point projection
        :param properties: list, the property fields of a road
        """
        assert not (geoms is None and shapefile is None)
        assert not (geoms is not None and shapefile is not None)

        self.geoms = []  # all input lines
        self._line_props = {}  # dict, line indexes to properties, all input lines
        self._line_cuts = {}  # dict, line index to cuts (point indexes) set, for major component and more

        if shapefile is not None:
            with fiona.open(shapefile) as source:
                for r in source:
                    s = shape(r['geometry'])
                    if s.geom_type == 'LineString':
                        props = {}
                        for p in properties:
                            props[p] = r['properties'][p]
                        self._line_props[len(self.geoms)] = props
                        self.geoms.append(s)
                        # for debugging
                        if geom_count is not None and len(self.geoms) == geom_count:
                            break
                    else:
                        logging.warning('Misc geometry type encountered and omit: {0}'.format(s.geom_type))
        else:
            if len(properties):
                logging.warning('Plain coordinates sequence detected, properties are ignored: {0}'.format(properties))
            [self.geoms.append(LineString(i)) for i in geoms]

        # Parameters
        self.resolution = resolution
        self.coord_type = coord_type
        self.point_buffer = point_buffer

        # Components
        self.connected = False
        self.connectivity = []  # list, pairs of line indices
        self.major_components = []  # list of list, line indices

        # Global edge info. Note: node_* are redundant.
        # DO NOT use these to iterate over graph.
        self.node_ids = {}  # dict, node (lon, lat) to ids, for major component
        self.node_xy = {}  # dict, node ids to (lon, lat), for major component
        self._edges = {}  # dict, key as self._edge_key, value as EdgeSegment, for major component
        self._pseudo_edges = []  # for major component and more
        self.nodes_counter = 0
        self.graph = nx.Graph()  # generated graph data

        if to_graph:
            self.gen_major_components()
            self.to_networkx()

    def edge_key(self, p1, p2=None):
        if p2 is None:
            p1, p2 = p1
        return tuple(sorted([self.node_ids[p1], self.node_ids[p2]]))

    @staticmethod
    def edge_key_nodes(n1, n2=None):
        node1, node2 = (n1, n2)
        if n2 is None:
            node1, node2 = n1
        return tuple(sorted([node1, node2]))

    def _register_node(self, p):
        if isinstance(p, Point):
            p = list(p.coords)[0]
        elif isinstance(p, tuple):
            pass
        else:
            raise TypeError('The point should be shapely::Point or a 2-tuple')

        if p not in self.node_ids:
            nid = self.nodes_counter
            self.node_ids[p] = nid
            self.node_xy[nid] = p
            self.nodes_counter += 1

        return self.node_ids[p]

    def _register_edge(self, edge, dist, raw_segment, line_index=None, cuts=None):
        n1 = self._register_node(edge[0])
        n2 = self._register_node(edge[1])

        if n1 < n2:
            edge = (n1, n2)
            edge_cuts = tuple(cuts) if cuts else None
        elif n1 > n2:
            edge = (n2, n1)
            edge_cuts = (cuts[1], cuts[0]) if cuts else None
        else:
            return None

        if edge[0] == 1186 and edge[1] == 1206:
            print 'bingo', line_index, cuts, edge

        if edge not in self._edges:
            es = EdgeInfo(edge, dist, raw_segment, line_index, edge_cuts)
            self._edges[edge] = es
        return edge

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
        edge = self.edge_key_nodes(id1, id2)
        del self._edges[edge]

    def _update_cut(self, line_index, cut):
        if line_index not in self._line_cuts:
            self._line_cuts[line_index] = set()
        self._line_cuts[line_index].add(cut)

    def validate_intersection(self, line_index, point):
        """
        Validate the intersection between given line and point.
        A spatial tolerance is considered to confirm the validation.
        :param line_index: index of specific line
        :param point: validated point
        :return: coordinates of intersected point on line, otherwise None
        """
        line = self.geoms[line_index]
        coords = list(line.coords)
        buffered_point = Point(point).buffer(self.point_buffer)
        touched = None
        if line.intersects(buffered_point):
            cut = point_projects_to_line(point, line)
            touched = coords[cut]
            self._update_cut(line_index, cut)
        return touched

    def validate_pairwise_connectivity(self, ainx, binx):
        """
        Validate the connectivity between two lines.

        In this algorithm, line A's connectivity to line B is confirmed
        by the validation of either of A's ends intersecting with B.
        :param ainx: index of line A
        :param binx: index of line B
        :return: boolean
        """
        a1, a2 = self.geoms[ainx].coords[0], self.geoms[ainx].coords[-1]
        b1, b2 = self.geoms[binx].coords[0], self.geoms[binx].coords[-1]
        valid = False

        touched = self.validate_intersection(ainx, b1)
        if touched is not None:
            self._pseudo_edges.append([(ainx, touched), (binx, b1)])
            valid = True

        touched = self.validate_intersection(ainx, b2)
        if touched is not None:
            self._pseudo_edges.append([(ainx, touched), (binx, b2)])
            valid = True

        touched = self.validate_intersection(binx, a1)
        if touched is not None:
            self._pseudo_edges.append([(binx, touched), (ainx, a1)])
            valid = True

        touched = self.validate_intersection(binx, a2)
        if touched is not None:
            self._pseudo_edges.append([(binx, touched), (ainx, a2)])
            valid = True

        return valid

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
                bar.update(k + 1)
                i, j = neighbors[k]
                if self.validate_pairwise_connectivity(i, j):
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
        if len(self.major_components) == 0:
            size = [0]
        else:
            size = [len(c) for c in self.major_components]
        print("\tComponent size: max {0}, median {1}, min {2}, average {3}" \
              .format(np.max(size), np.median(size), np.min(size), np.average(size)))
        print("\tTop comp. sizes: {0}".format(' '.join([str(i) for i in size[0:10]])))

        return self.connected, self.connectivity, self.major_components

    def major_component(self):
        """
        Get the largest self-connected component.
        :return: a list of lines consists of the largest component
        """
        if len(self.major_components) > 0:
            return list(self.major_components[0])
        else:
            return None

    def to_networkx(self):
        """Convert the major component to graph of NetworkX.
        """
        if len(self.major_components) is None:
            self.gen_major_components()

        major = self.major_component()
        if not major:
            return

        logging.info('Processing the largest component with {0} lines'.format(len(major)))

        # do line cutting
        logging.info('Cutting lines with specific resolution = {0} km'.format(self.resolution))
        with progressbar.ProgressBar(max_value=len(major)) as bar:
            for i in range(len(major)):
                bar.update(i + 1)
                line_index = major[i]
                line = self.geoms[line_index]
                coords = list(line.coords)
                intersects = self.line_cuts(line_index)
                cuts, dist = cut_line(line, self.resolution, intersects if intersects else set())

                # record cut-line info
                [self._update_cut(line_index, c) for c in cuts]

                # record edges info
                for j in range(1, len(cuts)):
                    scut = cuts[j - 1]
                    ecut = cuts[j]
                    self._register_edge((coords[scut], coords[ecut]),
                                        dist[j], coords[scut:ecut + 1],
                                        line_index, (scut, ecut))

        logging.info('Adding pseudo edges to eliminate gaps between edges')
        for pair in self._pseudo_edges:
            ainx, pa = pair[0]
            binx, pb = pair[1]
            if ainx in major and binx in major and pa != pb:
                d = great_circle_dist(pa, pb)
                self._register_edge((pa, pb), d, [pa, pb])

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
        major = self.major_component()
        for i in range(0, len(major)):
            line_index = major[i]
            line = self.geoms[line_index]
            if line.intersects(p_buf):
                cuts = self.line_cuts(line_index)
                if cuts is None:
                    continue
                for j in range(1, len(cuts)):
                    sinx = cuts[j - 1]
                    einx = cuts[j]
                    segment = line.coords[sinx:einx + 1]
                    ls = LineString(segment)
                    if ls.intersects(p_buf):
                        edge = self.edge_key(segment[0], segment[-1])
                        offset = ls.distance(Point(point))  # no buffer
                        projected_edges.append((edge, offset))
                        projected_segments.append(segment)
        result = sorted(projected_edges, key=lambda x: x[1], reverse=True)
        edges = list(set([i[0] for i in result]))
        return edges, projected_segments

    def subgraph_within_box(self, bounding_box):
        """
        Extract a subgraph bounded by a box.
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

    def lines_within_box(self, bounding_box):
        """
        Get all lines selected by a bounding box. Note that the whole line
        is selected when it partially falls into the box.
        :param bounding_box: the bounding coordinates in
            (minx, miny, maxx, maxy) or a Polygon instance
        :return: list, of LineInfo
        """
        subgraph = self.subgraph_within_box(bounding_box)
        lines = set()
        for edge in subgraph.edges():
            info = self.edge_info(edge)
            if info.line_index:
                lines.add(info.line_index)
        return [self.line_info(i) for i in lines]

    def line_info(self, line_index):
        """
        Get basic information of a line
        :param line_index: index of line in input
        :return: instance of LineInfo
        """
        return LineInfo(
            index=line_index,
            coords=list(self.geoms[line_index].coords),
            props=self._line_props.get(line_index),
            is_major=line_index in self.major_component(),
            cuts=self.line_cuts(line_index)
        )

    def edge_info(self, edge):
        """
        Get information of road segment attached to specific edge
        :param edge: tuple, of edge nodes
        :return: instance of EdgeInfo
        """
        n1, n2 = edge
        edge_key = self.edge_key_nodes(n1, n2)
        return self._edges.get(edge_key)

    def major_lines_info(self):
        """
        Return the info of lines in the major component
        :return: dict, line index to info
        """
        result = {}
        major = self.major_component()
        for line_index in major:
            result[line_index] = self.line_info(line_index)
        return result

    def line_edge_sequence(self, line_index):
        """
        Get a sequence of edges belonging to specific line
        :param line_index: the index of specific line in input
        :return: list, of 2-tuple edges

        Note: the resulted edges keep their order along line coordinates sequence.
        In other words, the 2-tuple edge does NOT assure uniqueness as returned by
        self.edge_key or self.edge_key_nodes.
        """
        cuts = self.line_cuts(line_index)
        coords = list(self.geoms[line_index].coords)
        edges = []
        if cuts is not None:
            for i in range(1, len(cuts)):
                sn = self.node_ids[coords[cuts[i - 1]]]
                en = self.node_ids[coords[cuts[i]]]
                edge = (sn, en)
                edges.append(edge)
        return edges

    def line_props(self, line_index):
        """
        Get properties of specific line.
        :param line_index: the index of lines in input
        :return: dict, property names and values
        """
        return self._line_props.get(line_index)

    def line_cuts(self, line_index):
        """
        Get cut sequence of specific line
        :param line_index: the index of line in input
        :return: list, of cut indexes of line
        """
        if line_index in self._line_cuts:
            return sorted(self._line_cuts.get(line_index))
        else:
            return None

    def edge_line_segment(self, edge):
        """
        Get original road segment in point sequence which is
        bridged by graph edges.
        :param edge: a tuple of two node ids
        :return: the point sequence of raw road segment
        """
        edge = self.edge_key_nodes(edge[0], edge[1])
        assert edge in self._edges
        return self._edges[edge].line_segment


if __name__ == '__main__':
    pass
