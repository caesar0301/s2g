#!/usr/env/bin python
# encoding: utf-8
import logging
import matplotlib
import sys
import pickle

import networkx as nx
import numpy as np

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from shapely.geometry import shape, Point, LineString
from shapely.ops import linemerge
from itertools import product

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

__all__ = [
    'plot_lines', 'great_circle_dist', 'line_distance', 'line_contains',
    'lines_touch', 'point_projects_to_line', 'point_projects_to_points',
    'cut_line', 'ShapeGraph'
]


def plot_lines(lines):
    def plot_line(ob):
        x, y = ob.xy
        plt.plot(x, y, linewidth=1, solid_capstyle='round', zorder=1)

    for u in lines:
        if u.geom_type in ['LineString', 'LinearRing', 'Point']:
            plot_line(u)
        elif u.geom_type is ['MultiLineString']:
            for p in u:
                plot_line(p)


def great_circle_dist(p1, p2):
    """Return the distance (in km) between two points in
    geographical coordinates.
    """
    lon0, lat0 = p1
    lon1, lat1 = p2
    EARTH_R = 6372.8
    lat0 = np.radians(float(lat0))
    lon0 = np.radians(float(lon0))
    lat1 = np.radians(float(lat1))
    lon1 = np.radians(float(lon1))
    dlon = lon0 - lon1
    y = np.sqrt(
        (np.cos(lat1) * np.sin(dlon)) ** 2
        + (np.cos(lat0) * np.sin(lat1)
           - np.sin(lat0) * np.cos(lat1) * np.cos(dlon)) ** 2)
    x = np.sin(lat0) * np.sin(lat1) + \
        np.cos(lat0) * np.cos(lat1) * np.cos(dlon)
    c = np.arctan2(y, x)
    return EARTH_R * c


def line_distance(coords):
    """Return total road distance in kilometers"""
    dist = []
    for i in range(0, len(coords) - 1):
        dist.append(great_circle_dist(coords[i], coords[i + 1]))
    return np.sum(dist)


def line_contains(line, point, buf=10e-5):
    p = Point(point).buffer(buf)
    return line.intersects(p)


def lines_touch(one, other, buf=10e-5):
    """Predict the connection of two lines"""

    def touches_ends(one, other, buf=10e-5):
        v1 = Point(other.coords[0]).buffer(buf)
        v2 = Point(other.coords[-1]).buffer(buf)
        if one.intersects(v1) or one.intersects(v2):
            return True
        return False

    return touches_ends(one, other, buf) or \
           touches_ends(other, one, buf)


def point_projects_to_line(point, line):
    """Get the nearest point index on line
    """
    nearest = None
    min_dist = -1
    for i in range(0, len(line.coords)):
        p = line.coords[i]
        # d = great_circle_dist(p, point)
        d = Point(p).distance(Point(point))
        if min_dist < 0:
            min_dist = d
            nearest = i
        else:
            if d < min_dist:
                min_dist = d
                nearest = i
    return nearest


def point_projects_to_points(point, others):
    """Get the nearest point given a group of points
    """
    p = Point(point)
    nearest = None
    min_dist = -1
    for other in others:
        dist = Point(other).distance(p)
        if min_dist < 0:
            min_dist = dist
            nearest = other
        else:
            if dist < min_dist:
                min_dist = dist
                nearest = other
    return nearest


def cut_line(line, resolution=1.0):
    assert line.geom_type == 'LineString'
    coords = line.coords
    sampled_points = [0]
    distances = [0]
    prev = coords[0]
    acc_dist = 0
    added = False
    for i in range(1, len(coords)):
        acc_dist += great_circle_dist(coords[i - 1], coords[i])
        if acc_dist >= resolution:
            added = True
            sampled_points.append(i)
            prev = coords[i]
            distances.append(acc_dist)
            acc_dist = 0
        else:
            added = False
    if not added:
        sampled_points.append(i)
    distances.append(acc_dist)
    return sampled_points, distances


class ShapeGraph(object):
    def __init__(self, geoms):
        self.geoms = [line for line in geoms]  # raw lines
        self.connected = False
        self.connectivity = None
        self.major_components = None
        # global edge info
        self.node_ids = {}
        self.node_xy = {}
        self._edges = {}
        # cut-line info
        self.line_cuts = {}
        self.nodes_counter = 0
        self.graph = nx.Graph()

    def _register_edge(self, p1, p2, dist):
        assert isinstance(p1, Point) or len(p1) == 2
        if p1 not in self.node_ids:
            self.node_ids[p1] = self.nodes_counter
            self.nodes_counter += 1
            self.node_xy[self.node_ids[p1]] = p1
        if p2 not in self.node_ids:
            self.node_ids[p2] = self.nodes_counter
            self.nodes_counter += 1
            self.node_xy[self.node_ids[p2]] = p2
        edge = tuple(sorted([self.node_ids[p1], self.node_ids[p2]]))
        if edge not in self._edges:
            self._edges[edge] = dist

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
        edge = tuple(sorted([id1, id2]))
        if edge in self._edges:
            del self._edges[edge]

    def gen_major_components(self):
        """Merge individual lines into groups of touching lines.
        """
        lines = self.geoms
        L = len(lines)
        pairs = []
        for i, j in product(range(0, L), range(0, L)):
            if j < i:
                continue
            if lines_touch(lines[i], lines[j]):
                pairs.append((i, j))

        # validate lines connectivity
        graph = nx.Graph()
        [graph.add_edge(i[0], i[1]) for i in pairs]

        # Pairs of lines which are connected (touched). These line pairs
        # are not guaranteed to be in the same major component.
        self.connectivity = [i for i in graph.edges() if i[0] != i[1]]

        # major connected components
        cc = nx.algorithms.components.connected_components(graph)
        self.major_components = sorted([i for i in cc], key=len, reverse=True)

        # predict if all lines are strongly connected
        self.connected = True if len(self.major_components) == 1 else False
        return \
            self.connected, self.connectivity, self.major_components

    def dump_major_components(self, ofile):
        """Dump lines' connectivity and major components to pickle file.
        """
        pickle.dump((self.connected, self.connectivity, self.major_components), ofile)

    def load_major_components(self, pickle_file):
        """Load lines' connectivity and major components from pickle file.
        """
        self.connected, self.connectivity, major_components = \
            pickle.load(pickle_file)
        self.major_components = sorted(major_components, key=len, reverse=True)

    def to_networkx(self, resolution=1.0, coord_type='lonlat', point_buffer=10e-5):
        """Convert the major component to graph of NetworkX.
        """
        if self.major_components is None:
            self.gen_major_components()

        major = list(self.major_components[0])
        logging.info('Use the largest component of {0} lines'.format(len(major)))

        # do line cutting
        logging.info('Cutting lines with specific resolution = {0} km'.format(resolution))
        for lid in major:
            line = self.geoms[lid]
            cuts, dist = cut_line(line, resolution)

            # record cut-line info
            self.line_cuts[lid] = (cuts, dist)

            # record edges info
            for j in range(1, len(cuts)):
                self._register_edge(line.coords[cuts[j - 1]],
                                    line.coords[cuts[j]],
                                    dist[j])

        # intersection interpolation for each connected edge pair
        logging.info('Intersection interpolation for connected (touched) edges')
        for s, e in self.connectivity:
            if s == e:
                continue
            if s not in major or e not in major:
                continue

            sline = self.geoms[s]
            eline = self.geoms[e]
            s1, s2 = sline.coords[0], sline.coords[-1]
            e1, e2 = eline.coords[0], eline.coords[-1]

            for l, p in zip([s, s, e, e], [e1, e2, s1, s2]):
                self._add_point_to_line(l, p)
                assert p in self.node_ids

        # assemble graph
        for edge, dist in self._edges.items():
            self.graph.add_edge(edge[0], edge[1], weight=dist)

        logging.info('Graph created with {0} nodes, {1} edges'.format(
            len(self.graph.edges()), len(self.graph.nodes())))

        return self.graph

    def _add_point_to_line(self, line_index, point, point_buffer=10e-5):
        """Project a point to line when they intersect within a buffer.
        """
        coords = self.geoms[line_index].coords
        cuts, distances = self.line_cuts[line_index]

        found = False
        point_result = None
        for i in range(1, len(cuts)):
            sp = cuts[i - 1]
            ep = cuts[i]

            segment = LineString(coords[sp:ep + 1])
            assert coords[ep] == segment.coords[-1]
            buffered_point = Point(point).buffer(point_buffer)

            if segment.intersects(buffered_point):
                found = True
                offset = point_projects_to_line(point, segment)
                point_result = sp + offset
                # assert segment.coords[offset] == coords[sp + offset]

                logging.debug('Insert new point: {0} at index {1} of line[{2}]' \
                              .format(point, sp + offset, line_index))

                # split original segment into two parts
                self._remove_edge(segment.coords[0], segment.coords[-1])
                # replace point with projected one in distance calculation
                self._register_edge(segment.coords[0], point,
                                    line_distance(segment.coords[0:offset + 1]))
                self._register_edge(point, segment.coords[-1],
                                    line_distance(segment.coords[offset:]))
            if found:
                break
        return point_result

    def show_graph(self):
        nx.draw_networkx(self.graph, pos=self.node_xy, node_size=10,
                         with_labels=False, arrows=False)
        plt.show()

    def point_projects_to_node(self, point, distance_tolerance=0.1):
        """
        Return the nearest graph node to given point.
        :param point: a tuple of (lon, lat) or shapely Point
        :param distance_tolerance: approximated buffer range in km
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

    def point_projects_to_edges(self, point, distance_tolerance=0.0):
        pass


if __name__ == '__main__':
    pass