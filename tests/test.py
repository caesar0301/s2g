#!/usr/bin/env python
import matplotlib
import os
import sys
import time
import unittest

import fiona
from shapely.geometry import shape, box, MultiLineString
import networkx as nx

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# add local s2g repo
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import s2g


class ShapeGraphCase(unittest.TestCase):

    def setUp(self):
        self.sg = ShapeGraphCase.sg
        self.show_plots = False

    @classmethod
    def setUpClass(cls):
        """ Creating ShapeGraph is slow, to avoid calling it for each test use setUpClass()
            and store the result as class variable
        """
        super(ShapeGraphCase, cls).setUpClass()
        # test shapefile input
        cls.shp = os.path.join(os.path.dirname(__file__), '../data/campus.shp')
        G1 = s2g.ShapeGraph(shapefile=cls.shp, to_graph=True, resolution=0.01,
                            properties=['osm_id'], geom_count=100)
        cls.sg = G1

    def test_coords_input(self):
        geoms = []
        with fiona.open(self.shp) as source:
            for r in source:
                s = shape(r['geometry'])
                geoms.append(s)
        s2g.ShapeGraph(geoms[0:30], to_graph=True, resolution=0.01,
                       properties=['osm_id'])

    def test_registered_edges(self):
        for edge, segment in self.sg._edges.items():
            assert edge[0] <= edge[1]
            if segment.line_index is None:
                assert segment.cuts is None
            else:
                s, e = segment.edge
                sc, ec = segment.cuts
                assert self.sg.node_xy[s] == self.sg.geoms[segment.line_index].coords[sc]

    def test_edge_properties(self):
        for edge, segment in self.sg._edges.items():
            if segment.line_index is not None:
                props = self.sg.line_props(segment.line_index)
                assert 'osm_id' in props

    def test_edge_key_nodes(self):
        assert self.sg.edge_key_nodes(2, 1) == (1, 2)
        assert self.sg.edge_key_nodes((2, 1)) == (1, 2)

    def test_edge_info(self):
        info = self.sg.edge_info((0, 1))
        assert info.line_index == 0
        assert info.cuts == (0, 1)

    def test_line_info(self):
        info = self.sg.line_info(0)
        assert 'osm_id' in info.props

    def test_major_lines_info(self):
        lines_info = self.sg.major_lines_info()
        major = self.sg.major_component()
        for line_index, info in lines_info.items():
            assert info.is_major is True
            assert info.index in major
        assert len(lines_info) == len(major)

    def test_line_edge_sequence(self):
        for line_index in self.sg.major_component():
            edges = self.sg.line_edge_sequence(line_index)
            for edge in edges:
                ekey = self.sg.edge_key_nodes(edge[0], edge[1])
                assert ekey in self.sg.graph.edges()

    def test_subgraph_within_box(self):
        bounding_box = box(121.428387, 31.027371, 121.430863, 31.030227)
        a = time.time()
        subgraph = self.sg.subgraph_within_box(bounding_box)
        print time.time() - a
        plt.figure()
        nx.draw(subgraph, pos=self.sg.node_xy, node_size=50)
        if self.show_plots:
            plt.show()

    def test_lines_within_box(self):
        bounding_box = box(121.428387, 31.027371, 121.430863, 31.030227)
        lines = self.sg.lines_within_box(bounding_box)
        for line in lines:
            assert line.is_major

    def test_point_projects_to_edge(self):
        # p = (114.83299055, 26.8892277)
        p = (121.428387, 31.027371)
        a = time.time()
        edges, segments = self.sg.point_projects_to_edges(p, 0.01)
        print time.time() - a

        plt.figure()
        s2g.plot_lines(MultiLineString(segments), color='orange')  # original roads
        for i in range(0, len(edges)):
            s, e = edges[i]
            sxy = self.sg.node_xy[s]
            exy = self.sg.node_xy[e]
            plt.plot([sxy[0], exy[0]], [sxy[1], exy[1]], color='green') # graph edges
        plt.plot(p[0], p[1], color='red', markersize=12, marker='o') # bridges
        if self.show_plots:
            plt.show()

    def test_edge_line_segment(self):
        pass

    def tearDown(self):
        pass


class BonusCase(unittest.TestCase):
    def setUp(self):
        shp = os.path.join(os.path.dirname(__file__), '../data/campus.shp')
        self.geoms = []
        with fiona.open(shp) as source:

            for r in source:
                s = shape(r['geometry'])
                self.geoms.append(s)

    def test_gcd(self):
        p1 = (114.83299055, 26.8892277)
        p2 = (121.428387, 31.027371)
        assert s2g.great_circle_dist(p1, p2) == s2g.gcd(p1, p2)
        assert s2g.great_circle_dist(p1, p1) == 0

    def test_bounded_segments(self):
        pass

    def tearDown(self):
        pass


def suite():
    suites = [ShapeGraphCase, BonusCase]
    suite = unittest.TestSuite()
    for s in suites:
        suite.addTest(unittest.makeSuite(s))
    return suite


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
