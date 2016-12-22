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

    @classmethod
    def setUpClass(cls):
        """ Creating ShapeGraph is slow, to avoid calling it for each test use setUpClass()
            and store the result as class variable
        """
        super(ShapeGraphCase, cls).setUpClass()
        # test shapefile input
        shp = os.path.join(os.path.dirname(__file__), '../data/campus.shp')
        G1 = s2g.ShapeGraph(shapefile=shp, to_graph=True, resolution=0.01,
                            properties=['osm_id'], geom_count=50)
        # test line or coordinates list input
        geoms = []
        with fiona.open(shp) as source:
            for r in source:
                s = shape(r['geometry'])
                geoms.append(s)
        G2 = s2g.ShapeGraph(geoms[0:50], to_graph=True, resolution=0.01,
                            properties=['osm_id'])
        cls.sg = G1

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
                props = self.sg.get_line_props(segment.line_index)
                assert 'osm_id' in props

    def test_edge_info(self):
        info = self.sg.edge_info((0, 1))
        assert info.line_index == 0
        assert info.cuts == (0, 1)

    def test_line_edge_sequence(self):
        for line_index in self.sg.largest_component():
            edges = self.sg.line_edge_sequence(line_index)
            for edge in edges:
                assert self.sg.edge_key_nodes(edge[0], edge[1]) in self.sg.graph.edges()

    def test_subgraph_within_box(self):
        bounding_box = box(114.572, 26.769, 114.971, 26.933)
        a = time.time()
        subgraph = self.sg.subgraph_within_box(bounding_box)
        print time.time() - a
        nx.draw(subgraph, pos=self.sg.node_xy)
        # plt.show()

    def test_point_projects_to_edge(self):
        p = (114.83299055, 26.8892277)
        a = time.time()
        edges, segments = self.sg.point_projects_to_edges(p, 1)
        print time.time() - a
        s2g.plot_lines(MultiLineString(segments), color='orange')  # original roads
        for i in range(0, len(edges)):
            s, e = edges[i]
            sxy = self.sg.node_xy[s]
            exy = self.sg.node_xy[e]
            plt.plot([sxy[0], exy[0]], [sxy[1], exy[1]], color='green') # graph edges
        plt.plot(p[0], p[1], color='red', markersize=12, marker='o') # bridges
        # plt.show()

    def test_road_segment_from_edge(self):
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

    def test_lines_within_box(self):
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
