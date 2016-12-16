#!/usr/bin/env python
import sys
import os
import fiona

from shapely.geometry import shape
import networkx as nx
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# add local s2g repo
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from s2g import ShapeGraph


if __name__ == '__main__':
    shp = os.path.join(os.path.dirname(__file__), '../data/campus.shp')

    with fiona.open(shp) as source:
        geoms = []
        for r in source:
            s = shape(r['geometry'])
            geoms.append(s)
    sg = ShapeGraph(geoms, to_graph=False)

    sg = ShapeGraph(shapefile=shp, to_graph=True, resolution=0.01)
    #nx.draw(sg.graph, pos=sg.node_xy, node_size=10, node_shape='.')
    #plt.show()
