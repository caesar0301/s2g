#!/usr/bin/env python
import sys
import os
import fiona

from shapely.geometry import shape

# add local s2g repo
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import s2g


if __name__ == '__main__':
    shp = os.path.join(os.path.dirname(__file__), '../data/city.shp')

    with fiona.open(shp) as source:
        geoms = []
        for r in source:
            s = shape(r['geometry'])
            geoms.append(s)
    sg = s2g.ShapeGraph(geoms, to_graph=False)

    sg = s2g.ShapeGraph(shapefile=shp, to_graph=True)
