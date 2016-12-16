# python-s2g

(S)hapefile (2) Graph/network converter in Python

[![Build Status](https://travis-ci.org/caesar0301/python-s2g.svg?branch=master)](https://travis-ci.org/caesar0301/python-s2g)

When we process GIS data, a non-trivial problem is the conversion from shape lines to graph or network data structure.
The latter may benefit from these out-of-box graphical libraries such as [networkx](http://networkx.github.io/)
and [igraph](http://igraph.org/python/). But the conversion is a headache to components open communities.
This mostly urges me to finish this tiny but useful library.
 
# Install

```$xslt
pip install -U s2g
```

# Usage

You have two alternative ways to construct the graph. One is reading from a raw shapefiles with `LineString` objects.
(Under the hood, I involve [fiona](https://pypi.python.org/pypi/Fiona/) to read geometries and
[shapely](https://pypi.python.org/pypi/Shapely) to analyze the data.).
Currently, this tool only supports conversion to *undirected graph*.

```python
from s2g import ShapeGraph
import networkx as nx

sg = ShapeGraph(shapefile='path/to/roads.shp', to_graph=True)
assert isinstance(sg.graph, nx.Graph)
```

The other way is designed for programmable usage or time-consuming process where intermediate data could be sniffed or
saved. Here is an example to read lines with [fiona]:

```python
from s2g import ShapeGraph
import fiona
from shapely.geometry import shape, LineString

shp = 'path/to/shapefile.shp'

with fiona.open(shp) as source:
    geoms = []
    for r in source:
        s = shape(r['geometry'])
        if isinstance(s, LineString):
            geoms.append(s)

# create ShapeGraph object from a list of lines
sg = ShapeGraph(geoms, to_graph=False)

# detect major components
mc = sg.gen_major_components()
# major components are mc[2]

# convert the largest component to networkx Graph
graph = sg.to_networkx()  # equivalently sg.graph
```

Dive into [source doc](https://github.com/caesar0301/python-s2g/blob/master/s2g/shapegraph.py) to discover other functionalities.

## References

* [shp2graph](https://cran.r-project.org/web/packages/shp2graph/index.html) in R by Binbin Lu, as well as [his talk](http://web.warwick.ac.uk/statsdept/user2011/TalkSlides/Contributed/17Aug_1600_FocusIV_2-Geospatial_1-Lu.pdf) on useR! 2011
* [A Tutorial on Topology Correction of Shapefiles Using GRASS](http://xiaming.me/posts/2015/08/29/a-tutorial-on-topology-correction-of-shapefiles-using-grass/)
