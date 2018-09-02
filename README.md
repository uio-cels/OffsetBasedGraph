[![Build Status](https://travis-ci.org/uio-cels/OffsetBasedGraph.svg?branch=v2.0)](https://travis-ci.org/uio-cels/OffsetBasedGraph)

**Note**: This is Offsetbasedgraph version 2.1, made to work with [Graph Peak Caller](https://github.com/uio-bmi/graph_peak_caller). For the version compatible with the analysis done in the paper *Coordinates and Intervals in Graph-based Reference Genomes*, go to [version *1.0.7*](https://github.com/uio-cels/OffsetBasedGraph/tree/release2).

# Offsetbasedgraph 
This is a python package for working with and representing genomic graphs without sequence in Python. Together with [pyvg](https://github.com/uio-bmi/pyvg), Offsetbasedgraph can be used to represent graphs created by [vg](https://github.com/vgteam/vg/).

# Installation
Offsetbasedgraph version 2.0 is compatible with Python 3 only. Install with pip:
```
pip3 install offsetbasedgraph
```
... or by cloning:
```
git clone https://github.com/uio-cels/OffsetBasedGraph.git
cd OffsetBasedGraph
pip3 install -e .
```

# Basic usage
Graphs can be created by sending a dict of Blocks (connected nodes) and dict of edges (one list of edges per key node) to the Graph object. 
```
from offsetbasedgraph import Graph, Block
graph = Graph({1: Block(10), 2: Block(5)}, {1: [2]})
```

An interval contain a start offset, end offset and a list of blocks it follows:
```
from offsetbasedgraph import Interval
interval = Interval(3, 6, [1, 2], graph)
```

Offsetbasedgraph can be used together with [pyvg](https://github.com/uio-bmi/pyvg) to represent huge graphs created by vg (requires vg json graph format):
```
from pyvg.conversion import json_file_to_obg_numpy_graph
offset_based_graph = json_file_to_obg_numpy_graph("vggraph.json")
```

