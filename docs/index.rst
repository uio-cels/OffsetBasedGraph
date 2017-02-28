.. OffsetBasedGraph documentation master file, created by
   sphinx-quickstart on Tue Oct  4 14:29:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Package OffsetBasedGraph
============================================

The Python package *OffsetBasedGraphs* is a simple implementation of offset-based sequence graphs.
The package lets you represent graphs as well of the *Translation* (difference) between them,
making it possible to e.g. convert interval from one graph to another.


Quickstart - how to create simple graphs
============================================

The only way to create a Graph is through the init method:

    >>> blocks = {"chr1-part1": Block(20), "chr1-part2": Block(10)}
    >>> adjencies = {"chr1-part1": ["chr1-part2"]}
    >>> graph = Graph(blocks, adjencies)
    >>> print(graph)

Graphs can be written to and read from files:

    >>> graph.to_file("mygraph.graph")
    >>> mygraph = Graph.from_file("mygraph.graph")

Offset-based graphs are more interesting when
one can represent intervals on them, and "translate"
intervals between different graphs.

Contents:

.. toctree::
   :maxdepth: 2

.. automodule:: offsetbasedgraph
 
.. autoclass:: Graph
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:

.. autoclass:: Interval
    :members:

.. autoclass:: Position
    :members:

.. autoclass:: AltLocus
    :members:

.. autoclass:: AltLoci
    :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

