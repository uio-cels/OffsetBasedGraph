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
--------------------------------------------

The only way to create a Graph is through the init method.
The following creates a very simple graph consisting of two connected blocks,
and saves that graph to file::

    from offsetbasedgraph import Block, Graph
    blocks = {"chr1-part1": Block(20), "chr1-part2": Block(10)}
    adjencies = {"chr1-part1": ["chr1-part2"]}
    graph = Graph(blocks, adjencies)
    graph.to_file("mygraph.graph")

The graph can later be read from file again::

    graph.to_file("mygraph.graph")

Translating between graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Offset-based graphs are more interesting when
one can represent intervals on them, and "translate"
intervals between different graphs.

Assume we have a simple graph with an alternative locus
and chrosome 1 as two blocks, not connected, and that we have a gene
on chromosome 1 (represented as a simple Interval)::

    from offsetbasedgraph import Block, Graph, Interval, Translation
    graph = Graph(
        {
            "chr1": Block(20),
            "chr1_alt1": Block(10),
        },
        {}
    )

    gene = Interval(0, 5, "chr1", graph)
    print(gene)

Now, assume we have a new graph where the first 5 base pairs of chr1 have
been merged with the first 5 base pairs of chr1_alt1. This operation
can easily be expressed by a *Translation* (see Translation class for details)::

    merge_translation = Translation(
        {
            "chr1_alt1": [Interval(0, 5, ["merged_part", "chr1_alt1-unmerged"])],
            "chr1": [Interval(0, 15, ["merged_part", "chr1-unmerged"])]
        },
        {
            "chr1_merged": [Interval(0, 5, ["chr1"]), Interval(0, 5, ["chr1_alt1"])],
            "chr1-unmerged": [Interval(5, 20, ["chr1"])],
            "chr1_alt1-unmerged": [Interval(5, 10, ["chr1_alt1"])]
        }
    )


Now, it is easy to translate the gene interval to our merged graph::

    translated_interval = merge_translation.translate(gene)
    print(translated_interval.get_single_path_intervals())

Example with real data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^




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

