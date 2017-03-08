.. OffsetBasedGraph documentation master file, created by
   sphinx-quickstart on Tue Oct  4 14:29:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Python package *OffsetBasedGraph* is a simple implementation of offset-based sequence graphs.
The package lets you represent graphs as well as the *Translation* (difference) between them,
making it possible to e.g. convert intervals from one graph to another.

Installation
=================
The package with dependencies can be installed with pip::

    pip3 install offsetbasedgraph

If you want to use the example data to run the examples using real data, you need to clone the repository::

    git clone git@github.com:uio-cels/OffsetBasedGraph.git
    cd OffsetBasedGraph
    pip3 install -e

Test that installation works, e.g. by trying to import the package::

    python3
    >>> import offsetbasedgraph



Quickstart
=================
The following are some simple examples to get you started. Alternatively, you can jump directly to the the
documentation below.

Create a simple graph
^^^^^^^^^^^^^^^^^^^^^^

Graphs can be created by passing a dict of blocks and a dict of edges to the init method.
In the following example, we create a very simple graph consisting of two connected blocks,
and save that graph to file::

    from offsetbasedgraph import Block, Graph
    blocks = {"chr1-part1": Block(20), "chr1-part2": Block(10)}
    adjencies = {"chr1-part1": ["chr1-part2"]}
    graph = Graph(blocks, adjencies)
    graph.to_file("mygraph.graph")

The graph can later be read from file again::

    graph = Graph.from_file("mygraph.graph")

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

    gene = Interval(0, 8, ["chr1"], graph)
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
            "merged_part": [Interval(0, 5, ["chr1"]), Interval(0, 5, ["chr1_alt1"], graph)],
            "chr1-unmerged": [Interval(5, 20, ["chr1"], graph)],
            "chr1_alt1-unmerged": [Interval(5, 10, ["chr1_alt1"], graph)]
        },
        graph=graph
    )


Now, it is easy to translate the gene interval to our merged graph::

    translated_interval = merge_translation.translate(gene)
    print(translated_interval)

Intervals and positions
^^^^^^^^^^^^^^^^^^^^^^^^^^
Positions consist of a region path ID and an offset::

    from offsetbasedgraph  import Position
    my_position = Position("chr1_some_alt", 100)


Singlepath intervals consist of a start position, end position and a list of region paths ids::

    from offsetbasedgraph import Interval
    start_position = Position("chr1", 100)
    end_position = Position("chr1_some_alt", 100)
    interval = Interval(start_position, end_position, ["chr1", "chr1_some_alt"])

The interval init method can also be initialized using only the start and end offsets instead of full start and end positions::

    interval = Interval(100, 100, ["chr1", "chr1_some_alt"])

The package also supports some multipath intervals. See the documentation.


Example using GRCh38
=============================================
This example assumes that you have cloned the repository (see Installation) and that you are positioned in the
examples directory.

Included in the package are some utility functions that makes it
easy to create GRCh38-like graphs.
We here create a simple graph covering an alt locus on GRCh38::

    from offsetbasedgraph.graphcreators import create_initial_grch38_graph, convert_to_numeric_graph
    graph = create_initial_grch38_graph("chrom.sizes.example")
    numeric_graph, name_translation = convert_to_numeric_graph(graph)
    print(numeric_graph)

*graph*  is now a simple graph with two blocks (an alt locus and chr1), and
*numeric_graph* is the same graph converted so that all blocks have numeric IDs.
*name_translation* is a translation from the original block names to the numeric format.
This numeric format makes it easier to work with the graph, but
we can always get back to the original graph later by using the translation object.

We now connect
the alt locus to chr1 and remove flanking regions of the alt locus
(regions identical to the main chromosome)::

    from offsetbasedgraph.graphcreators import connect_without_flanks
    connected_graph, connected_trans = connect_without_flanks(numeric_graph, "alt.loci.example", name_translation)
    print(connected_graph)

*connected_trans* now represents all the operations done when connecting the
alt locus. Using this translation object, we can translate intervals from an original
GRC3h8 graph (*graph*) to our *connected_graph*::

    from offsetbasedgraph.gene import GeneList
    gene_list = GeneList.from_file("genes.example")
    translated_genes = gene_list.translate(connected_trans)

We have now successfully represented genes on a graph based on GRCh38.

Create full GRCh38 graph
=============================================
The following shows a short snippet for building a graph from GRCh38, with merged flanks.
It assumes you are positioned in the examples directory.

NB: This code takes time to run, as remote sequence data needs to be downloaded::

    graph = create_initial_grch38_graph("grch38.chrom.sizes")
    numeric_graph, name_translation = convert_to_numeric_graph(graph)
    new_numeric_graph, numeric_translation = connect_without_flanks(
        numeric_graph, "grch38_alt_loci.txt", name_translation)
    name_graph, new_name_translation = convert_to_text_graph(
        new_numeric_graph, name_translation, numeric_translation)
    final_translation = name_translation + numeric_translation + new_name_translation
    final_translation.graph2 = name_graph

    grch38_graph_with_flanks = name_graph

    print(grch38_graph_with_flanks.summary())
    
    # Translate example genes
    gene_list = GeneList.from_file("genes.example")
    print(gene_list)

    translated_gene_list = gene_list.translate(final_translation)
    print(translated_gene_list)



Documentation
============================

.. toctree::
   :maxdepth: 2


.. automodule:: offsetbasedgraph

Graph
^^^^^^^^^^^^^^^^
.. autoclass:: Graph
    :members:
    :inherited-members:
    :show-inheritance:

Translation
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: Translation
    :members:

Interval and Position
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: Interval
    :members:

.. autoclass:: Position
    :members:

Genes
^^^^^^^^^^^^
.. automodule:: offsetbasedgraph.gene

.. autoclass:: Gene
    :members:

Utility functions
^^^^^^^^^^^^^^^^^^^^^^^^

Graphutils
"""""""""""""

.. automodule:: offsetbasedgraph.graphutils
    :members:

Graphcreators
"""""""""""""""""

.. automodule:: offsetbasedgraph.graphcreators
    :members:
