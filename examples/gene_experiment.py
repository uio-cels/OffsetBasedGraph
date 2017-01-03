"""
Creates a graph from GRCh38, represents genes on this graph and outputs some information about the genes

Usage example:
$ python3 gene_experiment.py grch38.chrom.sizes grch38_alt_loci.txt genes_chr1_GL383518v1_alt.txt

With only two alt loci:
python3 gene_experiment.py grch38.chrom.sizes-small grch38_alt_loci_small.txt genes_chr1_GL383518v1_alt.txt

"""

import sys
import csv
from offsetbasedgraph import Graph, Block, Translation
from genutils import flanks


def convert_to_numeric_graph(graph):
    a_to_b = {k: i for i, k in enumerate(graph.blocks.keys())}
    trans = Translation.make_name_translation(a_to_b, graph)
    new_graph = trans.translate_subgraph(graph)
    trans.graph2 = new_graph

    # ?:
    for intervals in trans._a_to_b.values():
        for interval in intervals:
            interval.graph = new_graph

    return new_graph, trans


def create_initial_grch38_graph(chrom_sizes_fn):
    """
    Creates an initial grch38 graph with no connected blocks
    :param chrom_sizes_fn: Filename of chrom sizes file
    :return: Returns a Graph with no connected blocks
    """
    blocks = {}
    with open(chrom_sizes_fn, 'r') as csvfile:
        chroms = csv.reader(csvfile, delimiter='\t')
        for i, chrom in enumerate(chroms):
            blocks[chrom[0]] = Block(int(chrom[1]))

    return Graph(blocks, {})


def convert_to_text_graph(graph, name_translation, numeric_translation):
    for key in numeric_translation._b_to_a.keys():
        assert key in graph.blocks, "%s not in %s" % (key, graph)
    new_dict = {}

    print("->", numeric_translation._a_to_b)
    print("<-", numeric_translation._b_to_a)
    # Set ids for rps in trans dict
    for i, key in enumerate(numeric_translation._b_to_a):
        rps = []

        # Get all region paths mapping to key
        for interval in numeric_translation._b_to_a[key]:
            rps.extend(interval.region_paths)

        new_id = str(i) + "".join(rps)

        new_dict[key] = new_id

    # Set ids for rps not in trans dict
    for n_id in name_translation._b_to_a:
        if n_id not in numeric_translation._a_to_b:
            new_dict[n_id] = name_translation._b_to_a[n_id][0].region_paths[0]

    a_to_b = new_dict
    print("-->", a_to_b)
    # b_to_a = {v[0]: [k] for k, v in a_to_b.items()}
    trans = Translation.make_name_translation(a_to_b, graph)
    new_graph = trans.translate_subgraph(graph)
    trans.graph2 = new_graph
    return new_graph, trans


def connect_without_flanks(graph, alt_loci_fn, name_translation):
    """
    Connects the alternative loci in the given file to the grch38 graph,
    without flanks.
    :param alt_loci_fn: Filename of file containing alternative loci.
    One alt locus on each line. Four columns: alt_locus_id  chr chr_start   chr_stop
    :return: Returns the new graph
    """
    f = open(alt_loci_fn)
    new_graph = graph
    orig_graph = graph.copy()
    final_trans = name_translation
    print("=== Final trans graph1  ===")
    print(final_trans.graph1)
    for line in f.readlines():
        print("== Iteration ==")
        print(line)
        l = line.split()
        alt_locus_id = l[0]
        main_chr = l[1]
        start = int(l[2])
        end = int(l[3])
        length = int(l[4])

        intervals = flanks.get_flanks(alt_locus_id, length, main_chr, start-1, end)
        print(intervals)
        #if final_trans is not None:

        # Merge start flank of alt locus with main
        merge_intervals = intervals[0:2]
        merge_intervals = [final_trans.translate(i) for i in merge_intervals]
        for intv in merge_intervals:
            intv.graph = new_graph
        prev_graph = new_graph
        new_graph, trans = new_graph.merge(merge_intervals)

        final_trans += trans

        #final_trans.graph2 = trans.graph2

        # Merge end flank of alt locus with main

        merge_intervals = intervals[2:4]
        merge_intervals = [final_trans.translate(i) for i in merge_intervals]
        for intv in merge_intervals:
            intv.graph = new_graph

        prev_graph = new_graph
        new_graph, trans = new_graph.merge(merge_intervals)

        assert final_trans.graph1 == list(final_trans._b_to_a.values())[-1][0].graph
        assert final_trans.graph2 == list(final_trans._a_to_b.values())[-1][0].graph
        assert trans.graph2 == list(trans._a_to_b.values())[-1][0].graph
        assert trans.graph1 == list(trans._b_to_a.values())[-1][0].graph

        final_trans += trans



    return new_graph, final_trans


def parse_genes_file(genes_fn):
    """
    Parses a file containing genes (on the format of UCSC), and returns a list of dicts
    :param genes_fn: File name
    :return: Returns a list of genes, one dict for each gene
    """
    genes = []
    with open(genes_fn) as f:
        g = csv.DictReader(f, delimiter='\t')
        for gene in g:
            genes.append(gene)

    return genes


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("""
            Usage: gene_experiment.py chrom_sizes_file alt_loci_file gene_file \n
            Example: python3 gene_experiment.py grch38.chrom.sizes grch38_alt_loci.txt genes_chr1_GL383518v1_alt.txt\n
            Smaller example: python3 gene_experiment.py grch38.chrom.sizes-small grch38_alt_loci_small.txt genes_chr1_GL383518v1_alt.txt
            """)
    else:
        graph = create_initial_grch38_graph(sys.argv[1])

        print("=== First graph===")
        print(graph)
        numeric_graph, name_translation = convert_to_numeric_graph(graph)

        print("=== Numeric graph ===")
        print(numeric_graph)
        print(name_translation)

        print("=== Connecting ===")
        new_numeric_graph, numeric_translation = connect_without_flanks(
            numeric_graph, sys.argv[2], name_translation)


        name_graph, new_name_translation = convert_to_text_graph(
            new_numeric_graph, name_translation, numeric_translation)
        final_translation = name_translation + numeric_translation + new_name_translation
        genes = parse_genes_file(sys.argv[3])
        print(genes)


