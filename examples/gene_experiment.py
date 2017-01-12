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

    # Set ids for rps in trans dict
    for i, key in enumerate(numeric_translation._b_to_a):
        rps = []

        # Get all region paths mapping to key
        for interval in numeric_translation._b_to_a[key]:
            rps.extend(interval.region_paths)
            new_id = str(i) + "".join((name_translation._b_to_a[rp][0].region_paths[0] for rp in rps))

        new_dict[key] = new_id

    # Set ids for rps not in trans dict
    for n_id in name_translation._b_to_a:
        if n_id not in numeric_translation._a_to_b:
            new_dict[n_id] = name_translation._b_to_a[n_id][0].region_paths[0]

    a_to_b = new_dict
    trans = Translation.make_name_translation(a_to_b, graph)
    new_graph = trans.translate_subgraph(graph)
    trans.graph2 = new_graph
    return new_graph, trans


def merge_flanks(intervals, final_trans, new_graph, name_translation):
    """
    Merges the start and end flanks represented in intervals on the given graph
    :param intervals: List of start and end flanks to merge [start, start, end, end]
    Intervals should be on first graph (i.e. name_translation.graph1)
    :param final_trans: Trans that will be updated
    :param new_graph: Current graph
    :param name_translation: Translation from human readable names to numeric IDs. Can be an empty translation
    :return: Returns the new graph and translation as a tuple
    """
    #if final_trans is not None:
    #print("=== Flanking intervals ===")
    #print(intervals)
    print('\n'.join([str(i) for i in intervals]))
    # Merge start flank of alt locus with main
    merge_intervals = intervals[0:2]
    merge_intervals = [name_translation.translate(i) for i in merge_intervals]
    merge_intervals = [final_trans.translate(i) for i in merge_intervals]
    #print("=== Flanking intervals after translate===")
    #print('\n'.join([str(i) for i in merge_intervals]))
    #print("=== Graph ===")
    #print(merge_intervals[0].graph)
    for intv in merge_intervals:
        intv.graph = new_graph
    if merge_intervals[0].length() > 0:
        for interval in merge_intervals:
            print(interval)
        new_graph, trans = new_graph.merge(merge_intervals)
        final_trans += trans
    else:
        # Only connect by edge
        new_graph, trans = new_graph.connect_postitions(
            new_graph.prev_position(merge_intervals[0].start_position),
            merge_intervals[1].start_position
        )
        final_trans += trans
        final_trans.graph2 = new_graph.copy()
        final_trans.graph2._update_a_b_graph(final_trans._a_to_b, new_graph)

    # Merge end flank of alt locus with main

    merge_intervals = intervals[2:4]

    if merge_intervals[0].length() > 0:
        merge_intervals = [final_trans.translate(name_translation.translate(i))
                           for i in merge_intervals]
        for intv in merge_intervals:
            intv.graph = new_graph

        new_graph, trans = new_graph.merge(merge_intervals)
        final_trans += trans
    else:
        # Change position 1 back for alt loci
        ig = name_translation.graph1
        merge_intervals[1].start_position = \
            ig.prev_position(merge_intervals[1].start_position)

        merge_intervals = [final_trans.translate(name_translation.translate(i))
                           for i in merge_intervals]

        for intv in merge_intervals:
            intv.graph = new_graph

        # Only connect by edge
        new_graph, trans = new_graph.connect_postitions(
            merge_intervals[1].start_position,
            new_graph.next_position(merge_intervals[0].start_position)
        )

        #print("=== end flank trans ===")
        #print(trans)

        final_trans += trans
        final_trans.graph2 = new_graph.copy()
        final_trans.graph2._update_a_b_graph(final_trans._a_to_b, new_graph)
        #print("No end flank")

    return new_graph, final_trans


def connect_without_flanks(graph, alt_loci_fn, name_translation):
    """
    Connects the alternative loci in the given file to the grch38 graph,
    without flanks.
    :param alt_loci_fn: Filename of file containing alternative loci.
    One alt locus on each line. Four columns: alt_locus_id  chr chr_start   chr_stop
    :return: Returns the new graph
    """
    f = open(alt_loci_fn)
    n_flanks = 0
    new_graph = graph
    print("===== first new graph ====")
    print(new_graph)
    orig_graph = graph.copy()
    final_trans = Translation(graph=graph)
    final_trans.graph2 = graph
    for line in f.readlines():
        if line.startswith("#"):
            print("Skipping line")
            continue
        print("================= Iteration ==")
        print(line)
        l = line.split()
        alt_locus_id = l[0]
        main_chr = l[1]
        start = int(l[2])
        end = int(l[3])
        length = int(l[4])

        intervals = flanks.get_flanks(alt_locus_id, length,
                                      main_chr, start-1, end)
        new_graph, final_trans = merge_flanks(intervals, final_trans,
                                              new_graph, name_translation)

    f.close()
    print("NUMBER OF FLANKS: %d" % n_flanks)
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


def get_genes_as_intervals(fn):
    """
    Returns a dict. Keys are gene names and values are intervals representing the gene
    :param fn: File name of file containing genes (on the format of UCSC)
    :return:
    """
    genes = parse_genes_file(fn)
    out = {}
    for gene in genes:
        print(gene)

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
        #print(graph)
        numeric_graph, name_translation = convert_to_numeric_graph(graph)

        print("=== Numeric graph ===")
        #print(numeric_graph)
        #print(name_translation)

        print("=== Connecting ===")
        new_numeric_graph, numeric_translation = connect_without_flanks(
            numeric_graph, sys.argv[2], name_translation)

        name_graph, new_name_translation = convert_to_text_graph(
            new_numeric_graph, name_translation, numeric_translation)
        final_translation = name_translation + numeric_translation + new_name_translation
        print(final_translation)
        genes = get_genes_as_intervals(sys.argv[3]) # parse_genes_file(sys.argv[3])
        # print(genes)
    print("SUCCESS!")

