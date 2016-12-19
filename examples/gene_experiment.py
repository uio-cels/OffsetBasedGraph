"""
Creates a graph from GRCh38, represents genes on this graph and outputs some information about the genes

Usage example:
$ python3 gene_experiment.py grch38.chrom.sizes grch38_alt_loci.txt genes_chr1_GL383518v1_alt.txt

"""

import sys
import csv
from offsetbasedgraph import Graph, Block
from genutils import flanks

def create_initial_grch38_graph(chrom_sizes_fn):
    """
    Creates an initial grch38 graph with no connected blocks
    :param chrom_sizes_fn: Filename of chrom sizes file
    :return: Returns a Graph with no connected blocks
    """
    blocks = {}
    with open(chrom_sizes_fn, 'r') as csvfile:
        chroms = csv.reader(csvfile, delimiter='\t')
        for chrom in chroms:
            blocks[chrom[0]] = int(chrom[1])
    return Graph(blocks, {})


def connect_without_flanks(graph, alt_loci_fn):
    """
    Connects the alternative loci in the given file to the grch38 graph,
    without flanks.
    :param alt_loci_fn: Filename of file containing alternative loci.
    One alt locus on each line. Four columns: alt_locus_id  chr chr_start   chr_stop
    :return: Returns the new graph
    """
    f = open(alt_loci_fn)
    new_graph = graph
    final_trans = None
    for line in f.readlines():
        l = line.split()
        alt_locus_id = l[0]
        main_chr = l[1]
        start = int(l[2])
        end = int(l[3])

        intervals = flanks.get_flanks(alt_locus_id, main_chr, start, end)
        if final_trans is not None:
            intervals = [final_trans.translate_interval(i) for i in intervals]

        new_graph, trans = new_graph.merge(intervals[0:2])
        if final_trans is None:
            final_trans = trans
        else:
            final_trans += trans

        new_graph, trans = new_graph.merge(intervals[2:4])
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
            Example: python3 gene_experiment.py grch38.chrom.sizes grch38_alt_loci.txt genes_chr1_GL383518v1_alt.txt
            """)
    else:
        graph = create_initial_grch38_graph(sys.argv[1])
        graph = connect_without_flanks(graph, sys.argv[2])
        genes = parse_genes_file(sys.argv[3])
        print(genes)


