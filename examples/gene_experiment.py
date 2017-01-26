"""
Creates a graph from GRCh38, represents genes on this graph and outputs some information about the genes

Usage example:
$ python3 gene_experiment.py grch38.chrom.sizes grch38_alt_loci.txt genes_chr1_GL383518v1_alt.txt

With only two alt loci:
python3 gene_experiment.py grch38.chrom.sizes-small grch38_alt_loci_small.txt genes_chr1_GL383518v1_alt.txt

"""

import sys
import csv
from offsetbasedgraph import Graph, Block, Translation, Interval, Position

from offsetbasedgraph.graphutils import Gene, convert_to_numeric_graph, connect_without_flanks, \
    convert_to_text_graph, merge_flanks, connect_without_flanks, parse_genes_file, \
    get_genes_as_intervals, get_gene_objects_as_intervals, find_exon_duplicates, \
    create_initial_grch38_graph

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("""
            Usage: gene_experiment.py chrom_sizes_file alt_loci_file gene_file \n
            Example: python3 gene_experiment.py grch38.chrom.sizes grch38_alt_loci.txt genes_chr1_GL383518v1_alt.txt\n
            Smaller example: python3 gene_experiment.py grch38.chrom.sizes-small grch38_alt_loci_small.txt genes_chr1_GL383518v1_alt.txt
            """)
    else:
#
#
#        graph = create_initial_grch38_graph(sys.argv[1])
#        numeric_graph, name_translation = convert_to_numeric_graph(graph)
#        new_numeric_graph, numeric_translation = connect_without_flanks(
#            numeric_graph, sys.argv[2], name_translation)
#        name_graph, new_name_translation = convert_to_text_graph(
#            new_numeric_graph, name_translation, numeric_translation)
#
#        # To save to disk:
#        graph.to_file("initial_grch38_graph")
#        name_translation.to_file("name_translation")
#        numeric_graph.to_file("numeric_graph")
#        new_numeric_graph.to_file("new_numeric_graph")
#        numeric_translation.to_file("numeric_translation")
#        new_name_translation.to_file("new_name_translation")
#
        # To read from disk:

        graph = Graph.from_file("initial_grch38_graph")
        numeric_graph = Graph.from_file("numeric_graph")
        name_translation = Translation.from_file("name_translation")
        new_numeric_graph = Graph.from_file("new_numeric_graph")
        numeric_translation = Translation.from_file("numeric_translation")
        new_name_translation = Translation.from_file("new_name_translation")

        #print(numeric_translation.graph1)


        final_translation = name_translation + numeric_translation
        final_translation = final_translation + new_name_translation
        genes = get_gene_objects_as_intervals(sys.argv[3], graph)
        find_exon_duplicates(genes, final_translation)
        exit()

        genes = get_genes_as_intervals(sys.argv[3], graph)
        genes_compact_graph = {}
        for g in genes:
            genes_compact_graph[g] = \
                numeric_translation.translate(
                    name_translation.translate(genes[g])
                )

        # To human readable graph
        genes_readable = {}
        for g in genes_compact_graph:
            genes_readable[g] = new_name_translation.translate(
                genes_compact_graph[g])

        print("\nOriginal genes")
        for g in genes:
            print("%s: %s" % (g, genes[g]))

        print("\nTranslated: ")
        for g in genes_compact_graph:
            print("%s: %s" % (g, genes_compact_graph[g]))

        print("\nReadable")
        for g in genes_readable:
            print("%s: %s" % (g, genes_readable[g]))

    print("SUCCESS!")

