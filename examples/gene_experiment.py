"""
Creates a graph from GRCh38, represents genes on this graph and outputs some information about the genes

Usage example:
$ python3 gene_experiment.py grch38.chrom.sizes grch38_alt_loci.txt genes_chr1_GL383518v1_alt.txt

With only two alt loci:
python3 gene_experiment.py grch38.chrom.sizes-small grch38_alt_loci_small.txt genes_chr1_GL383518v1_alt.txt

"""

import sys
import argparse
import csv
from offsetbasedgraph import Graph, Block, Translation, Interval, Position

from offsetbasedgraph.graphutils import Gene, convert_to_numeric_graph, connect_without_flanks, \
    convert_to_text_graph, merge_flanks, connect_without_flanks, parse_genes_file, \
    get_genes_as_intervals, get_gene_objects_as_intervals, find_exon_duplicates, \
    create_initial_grch38_graph, blast_test


def create_graph(args):
    graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    numeric_graph, name_translation = convert_to_numeric_graph(graph)
    new_numeric_graph, numeric_translation = connect_without_flanks(
        numeric_graph, args.alt_locations_file_name, name_translation)
    name_graph, new_name_translation = convert_to_text_graph(
        new_numeric_graph, name_translation, numeric_translation)

    final_translation = name_translation + numeric_translation + new_name_translation
    final_translation.to_file(args.out_file_name)
    print("Graph and translation object stored in %s" % (args.out_file_name))


def check_duplicate_genes(args):
    genes_file_name = args.genes_file_name
    final_trans = Translation.from_file(args.translation_file_name)
    genes = get_gene_objects_as_intervals(genes_file_name, final_trans.graph1)
    find_exon_duplicates(genes, final_trans)
    # print(genes_file_name)


def merge_alignment(args):
    trans = Translation.from_file(args.translation_file_name)
    graph = trans.graph1
    from offsetbasedgraph.graphutils import merge_alt_using_cigar, grch38_graph_to_numeric
    graph, trans = grch38_graph_to_numeric(graph)
    merge_alt_using_cigar(graph, trans, args.alt_locus_id)


def merge_all_alignments(args):
    from offsetbasedgraph.graphutils import merge_alt_using_cigar, grch38_graph_to_numeric
    text_graph = create_initial_grch38_graph(args.chrom_sizes_file_name)  # Text ids (chrom names and alt names)
    graph, numeric_trans = grch38_graph_to_numeric(text_graph)

    # Go through all alts in this graph
    new_graph = graph.copy()
    trans = numeric_trans
    #trans = Translation({}, {}, graph=graph)
    i = 0
    for b in text_graph.blocks:
        if "alt" in b:
            print("Merging %s" % b)
            trans, new_graph = merge_alt_using_cigar(new_graph, trans, b)
            i += 1
            if i > 2000:
                break

    print("To file")
    trans.to_file(args.out_file_name)

def visualize_genes(args):
    trans = Translation.from_file(args.translation_file_name)
    graph = trans.graph2

    # Find blocks that genes cover, create subgraph using them
    from offsetbasedgraph.graphutils import GeneList
    genes = GeneList.from_file(args.genes_file_name).gene_list
    blocks = []

    """
    for gene in genes:
        blocks.extend(gene.transcription_region.region_paths)
        #for i in gene.transcription_region.region_paths:
        print("adding %s" % gene.transcription_region.region_paths)

    blocks = set(blocks)
    """
    trans_regions = [g.transcription_region for g in genes]
    subgraph, trans = graph.create_subgraph_from_intervals(trans_regions, 20000)
    #for g in genes:
    #    g.trans = trans.translate(g)

    # NB: remember to translate genes using trans
    for g in genes:
        g.transcription_region = trans.translate(g.transcription_region)
        for exon in g.exons:
            exon = trans.translate(exon)

    levels = Graph.level_dict(subgraph.blocks)

    #print(levels)
    #return

    # Find start block by choosing a block having no edges in
    start = None
    for b in subgraph.blocks:
        if len(subgraph.reverse_adj_list[b]) == 0:
            start = b
            break

    print("== Subgraph ==")
    print(subgraph)

    assert start is not None

    #return

    from offsetbasedgraph import VisualizeHtml
    subgraph.start_block = start
    max_offset = sum([subgraph.blocks[b].length() for b in subgraph.blocks])
    v = VisualizeHtml(subgraph, 0, max_offset , 0, levels, "", 800, genes)
    print(v.get_wrapped_html())




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Interact with a graph created from GRCh38')
    subparsers = parser.add_subparsers(help='Subcommands')

    # Subcommand for create graph
    parser_create_graph = subparsers.add_parser('create_graph', help='Create graph')
    parser_create_graph.add_argument('chrom_sizes_file_name',
                    help='Tabular file containing two columns, chrom/alt name and size')
    parser_create_graph.add_argument('alt_locations_file_name',
                    help='File containing alternative loci')
    parser_create_graph.add_argument('out_file_name',
                    help='Name of file to store graph and translation objects insize')
    parser_create_graph.set_defaults(func=create_graph)

    # Subcommand for genes
    parser_genes = subparsers.add_parser('check_duplicate_genes', help='Check duplicate genes')
    parser_genes.add_argument('translation_file_name',
                                help='Translation file created by running create_graph')
    parser_genes.add_argument('genes_file_name', help='Genes')
    parser_genes.set_defaults(func=check_duplicate_genes)

    # Subcommand for merge alt loci using alignments
    parser_merge_alignments = subparsers.add_parser('merge_alignment', help='Merge graph using alignments of alt locus')
    parser_merge_alignments.add_argument('translation_file_name',
                                help='Translation file created by running create_graph')
    parser_merge_alignments.add_argument('alt_locus_id', help='Id of alt locus (e.g. chr2_KI270774v1_alt')
    parser_merge_alignments.set_defaults(func=merge_alignment) # Subcommand for merge alt loci using alignments

    # Merge all alignments
    parser_merge_all_alignments = subparsers.add_parser('merge_all_alignments', help='Merge graph using alignments of ALL alt loci')
    parser_merge_all_alignments.add_argument('chrom_sizes_file_name',
                    help='Tabular file containing two columns, chrom/alt name and size')
    parser_merge_all_alignments.add_argument('out_file_name',
                                help='File name to store translation object for new graph')
    parser_merge_all_alignments.set_defaults(func=merge_all_alignments)

    # Visualize genes
    parser_visualize_genes = subparsers.add_parser('visualize_genes', help='Produce html visualization (that can be saved and opened in a browser)')
    parser_visualize_genes.add_argument('translation_file_name',
                    help='')
    parser_visualize_genes.add_argument('genes_file_name',
                                help='Pickled genes file')
    parser_visualize_genes.set_defaults(func=visualize_genes)



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)


    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.help()

    sys.exit()


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

