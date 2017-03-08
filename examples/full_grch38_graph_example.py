# NB: This code takes time to run, as remote sequence data needs to be downloaded.

from offsetbasedgraph.graphutils import *
from offsetbasedgraph.gene import GeneList, GeneIO

graph = create_initial_grch38_graph("grch38.chrom.sizes")
numeric_graph, name_translation = convert_to_numeric_graph(graph)
new_numeric_graph, numeric_translation = connect_without_flanks(
    numeric_graph, "grch38_alt_loci.txt", name_translation)
name_graph, new_name_translation = convert_to_text_graph(
    new_numeric_graph, name_translation, numeric_translation)
final_translation = name_translation + numeric_translation + new_name_translation
final_translation.graph2 = name_graph
final_translation.to_file("grch38.graph")

grch38_graph_with_flanks = name_graph
print(grch38_graph_with_flanks.summary())

# Translate example genes
genes = GeneList.from_file("genes.example")
translated_genes = genes.translate(final_translation)
for gene in translated_genes.gene_list:
    print(GeneIO(gene).to_file_line())
