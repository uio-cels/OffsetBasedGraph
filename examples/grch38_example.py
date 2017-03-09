from offsetbasedgraph.graphcreators import create_initial_grch38_graph, convert_to_numeric_graph
graph = create_initial_grch38_graph("chrom.sizes.example")
numeric_graph, name_translation = convert_to_numeric_graph(graph)
print(numeric_graph)

from offsetbasedgraph.graphcreators import connect_without_flanks
connected_graph, connected_trans = connect_without_flanks(numeric_graph, "alt.loci.example", name_translation)
print(connected_graph)

from offsetbasedgraph.gene import GeneList
gene_list = GeneList.from_file("genes.example")
translated_gene_list = gene_list.translate(connected_trans)
