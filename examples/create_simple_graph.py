from offsetbasedgraph import Block, Graph
blocks = {"chr1-part1": Block(20), "chr1-part2": Block(10)}
adjencies = {"chr1-part1": ["chr1-part2"]}
graph = Graph(blocks, adjencies)
graph.to_file("mygraph.graph")


graph = Graph.from_file("mygraph.graph")