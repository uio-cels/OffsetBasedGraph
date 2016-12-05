from offsetbasedgraph import Graph, Block

def get_simple_graph():
    blocks = {1: Block(10), 2: Block(20), 3: Block(10), 4: Block(15)}
    block_edges = {1: [2], 1: [3], 2: [4], 3: [4]}
    graph = Graph(blocks, block_edges)
    return graph



def get_medium_complex_graph():
    blocks = {1: Block(10), 2: Block(20), 3: Block(10), 4: Block(15)}

