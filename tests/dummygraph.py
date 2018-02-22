from offsetbasedgraph import Graph, Block, Interval, Position


def get_name_graph():
    blocks = {"A": Block(20), "B": Block(10)}
    return Graph(blocks, {})


def get_simple_graph():
    blocks = {1: Block(10), 2: Block(20), 3: Block(10), 4: Block(15)}
    block_edges = {1: [2, 3], 2: [4], 3: [4]}
    graph = Graph(blocks, block_edges)
    return graph


def get_mergable_graph():
    blocks = {1: Block(10), 2: Block(20), 3: Block(20), 4: Block(15)}
    block_edges = {1: [2, 3], 2: [4], 3: [4]}
    graph = Graph(blocks, block_edges)
    return graph


def get_insulated_graph():
    blocks = {1: Block(10),
              2: Block(10), 3: Block(10), 4: Block(10),
              5: Block(5), 6: Block(15), 7: Block(10),
              8: Block(10)}

    edges = {1: [2, 5],
             2: [3], 3: [4], 4: [8],
             5: [6], 6: [7], 7: [8]}

    return Graph(blocks, edges)


def get_disjoint_graph():
    blocks = {1: Block(10), 2: Block(20), 3: Block(30)}
    block_edges = {}
    return Graph(blocks, block_edges)


def get_realistic_graph():
    blocks = {i: Block(10) for i in range(1, 9)}
    edges = {i: [i+1] for i in range(1,6)}
    edges[1].append(7)
    edges[2].append(8)
    edges[7] = [4]
    edges[8] = [5]
    return Graph(blocks, edges)


def get_medium_complex_graph():
    blocks = {1: Block(10), 2: Block(20), 3: Block(10), 4: Block(15)}
