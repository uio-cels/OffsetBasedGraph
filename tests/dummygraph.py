from offsetbasedgraph import Graph, Block, Translation, Interval, Position


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
    blocks = {i: Block(10) for i in range(9)}
    edges = {i: [i+1] for i in range(6)}
    edges[1].append(7)
    edges[2].append(8)
    edges[7] = [4]
    edges[8] = [5]
    return Graph(blocks, edges)


def get_medium_complex_graph():
    blocks = {1: Block(10), 2: Block(20), 3: Block(10), 4: Block(15)}


# Graphs and translations
def get_translation_single_block():
    graph = Graph({1: Block(10)}, {})   # Simple graph with one block
    graph2 = Graph({2: Block(5), 3: Block(5)}, {2: [3]})

    # Intervals on graph 1
    interval1 = Interval(Position(1, 0), Position(1, 5), [1], graph)  # Sub of first block
    interval2 = Interval(Position(1, 5), Position(1, 10), [1], graph)  # Sub of first block

    # Interval for whole of graph 2
    interval3 = Interval(Position(2, 0), Position(3, 5), [2, 3], graph2)

    trans = Translation({1: [interval3]}, {2: [interval1], 3: [interval2]})

    return graph, graph2, trans


def get_merged_translation():
    graph1 = Graph({1: Block(10), 2: Block(10)}, {})  # Two disjoint blocks
    graph2 = Graph({3: Block(10)}, {})   # Simple graph with one block

    # Translation: Merge from the two blocks in graph 1 into one in graph2
    trans = Translation(
        {
            1: [Interval(0, 10, [3], graph2)],
            2: [Interval(0, 10, [3], graph2)]
        },
        {
            3: [Interval(0, 10, [1], graph1),
                Interval(0, 10, [2], graph1)]
        }
    )

    return graph1, graph2, trans


def get_merged_middle_translation():
    graph1 = Graph({1: Block(3), 2: Block(3)}, {})
    graph2 = Graph(
                {
                    1: Block(1),
                    2: Block(1),
                    3: Block(1),
                    4: Block(1),
                    5: Block(1)
                },
                {   1: [5],
                    2: [5],
                    5: [3, 4]
                }
    )
    trans = Translation(
        {
            1: [Interval(0, 1, [1, 5, 3], graph2)],
            2: [Interval(0, 1, [2, 5, 4], graph2)]
        },
        {
            1: [Interval(0, 1, [1], graph1)],
            2: [Interval(0, 1, [2], graph1)],
            5: [Interval(1, 2, [1], graph1), Interval(1, 2, [2], graph1)],
            3: [Interval(2, 3, [1], graph1)],
            4: [Interval(2, 3, [2], graph1)]
        }
    )
    return graph1, graph2, trans
