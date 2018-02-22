import unittest
import dummygraph
from offsetbasedgraph import Graph, Block, Interval, Position, BlockCollection


DEBUG = False


def simple_graph():
    blocks = {_id: block for _id, block in
              enumerate([10, 20, 30, 40])}
    adj_list = {0: [1], 1: [2], 2: [3], 0: [4], 4: [3]}
    return blocks, adj_list


def disjoint_graph():
    blocks = {_id: Block(block) for _id, block in
              enumerate([30, 40])}
    adj_list = {}
    return Graph(blocks, adj_list)


class TestGraph(unittest.TestCase):

    def test_init(self):
        blocks = {1: Block(5), 2: Block(3)}
        graph = Graph(blocks, {})

    def test_block_in_graph(self):
        blocks = {1: Block(3), 2: Block(3)}
        graph = Graph(blocks, {})

        self.assertTrue(graph.block_in_graph(1))
        self.assertTrue(graph.block_in_graph(2))
        self.assertTrue(not graph.block_in_graph(3))


    def assert_graph_equals(self, graph, blocks, adj_list):
        new_graph = Graph(blocks, adj_list)
        self.assertEqual(new_graph, graph)

    def test_graph_equals(self):
        # Case 1
        graph1 = Graph({1: Block(10)}, {})
        graph2 = Graph({1: Block(10)}, {})
        self.assertEqual(graph1, graph2)

        # Case 2
        graph1 = Graph({1: Block(10), 2: Block(1)}, {})
        graph2 = Graph({1: Block(10)}, {})
        self.assertTrue(graph1 != graph2)

        # Case 3
        graph1 = Graph({1: Block(10), 2: Block(1)}, {1: []})
        graph2 = Graph({1: Block(10), 2: Block(1)}, {})
        self.assertEqual(graph1, graph2)

        # Case 4
        graph1 = Graph({1: Block(10), 2: Block(1)}, {1: [2]})
        graph2 = Graph({1: Block(10), 2: Block(1)}, {1: [2]})
        self.assertEqual(graph1, graph2)

    def test_init(self):
        blocks, adj_list = simple_graph()
        graph = Graph(blocks, adj_list)
        self.assert_graph_equals(graph, blocks, adj_list)

    def test_reverse_edges(self):
        adj_list = {1: [10, 20],
                    2: [100],
                    3: [100]}

        reverse_list = Graph._get_reverse_edges(adj_list)
        fasit = {-10: [-1], -20: [-1], -100: [-2, -3]}
        self.assertEqual(reverse_list, fasit)

    def _setup_merge(self):
        graph = disjoint_graph()
        interval_a = Interval(
            Position(0, 0),
            Position(0, 10), graph=graph)

        interval_b = Interval(
            Position(1, 0),
            Position(1, 10), graph=graph)

        return graph, interval_a, interval_b

    def _test_get_all_block_borders(self):
        blocks = {1: Block(20), 2: Block(20),
                  11: Block(20), 12: Block(20)}

        adj_list = {1: [2], 11: [12]}
        graph = Graph(blocks, adj_list)

        interval_a = Interval(Position(1, 10), Position(2, 11))
        interval_b = Interval(Position(11, 15), Position(12, 16))
        border_list = graph._get_all_block_borders(interval_a, interval_b)
        self.assertEqual(border_list, [5, 10, 21])

    def _test_split_blocks_at_starts_and_ends(self):
        graph = dummygraph.get_disjoint_graph()
        intervals = [Interval(Position(1, 0), Position(1, 5), graph=graph),
                     Interval(Position(2, 5), Position(2, 10), graph=graph),
                     Interval(Position(3, 25), Position(3, 30), graph=graph)]

        trans = graph._split_blocks_at_starts_and_ends(intervals)

        rps1 = trans.translate_rp(1).region_paths
        rps2 = trans.translate_rp(2).region_paths
        rps3 = trans.translate_rp(3).region_paths

        # Check that starts and ends have not been split
        self.assertEqual(len(rps1), 2)
        self.assertEqual(len(rps2), 3)
        self.assertEqual(len(rps3), 2)

        true_blocks = dict(zip(rps1, [Block(5), Block(5)]))
        true_blocks.update(dict(zip(rps2, [Block(5), Block(5), Block(10)])))
        true_blocks.update(dict(zip(rps3, [Block(25), Block(5)])))

        true_adjs = {rps1[0]: [rps1[1]],
                     rps2[0]: [rps2[1]], rps2[1]: [rps2[2]],
                     rps3[0]: [rps3[1]]}
        self.assert_graph_equals(graph, true_blocks, true_adjs)

    def test_connect_intervals(self):
        pass

    def test_from_file(self):
        pass

    def test_has_identical_structure(self):
        # Case 1
        g1 = Graph(
            {
                1: Block(1),
                2: Block(10)
            },
            {
                1: [2]
            }
        )

        g2 = Graph(
            {
                5: Block(10),
                2: Block(1)
            },
            {
                5: [2]
            }
        )

        self.assertTrue(g1.has_identical_structure(g2))

        # Case 2
        g1 = Graph(
            {
                1: Block(1),
                2: Block(1),
                3: Block(1),
                4: Block(1)
            },
            {
                1: [2, 3],
                3: [4],
                2: [4]
            }
        )

        g2 = Graph(
            {
                10: Block(2),
                20: Block(2),
                30: Block(2),
                40: Block(2)
            },
            {
                10: [20, 30],
                30: [40],
                20: [40]
            }
        )

        self.assertTrue(g1.has_identical_structure(g2))

    def test_find_critical_blocks(self):
        pass
        #graph = dummygraph.get_realistic_graph()
        #critical_blocks = graph.find_critical_blocks(0)
        #self.assertEqual(critical_blocks, [1, 5, 6])

    def _test_to_from_file(self):
        for graph in [dummygraph.get_simple_graph(),
                      dummygraph.get_disjoint_graph()]:
            graph.to_file("test_graph")

            graph2 = Graph.from_file("test_graph")
            self.assertEqual(graph, graph2)

    def test_block_collection(self):
        blocks = BlockCollection({1: Block(10), 2: Block(12)})

        self.assertTrue(1 in blocks)
        self.assertTrue(2 in blocks)
        self.assertTrue(3 not in blocks)

if __name__ == "__main__":
    unittest.main()
