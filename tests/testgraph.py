import unittest
from offsetbasedgraph import Graph, Block, Interval, Position, Translation


def simple_graph():
    blocks = {_id: block for _id, block in
              enumerate([10, 20, 30, 40])}
    adj_list = {0: [1], 1: [2], 2: [3], 0: [4], 4: [3]}
    return blocks, adj_list


def disjoint_graph():
    blocks = {_id: block for _id, block in
              enumerate([30, 40])}
    adj_list = {}
    return Graph(blocks, adj_list)


class TestGraph(unittest.TestCase):

    def assert_graph_equals(self, graph, blocks, adj_list):
        self.assertEqual(graph.blocks, blocks)
        self.assertEqual(graph.adj_list, adj_list)
        self.assertEqual(graph.reverse_adj_list,
                         graph._get_reverse_edges(adj_list))

    def test_init(self):
        blocks, adj_list = simple_graph()
        graph = Graph(blocks, adj_list)
        self.assert_graph_equals(graph, blocks, adj_list)

    def test_reverse_edges(self):
        adj_list = {1: [10, 20],
                    2: [100],
                    3: [100]}

        reverse_list = Graph._get_reverse_edges(adj_list)
        fasit = {10: [1], 20: [1], 100: [2, 3]}
        self.assertEqual(reverse_list, fasit)

    def _setup_merge(self):
        graph = disjoint_graph()
        interval_a = Interval(
            Position(0, 0),
            Position(0, 10))

        interval_b = Interval(
            Position(1, 0),
            Position(1, 10))
        return graph, interval_a, interval_b

    def test_merge_translation(self):
        graph, interval_a, interval_b = self._setup_merge()
        translation = graph.merge_intervals(
            interval_a,
            interval_b)

        A = 0
        B = 1
        a_first = translation._a_to_b[A].region_paths[0]
        a_last = translation._a_to_b[A].region_paths[1]
        b_first = translation._a_to_b[B].region_paths[0]
        b_last = translation._a_to_b[B].region_paths[1]

        # Sanity check on merge
        self.assertEqual(a_first, b_first)
        self.assertEqual(a_first.length() == interval_a.length())

        # Check that translation object is correct
        a_to_b = {
            A: Interval(Position(a_first, 0),
                        Position(a_last, 20)),
            B: Interval(Position(b_first, 0),
                        Position(b_last, 30))
            }
        b_to_a = {
            a_first: set(interval_a, interval_b),
            a_last: set(
                Interval(
                    Position(A, 10),
                    Position(A, 30)),
                Interval(
                    Position(B, 10),
                    Position(B, 40))
                )
            }
        self.assertEqual(
            Translation(a_to_b, b_to_a),
            translation)

        # Check graph
        blocks = {a_first: Block(10),
                  a_last: Block(20),
                  b_last: Block(30)
                  }
        adj_list = {a_first: [a_last,  b_last]}
        self.assert_graph_equals(
            graph, Graph(blocks, adj_list))

    def test_get_all_block_borders(self):
        blocks = {1: Block(20), 2: Block(20),
                  11: Block(20), 12: Block(20)}

        adj_list = {1: [2], 11: [12]}
        graph = Graph(blocks, adj_list)

        interval_a = Interval(Position(1, 10), Position(2, 11))
        interval_b = Interval(Position(11, 15), Position(12, 16))
        border_list = graph._get_all_block_borders(interval_a, interval_b)
        self.assertEqual(border_list, [5, 10, 21])

    def test_connect_intervals(self):
        pass

    def test_from_file(self):
        pass


if __name__ == "__main__":
    unittest.main()
