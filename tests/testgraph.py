import unittest
from offsetbasedgraph import Graph, Block, Interval, Position, Translation


def simple_graph():
    blocks = {_id: block for _id, block in
              enumerate([10, 20, 30, 40])}
    adj_list = {0: 1, 1: 2, 2: 3, 0: 4, 4: 3}
    return blocks, adj_list


def disjoint_graph():
    blocks = {_id: block for _id, block in
              enumerate([10, 30, 20, 40])}
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

    def test_merge_intervals(self):
        graph = disjoint_graph()
        interval_a = Interval(
            Position(1, 0),
            Position(1, 10))

        interval_b = Interval(
            Position(4, 0),
            Position(4, 10))
        translation = graph.merge_intervals(
            interval_a,
            interval_b)

        true_translation = Translation(
            {interval_a: -1,
             interval_b: -1},
            {-1: [interval_a, interval_b]}
            )

        self.assertEqual(translation,
                         true_translation)

    def test_connect(self):
        pass

    def test_from_file(self):
        pass
    
