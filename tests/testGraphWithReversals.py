import unittest

from offsetbasedgraph import Graph, Block, Interval, Position, Translation, BlockCollection, GraphWithReversals

class TestGraphWithReversals(unittest.TestCase):

    def test_init(self):
        blocks = {1: Block(5), 2: Block(3)}
        graph = GraphWithReversals(blocks, {})

    def test_block_in_graph(self):
        blocks = {1: Block(3), 2: Block(3)}
        graph = GraphWithReversals(blocks, {})

        self.assertTrue(graph.block_in_graph(1))
        self.assertTrue(graph.block_in_graph(2))
        self.assertTrue(not graph.block_in_graph(3))

    def test_graph_without_reversals_identical_to_graph(self):
        blocks = {1: Block(3), 2: Block(3)}
        graph_with_reversals = GraphWithReversals(blocks, {2: [3]})
        graph = Graph(blocks, {2: [3]})

        self.assertEqual(graph_with_reversals, graph)


if __name__ == "__main__":
    unittest.main()
