import unittest
from offsetbasedgraph import Graph, Block, Interval
from offsetbasedgraph.sequencegraph import SequenceGraph

class TestSequenceGraph(unittest.TestCase):

    def setUp(self):
        self.obgraph = Graph(
            {3: Block(4),
             4: Block(2),
             5: Block(5)}
            ,
            {
                3: [4],
                4: [5]
            }
        )
        self.obgraph.convert_to_numpy_backend()

        self.sequencegraph = SequenceGraph.create_empty_from_ob_graph(self.obgraph)

        self.obgraph2 = Graph(
            {1: Block(3),
             2: Block(3),
             3: Block(3)}
            ,
            {
                1: [2],
                2: [3]
            }
        )
        self.obgraph2.convert_to_numpy_backend()
        self.sequencegraph2 = SequenceGraph.create_empty_from_ob_graph(self.obgraph2)
        self.sequencegraph2.set_sequence(1, "aag")
        self.sequencegraph2.set_sequence(2, "gaa")
        self.sequencegraph2.set_sequence(3, "aga")
        self.retriever = self.sequencegraph2

    def test_set_get_sequence(self):
        self.sequencegraph.set_sequence(3, "aacc")
        seq = self.sequencegraph.get_sequence(3)
        self.assertEqual(seq, "aacc")

    def test_set_get_subsequence(self):
        self.sequencegraph.set_sequence(3, "aacc")
        seq = self.sequencegraph.get_sequence(3, 1)
        self.assertEqual(seq, "acc")

        seq = self.sequencegraph.get_sequence(3, 1, 2)
        self.assertEqual(seq, "a")

        seq = self.sequencegraph.get_sequence(3, 2, 3)
        self.assertEqual(seq, "c")

    def test_single_node_interval(self):
        interval = Interval(0, 3, [1])
        self.assertEqual(self.retriever.get_interval_sequence(interval),
                         "aag")

    def test_single_node_interval_with_offset(self):
        interval = Interval(1, 3, [1])
        self.assertEqual(self.retriever.get_interval_sequence(interval),
                         "ag")

    def test_single_node_interval_with_dual_offset(self):
        interval = Interval(1, 2, [1])
        self.assertEqual(self.retriever.get_interval_sequence(interval),
                         "a")

    def test_reversed_single_node_interval(self):
        interval = Interval(0, 3, [-1])
        self.assertEqual(self.retriever.get_interval_sequence(interval),
                         "ctt")

    def test_reversed_single_node_interval_with_dual_offsetl(self):
        interval = Interval(1, 2, [-3])
        self.assertEqual(self.retriever.get_interval_sequence(interval),
                         "c")

    def test_multiple_nodes_interval(self):
        interval = Interval(0, 3, [1, 2])
        self.assertEqual(self.retriever.get_interval_sequence(interval),
                         "aaggaa")

    def test_multiple_nodes_interval_second_rp_reversed(self):
        interval = Interval(0, 3, [1, -2])
        self.assertEqual(self.retriever.get_interval_sequence(interval),
                         "aagttc")

    def test_long_interval(self):
        interval = Interval(1, 1, [1, 2, 3, -3, -2, -1])
        self.assertEqual(self.retriever.get_interval_sequence(interval),
                         "aggaaagatctttcc")

    def test_to_from_file(self):
        self.sequencegraph.to_file("test.sequencegraph")
        other = SequenceGraph.from_file("test.sequencegraph")

        self.assertEqual(self.sequencegraph.get_sequence(3), other.get_sequence(3))
        self.assertEqual(self.sequencegraph.get_sequence(4), other.get_sequence(4))
        self.assertEqual(self.sequencegraph.get_sequence(5), other.get_sequence(5))

    def test_set_invalid_sequence(self):
        self.sequencegraph.set_sequence(3, "aryu")
        self.assertEqual(self.sequencegraph.get_sequence(3, 0, 4), "annn")

if __name__ == "__main__":
    unittest.main()