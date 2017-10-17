from offsetbasedgraph.graphtraverser import GraphTraverser, GraphTraverserUsingSequence
from offsetbasedgraph import Graph, Block, Interval, GraphWithReversals
from pyvg.sequences import SequenceRetriever
import unittest

class TestGraphTraverserUsingSequence(unittest.TestCase):

    def setUp(self):
        self.simple_graph = GraphWithReversals({
            1: Block(3),
            2: Block(3),
            3: Block(3),
            4: Block(3)
        },
        {
            1: [2, 4],
            2: [3]
        })

        self.search_sequence = "AAATTTGGG"
        self.sequence_retriever = SequenceRetriever(
            {1: "AAA",
             2: "TTT",
             3: "GGG",
             4: "CCC"
             }
        )

    def test_simple(self):

        traverser = GraphTraverserUsingSequence(self.simple_graph, self.search_sequence, self.sequence_retriever)
        traverser.search_from_node(1)
        self.assertEqual(traverser.get_nodes_found(), [1, 2, 3])

        traverser = GraphTraverserUsingSequence(self.simple_graph, "AAATTT", self.sequence_retriever)
        traverser.search_from_node(1)
        self.assertEqual(traverser.get_nodes_found(), [1, 2])

    def test_complex(self):
        graph = GraphWithReversals(
            {1: Block(3),
             2: Block(3),
             3: Block(3),
             4: Block(3),
             5: Block(3)
             },
            {1: [2, 4],
             2: [3],
             4: [5]
             }
        )
        search_sequence = "AAATTTGGG"
        sequence_retriever = SequenceRetriever(
            {1: "AAA",
             2: "TTT",
             3: "GGG",
             4: "TTT",
             5: "GGG"
             }
        )

        traverser = GraphTraverserUsingSequence(graph, search_sequence, sequence_retriever)
        traverser.search_from_node(1)
        self.assertEqual(traverser.get_nodes_found(), [1, 2, 3])




if __name__ == "__main__":
    unittest.main()