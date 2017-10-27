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
        print(self.simple_graph.reverse_adj_list)
        traverser = GraphTraverserUsingSequence(self.simple_graph, self.search_sequence, self.sequence_retriever)
        traverser.search_from_node(1)
        self.assertEqual(traverser.get_nodes_found(), [1, 2, 3])

    def test_simple2(self):
        traverser = GraphTraverserUsingSequence(self.simple_graph, "AAATTT", self.sequence_retriever)
        traverser.search_from_node(1)
        self.assertEqual(traverser.get_nodes_found(), [1, 2])

    def test_complex(self):
        graph = GraphWithReversals(
            {1: Block(3),
             2: Block(3),
             3: Block(3),
             4: Block(3),
             5: Block(3),
             6: Block(3)
             },
            {1: [2, 4],
             2: [3],
             4: [5],
             5: [6]
             }
        )
        search_sequence = "AAATTTGGG"
        sequence_retriever = SequenceRetriever(
            {1: "AAA",
             2: "TTT",
             3: "GAG",
             4: "TTT",
             5: "GGG",
             6: "CGG"
             }
        )

        traverser = GraphTraverserUsingSequence(graph, search_sequence, sequence_retriever)
        traverser.search_from_node(1)
        self.assertEqual(traverser.get_nodes_found(), [1, 4, 5])

    def test_many_paths_to_same_node(self):
        blocks= {i: Block(1) for i in range(1, 6)}
        edges = {
                    1: [2, 5],
                    2: [3, 5],
                    3: [4, 5],
                    4: [5]
        }
        graph = GraphWithReversals(blocks, edges)

        search_sequence = "AAAG"
        sequence_retriever = SequenceRetriever(
            {1: "A",
             2: "A",
             3: "A",
             4: "A",
             5: "G"
             }
        )

        traverser = GraphTraverserUsingSequence(graph, search_sequence, sequence_retriever)
        traverser.search_from_node(1)
        self.assertEqual(traverser.get_nodes_found(), [1, 2, 3, 5])


    def test_reversal(self):
        blocks= {i: Block(1) for i in range(1, 6)}
        edges = {
                    1: [2, 5],
                    2: [-2, 3],
                    -2: [5],
                    3: [4],
                    4: [5]
        }
        graph = GraphWithReversals(blocks, edges)

        search_sequence = "AAAG"
        sequence_retriever = SequenceRetriever(
            {1: "A",
             2: "A",
             3: "T",
             4: "A",
             5: "G"
             }
        )

        traverser = GraphTraverserUsingSequence(graph, search_sequence, sequence_retriever)
        traverser.search_from_node(1)
        self.assertEqual(traverser.get_nodes_found(), [1, 2, -2, 5])

if __name__ == "__main__":
    unittest.main()