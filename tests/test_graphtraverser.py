from offsetbasedgraph.graphtraverser import GraphTraverser, GraphTraverserUsingSequence, GraphTravserserBetweenNodes
from offsetbasedgraph import Graph, Block, Interval, GraphWithReversals
from pyvg.sequences import SequenceRetriever
from offsetbasedgraph import SequenceGraph
import unittest

class TestGraphTravserserBetweenNodes(unittest.TestCase):

    def setUp(self):
        self.simple_graph = Graph(
            {i: Block(3) for i in range(1, 9)},
            {
                1: [2, 3],
                 2: [4],
                 3: [4],
                 4: [5],
                 5: [6, 7],
                 6: [8],
                 7: [8]
             })

        self.hierarchical_graph = Graph(
            {i: Block(3) for i in range(1, 13)},
            {
                11: [1],
                1: [2, 3],
                2: [7, 8],
                3: [4, 5],
                4: [6],
                5: [6],
                6: [10],
                7: [9],
                8: [9],
                9: [10],
                10: [12]
             })

    def test_snarl_traverser(self):
        graph = Graph(
            {
                1: Block(3),
                2: Block(3),
                3: Block(3),
                4: Block(3)
            },
            {
                1: [2, 4],
                2: [3],
                4: [3]
            }
        )

        traverser = GraphTravserserBetweenNodes(graph)
        subgraph = traverser.get_snarl_subgraph(1, 3, print_debug=True)
        correct_subgraph = Graph({2: Block(3), 4: Block(3)}, {})
        self.assertEqual(correct_subgraph, subgraph)

    def test_snarl_traverser2(self):


        traverser = GraphTravserserBetweenNodes(self.hierarchical_graph)
        subgraph = traverser.get_snarl_subgraph(3, 6, print_debug=True)
        correct_subgraph = Graph({4: Block(3), 5: Block(3)}, {})
        self.assertEqual(correct_subgraph, subgraph)

    def test_snarl_traverser_hierarchical(self):
        traverser = GraphTravserserBetweenNodes(self.hierarchical_graph)
        subgraph = traverser.get_snarl_subgraph(1, 10, include_start_and_end=True, print_debug=True)
        correct_subgraph = Graph({i: Block(3) for i in range(1, 11)},
                                 {
                                     1: [2, 3],
                                     2: [7, 8],
                                     3: [4, 5],
                                     4: [6],
                                     5: [6],
                                     7: [9],
                                     8: [9],
                                     6: [10],
                                     9: [10]
                                 })

        print("Subgraph")
        print(subgraph)

        self.assertEqual(correct_subgraph, subgraph)



    def test_simple_snarl_traverser(self):
        graph = Graph({i: Block(3) for i in range(1, 7)},
                      {
                          1: [2, 3, 5],
                          2: [4],
                          3: [4],
                          5: [6]
                      })

        traverser = GraphTravserserBetweenNodes(graph)
        subgraph = traverser.get_snarl_subgraph(1, 4)
        correct_subgraph = Graph({2: Block(3), 3: Block(3)}, {})
        self.assertEqual(correct_subgraph, subgraph)

        traverser = GraphTravserserBetweenNodes(graph)
        subgraph = traverser.get_snarl_subgraph(1, 6, include_start_and_end=True)
        correct_subgraph = Graph({1: Block(3), 5: Block(3), 6: Block(3)},
                                 {1: [5],
                                  5: [6]
                                  })

        self.assertEqual(correct_subgraph, subgraph)

    def test_simple(self):
        traverser = GraphTravserserBetweenNodes(self.simple_graph)
        subgraph = traverser.get_greedy_subgraph_between_nodes(1, 4)
        correct_subgraph = Graph({i: Block(3) for i in range(1, 5)},
                                 {
                                    1: [2, 3],
                                    2: [4],
                                    3: [4]
                                 })
        self.assertEqual(correct_subgraph, subgraph)

    def test_simple2(self):

        traverser = GraphTravserserBetweenNodes(self.simple_graph)
        subgraph = traverser.get_greedy_subgraph_between_nodes(1, 8)
        self.assertEqual(self.simple_graph, subgraph)


# Must be rewritten to use sequenceGraph
"""
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
        self.simple_graph.convert_to_numpy_backend()

        self.search_sequence = "aaatttggg"
        self.sequence_retriever = SequenceGraph.create_empty_from_ob_graph(self.simple_graph)
        self.sequence_retriever.set_sequence(1, "AAA")
        self.sequence_retriever.set_sequence(2, "TTT")
        self.sequence_retriever.set_sequence(3, "CCC")
        self.sequence_retriever.set_sequence(4, "GGG")


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

        search_sequence = "AATG"
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
"""
if __name__ == "__main__":
    unittest.main()