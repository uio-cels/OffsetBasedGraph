import unittest
from offsetbasedgraph import GraphWithReversals, Block
from offsetbasedgraph.graph import BlockArray, AdjListAsNumpyArrays

class TestBackendConversion(unittest.TestCase):

    def test_simple(self):

        graph = GraphWithReversals({
            1: Block(5),
            2: Block(10),
            3: Block(15)
        }, {1: [2], 2: [3]})

        graph2 = GraphWithReversals({
            1: Block(5),
            2: Block(10),
            3: Block(15)
        }, {1: [2], 2: [3]})

        graph2.convert_to_numpy_backend()
        self.assertTrue(1 in graph2.blocks)
        self.assertTrue(isinstance(graph2.blocks, BlockArray))
        self.assertTrue(isinstance(graph2.adj_list, AdjListAsNumpyArrays))
        self.assertEqual(graph2, graph)

        graph2.convert_to_dict_backend()
        self.assertEqual(graph, graph2)

class TestToFromSingleNumpyFile(unittest.TestCase):

    def test_simple(self):
        graph = GraphWithReversals({
            1: Block(5),
            2: Block(10),
            3: Block(15)
        }, {1: [2], 2: [3]})

        graph.convert_to_numpy_backend()
        graph.to_file("test.nobg")

        graph2 = GraphWithReversals.from_file("test.nobg")

        self.assertEqual(graph, graph2)

if __name__ == "__main__":
    unittest.main()
