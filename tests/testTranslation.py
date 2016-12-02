import unittest
import numpy as np
from . import dummygraph
from offsetbasedgraph import Interval, Position, Graph, Translation, Block

def translation_single_block():
    graph = Graph({1: Block(10)}, {})   # Simple graph with one block
    graph2 = Graph({2: Block(5), 3: Block(5)}, {2: [3]})

    # Intervals on graph 1
    interval1 = Interval(Position(1, 0), Position(1, 5), [1])  # Sub of first block
    interval2 = Interval(Position(1, 5), Position(1, 10), [1])  # Sub of first block

    # Interval for whole of graph 2
    interval3 = Interval(Position(2, 0), Position(3, 5), [2, 3])

    trans = Translation({1: interval3}, {2: interval1, 3: interval2})

    return graph, graph2, trans

class TestTranslation(unittest.TestCase):

    def testCreateTranslation(self):
        graph, graph2, trans = translation_single_block()

        self.assertTrue(len(trans._a_to_b) == 1)
        self.assertTrue(len(trans._b_to_a) == 2)

    def testTranslatePosition(self):
        pos_graph1 = Position(1, 4)
        pos_graph2 = trans.translate_position()



if __name__ == "__main__":
    unittest.main()

