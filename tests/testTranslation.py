import unittest
import numpy as np
from . import dummygraph
from offsetbasedgraph import Interval, Position, Graph, Translation, Block

def get_translation_single_block():
    graph = Graph({1: Block(10)}, {})   # Simple graph with one block
    graph2 = Graph({2: Block(5), 3: Block(5)}, {2: [3]})

    # Intervals on graph 1
    interval1 = Interval(Position(1, 0), Position(1, 5), [1])  # Sub of first block
    interval2 = Interval(Position(1, 5), Position(1, 10), [1])  # Sub of first block

    # Interval for whole of graph 2
    interval3 = Interval(Position(2, 0), Position(3, 5), [2, 3])

    trans = Translation({1: [interval3]}, {2: [interval1], 3: [interval2]})

    return graph, graph2, trans


class TestTranslation(unittest.TestCase):

    def testCreateTranslation(self):
        graph, graph2, trans = get_translation_single_block()

        self.assertTrue(len(trans._a_to_b) == 1)
        self.assertTrue(len(trans._b_to_a) == 2)

    def testSimpleTranslatePosition(self):
        graph, graph2, trans = get_translation_single_block()
        pos_graph1 = Position(1, 4)

        # Case 1
        pos_graph2 = trans.translate_position(pos_graph1)
        self.assertTrue(pos_graph2, Position(2, 4))

        pos_graph1_back = trans.translate_position(pos_graph2, True)
        self.assertEqual(pos_graph1_back, pos_graph1)

        # Run some other cases back and forth
        for i in range(0, 10):
            self.assertEqual(trans.translate_position(
                                trans.translate_position(pos_graph1),
                                True
                            ), pos_graph1,
                            "Position %s not translated correctly back and forth"
            )

    def testSimpleTranslateInterval(self):
        graph, graph2, trans = get_translation_single_block()
        interval_graph1 = Interval(Position(1, 3), Position(1, 8), [1])
        interval_graph2 = Interval(Position(2, 3), Position(3, 3), [2, 3])
        translated = trans.translate_interval(interval_graph1)
        self.assertEqual(translated, interval_graph1,
                    "Translated interval %s not equal to %s" % (translated,
                                                               interval_graph2))
        translated_back = trans.translate_interval(translated, True)
        self.assertEqual(translated_back, interval_graph1,
                    "Translated back interval %s not equal to %s" %
                         (translated_back, interval_graph1))

    def testAddTranslation(self):
        # Scenario: Splitting one region path two times
        graph, graph2, trans = get_translation_single_block()
        # Split first region path of graph2 again
        intervalgraph3 = Interval(0, 2, [4, 5])  # Block 3 and 4 have length 5 in total
        trans2 = Translation({2: [intervalgraph3]},
                             {4: [Interval(0, 3, [2])],
                              5: [Interval(3, 5, [2])]}
                            )
        trans = trans + trans2
        correct_trans = Translation({1: [Interval(0, 5, [4, 5, 3])]},
                                    {4: [Interval(0, 3, [1])],
                                     5: [Interval(3, 5, [1])],
                                     3: [Interval(5, 10,[1])]
                                    })

        self.assertTrue(trans, correct_trans)







if __name__ == "__main__":
    unittest.main()

