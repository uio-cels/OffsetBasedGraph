import unittest
import dummygraph
from offsetbasedgraph import Interval, Position, Graph, Block, Translation


class TestTranslation(unittest.TestCase):

    def testCreateTranslation(self):
        graph, graph2, trans = dummygraph.get_translation_single_block()

        self.assertTrue(len(trans._a_to_b) == 1)
        self.assertTrue(len(trans._b_to_a) == 2)

        self.assertEqual(trans.graph2, graph2)
        self.assertEqual(trans.graph1, graph)

    def testSimpleTranslatePosition(self):
        graph, graph2, trans = dummygraph.get_translation_single_block()
        pos_graph1 = Position(1, 4)

        # Case 1
        pos_graph2 = trans.translate_position(pos_graph1)
        self.assertTrue(pos_graph2[0], Position(2, 4))

        pos_graph1_back = trans.translate_position(pos_graph2[0], True)
        self.assertEqual(pos_graph1_back[0], pos_graph1)

        # Run some other cases back and forth
        for i in range(0, 10):
            self.assertEqual(trans.translate_position(
                                trans.translate_position(pos_graph1)[0],
                                True
                            )[0], pos_graph1,
                            "Position %s not translated correctly back and forth"
                            % pos_graph1
            )

    def testSimpleTranslateInterval(self):
        graph, graph2, trans = dummygraph.get_translation_single_block()
        interval_graph1 = Interval(Position(1, 3), Position(1, 8), [1], graph)
        interval_graph2 = Interval(Position(2, 3), Position(3, 3), [2, 3], graph2)
        translated = trans.translate_interval(interval_graph1)
        #print("Translated to %s" % translated)
        single_path_intervals = translated.get_single_path_intervals()

        self.assertEqual(len(single_path_intervals), 1,
                         "Interval should be translated to 1 other interval")

        translated_intervals = translated.get_single_path_intervals()
        t = translated_intervals[0]
        self.assertEqual(t, interval_graph2,
                "Translated interval %s not equal to %s" % (t, interval_graph2))
        translated_back = trans.translate_interval(t, True)
        t_back = translated_back.get_single_path_intervals()[0]

        self.assertEqual(t_back, interval_graph1,
                "Translated back interval %s != to %s" % (t_back, interval_graph1))

    def test_translate_interval_special_case(self):
        print("=== SPecial case === ")
        graph = Graph({1: Block(4)}, {})
        graph2 = Graph({2: Block(2), 3: Block(2)}, {2: [3]})
        trans = Translation(
                    {1: [Interval(0, 2, [2, 3], graph2)]},
                    {2: [Interval(0, 2, [1], graph)],
                     3: [Interval(2, 4, [1], graph)]})
        interval = Interval(3, 4, [1], graph)

        translated = trans.translate_interval(interval).get_single_path_intervals()[0]
        correct = Interval(1, 2, [3], graph2)

        self.assertEqual(correct, translated, "correct %s != %s")



    def test_translate_on_merged_graph(self):
        graph1, graph2, trans = dummygraph.get_merged_translation()

        interval1 = Interval(0, 10, [1], graph1)
        correct_tran = Interval(0, 10, [3], graph2)

        translated = trans.translate_interval(interval1).get_single_path_intervals()[0]
        self.assertEqual(translated, correct_tran)

        # translate back
        back = trans.translate_interval(translated, True).get_single_path_intervals()
        self.assertEqual(len(back), 2)

        self.assertTrue(back[0] == interval1 or back[1] == interval1)

    def test_add_translation_on_merged_graph(self):
        graph1, graph2, trans = dummygraph.get_merged_translation()
        interval1 = Interval(0, 10, [3], graph2)

        # Translate back to original graph (block 1)

        #graph3 = Graph({1: Block(10)}, {})
        trans_back = Translation(
            {3: [Interval(0, 10, [1], graph1)]},
            {1: [Interval(0, 10, [3], graph2)]}
        )

        trans_sum = trans + trans_back
        interval_back = trans_sum.translate_interval(interval1).get_single_path_intervals()
        self.assertTrue(interval_back, interval1)


    def test_translate_intervals_forth_and_back(self):
        pass

    def testAddTranslation(self):
        # Scenario: Splitting one region path two times
        graph, graph2, trans = dummygraph.get_translation_single_block()
        graph3 = Graph({4: Block(3), 5: Block(2), 3: Block(5)}, {4: [5], 5: [3]})

        # Split first region path of graph2 again
        intervalgraph3 = Interval(0, 2, [4, 5], graph3)  # Block 3 and 4 have length 5 in total
        trans2 = Translation({2: [intervalgraph3]},
                             {4: [Interval(0, 3, [2], graph2)],
                              5: [Interval(3, 5, [2], graph2)]}
                            )
        trans3 = trans + trans2
        correct_trans = Translation({1: [Interval(0, 5, [4, 5, 3], graph3)]},
                                    {4: [Interval(0, 3, [1], graph)],
                                     5: [Interval(3, 5, [1], graph)],
                                     3: [Interval(5, 10,[1], graph)]
                                    })

        self.assertTrue(trans3, correct_trans)

    def test_translate_subgraph(self):
        # Case 1
        graph1, graph2, trans = dummygraph.get_translation_single_block()
        translated_graph = trans.translate_subgraph(graph1)
        self.assertEqual(graph2, translated_graph)

        # Case 2
        graph1, graph2, trans = dummygraph.get_merged_translation()
        translated_graph = trans.translate_subgraph(graph1)
        self.assertEqual(graph2, translated_graph)

        # Case 3
        graph1, graph2, trans = dummygraph.get_merged_middle_translation()
        translated_graph = trans.translate_subgraph(graph1)
        self.assertEqual(graph2, translated_graph)



if __name__ == "__main__":
    unittest.main()

