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

    def test_translate_rp(self):
        graph, graph2, trans = dummygraph.get_translation_single_block()
        self.assertEqual(trans.translate_rp(1)[0], trans._a_to_b[1][0])
        graph.blocks[5] = Block(2)
        self.assertEqual(trans.translate_rp(5)[0], Interval(0, 2, [5]))

    def test_translate_method(self):
        # Tests the wrapper method translate
        # Only simple cases to test that the wrapping works.
        # Testing of translation is not done here
        graph, graph2, trans = dummygraph.get_translation_single_block()
        position = Position(1, 0)
        self.assertEqual(trans.translate(position), Position(2, 0))
        interval = Interval(0, 1, [1])
        self.assertEqual(trans.translate(interval), Interval(0, 1, [2]))
        not_interval_or_pos = ""
        self.assertRaises(ValueError, trans.translate, not_interval_or_pos)

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

        # Case 4: Subgraph is not full graph
        graph1 = Graph(
            {
                1: Block(1),
                2: Block(1),
                3: Block(1)},
            {
                1: [2],
                2: [3]
            }
        )
        graph2 = Graph(
            {
                10: Block(1),
                20: Block(1),
                30: Block(1)},
            {
                10: [20],
                20: [30]
            }
        )
        trans = Translation(
            {
                1: [Interval(0, 1, [10], graph2)],
                2: [Interval(0, 1, [20], graph2)],
                3: [Interval(0, 1, [30], graph2)],
            },
            {
                10: [Interval(0, 1, [1], graph2)],
                20: [Interval(0, 1, [2], graph2)],
                30: [Interval(0, 1, [3], graph2)],
            }
        )
        subgraph = graph1 = Graph(
            {
                2: Block(1),
                3: Block(1)},
            {
                2: [3]
            }
        )
        correct_translated = Graph(
            {
                20: Block(1),
                30: Block(1)},
            {
                20: [30]
            }
        )

        translated_graph = trans.translate_subgraph(subgraph)
        self.assertEqual(correct_translated, translated_graph)

    def test_empty_start(self):
        pass  # TODO

    def test_overlapping_merge(self):
        graph = dummygraph.get_disjoint_graph()

        intervals_a = [Interval(5, 7, [2], graph=graph),
                       Interval(14, 16, [3], graph=graph)]

        intervals_b = [Interval(2, 8, [1], graph=graph),
                       Interval(12, 18, [3], graph=graph)]

        new_graph, trans = graph.merge(intervals_a)

        trans.graph2 = new_graph
        intervals_t = [trans.translate(i) for i in intervals_b]

        last_graph, new_trans = new_graph.merge(intervals_t)
        new_trans.graph2 = last_graph
        full_trans = trans + new_trans

        a, b, c, d, e = full_trans.translate(Interval(0, 30, [3])).region_paths
        f, b2, c2, d2, g = full_trans.translate(Interval(0, 10, [1])).region_paths
        h, c3, i = full_trans.translate(Interval(0, 20, [2])).region_paths

        self.assertEqual(b, b2)
        self.assertEqual(c, c2)
        self.assertEqual(c, c3)
        self.assertEqual(d, d2)

        t_blocks = {a: Block(12),
                    b: Block(2),
                    c: Block(2),
                    d: Block(2),
                    e: Block(12),
                    f: Block(2),
                    g: Block(2),
                    h: Block(5),
                    i: Block(13)}

        t_edges = {a: [b], b: [c], c: [d, i], d: [e, g],
                   f: [b], h: [c]}

        self.assertEqual(Graph(t_blocks, t_edges), last_graph)

    def test_overlapping_merge2(self):
        graph = dummygraph.get_disjoint_graph()
        intervals_a = [Interval(4, 7, [2], graph=graph),
                       Interval(11, 14, [3], graph=graph)]

        intervals_b = [Interval(2, 8, [1], graph=graph),
                       Interval(12, 18, [3], graph=graph)]

        new_graph, trans = graph.merge(intervals_a)

        trans.graph2 = new_graph
        intervals_t = [trans.translate(i) for i in intervals_b]

        last_graph, new_trans = new_graph.merge(intervals_t)
        new_trans.graph2 = last_graph
        full_trans = trans + new_trans
        a, b, c, d, e = full_trans.translate(Interval(0, 30, [3])).region_paths
        f, b2, c2, g = full_trans.translate(Interval(0, 20, [2])).region_paths
        h, c3, d3,  i = full_trans.translate(Interval(0, 10, [1])).region_paths

        self.assertEqual(b, b2)
        self.assertEqual(c, c2)
        self.assertEqual(c, c3)
        self.assertEqual(d, d3)

        t_blocks = {a: Block(11),
                    b: Block(1),
                    c: Block(2),
                    d: Block(4),
                    e: Block(12),
                    f: Block(4),
                    g: Block(13),
                    h: Block(2),
                    i: Block(2)}

        t_edges = {a: [b], b: [c], c: [d, g], d: [e, i],
                   f: [b], h: [c]}

        self.assertEqual(Graph(t_blocks, t_edges), last_graph)


    def test_name_translation(self):
        graph = dummygraph.get_name_graph()
        a_to_b = {k: i for i, k in enumerate(graph.blocks.keys())}
        trans = Translation.make_name_translation(a_to_b, graph)
        for k, v in a_to_b.items():
            self.assertEqual(trans._b_to_a[v][0].region_paths[0],
                             k)

        num_graph = trans.translate_subgraph(graph)
        trans.graph2 = num_graph
        intervals = [Interval(5, 15, ["A"], graph),
                     Interval(0, 10, ["B"], graph)]

        for interval in intervals:
            self.assertEqual(
                trans.translate(interval),
                Interval(interval.start_position.offset,
                         interval.end_position.offset,
                         [a_to_b[rp] for rp in interval.region_paths]))

        num_graph2, num_trans = num_graph.merge(
            [trans.translate(interval) for interval in intervals])

    def test_to_from_file(self):
        graph1, graph2, trans = dummygraph.get_merged_middle_translation()
        trans.to_file("test_trans")
        t2 = Translation.from_file("test_trans")
        self.assertEqual(trans, t2)
        self.assertEqual(trans.graph1, t2.graph1)
        self.assertEqual(trans.graph2, t2.graph2)

    def test_translate_subgraph_steps(self):

        graph = Graph({
            0: Block(1),
            1: Block(1),
            2: Block(1),
            3: Block(1),
            4: Block(1),
        },
            {
                0: [1],
                1: [2],
                2: [3],
                3: [4]
            })

        ab = {
            1: [Interval(0, 1, [10])],
            2: [Interval(0, 1, [11])]
        }
        ba = {
            10: [Interval(0, 1, [1])],
            11: [Interval(0, 1, [2])]
        }

        trans = Translation(ab, ba, graph=graph)
        # Old edges
        old_edges, new = trans.get_old_edges(graph)
        self.assertTrue(3 in old_edges)
        self.assertTrue(4 in old_edges[3])


        trans.get_external_edges(graph, old_edges, new)
        self.assertEqual(old_edges[0], [10  ])
        self.assertTrue(11 in old_edges[10])
        self.assertEqual(old_edges[11], [3])
        self.assertEqual(old_edges[3], [4])

    def test_translation_equal(self):
        graph1, graph2, trans = dummygraph.get_merged_middle_translation()

        self.assertEqual(trans, trans.copy())
        trans2 = Translation({}, {})
        self.assertTrue(trans != trans2)

        trans2 = trans.copy()
        trans2._a_to_b[1] = [Interval(1, 1, [1, 5, 3], graph2)]

        self.assertTrue(trans != trans2)

        trans2 = trans.copy()
        trans2._a_to_b[1] = [Interval(0, 1, [1, 5, 3], graph2)]
        self.assertTrue(trans == trans2)

        new_ab =  {
            2: [Interval(0, 1, [2, 5, 4], graph2)],
            1: [Interval(0, 1, [1, 5, 3], None)]
        }

        trans2 = trans.copy()
        trans2._a_to_b = new_ab

        self.assertEqual(trans, trans2)


    def test_translate_interval_special_case(self):
        # Special case with start of interval not being in start of RP
        graph = Graph(
            {1: Block(5),
             2: Block(5),
             3: Block(5)
            },
            {
             1: [2],
             2: [3],
            })

        ab = {
            1: [Interval(0, 1, [11, 12, 13, 14, 15])],
            2: [Interval(0, 5, [16])],
            3: [Interval(0, 1, [21, 22, 23, 24, 25])]
        }

        ba = {
            11: [Interval(0, 1, [1])],
            12: [Interval(1, 2, [1])],
            13: [Interval(2, 3, [1])],
            14: [Interval(3, 4, [1])],
            15: [Interval(4, 5, [1])],
            16: [Interval(0, 5, [2])],
            21: [Interval(0, 1, [3])],
            22: [Interval(1, 2, [3])],
            23: [Interval(2, 3, [3])],
            24: [Interval(3, 4, [3])],
            25: [Interval(4, 5, [3])]
        }

        trans = Translation(ab, ba, graph=graph)

        graph2 = trans.translate_subgraph(graph)
        trans.graph2 = graph2

        # Translate intervals starting in rp 1 with different offsets
        for i in range(0, 5):
            interval = Interval(i, 3, [1, 2, 3], graph)
            correct_length = (5-i) + 5 + 3
            self.assertEqual(interval.length(), correct_length)

            translated = trans.translate(interval)
            self.assertEqual(translated.length(), correct_length)
            self.assertEqual(translated.start_position.region_path_id, 11 + i)

        # Change end pos
        for i in range(0, 5):
            interval = Interval(2, i+1, [1, 2, 3], graph)
            correct_length = i + 1 + 5 + 3
            self.assertEqual(interval.length(), correct_length)

            translated = trans.translate(interval)
            self.assertEqual(translated.length(), correct_length)
            self.assertEqual(translated.start_position.region_path_id, 13)
            self.assertEqual(translated.end_position.region_path_id, 21 + i)

    def test_str(self):
        graph1, graph2, trans = dummygraph.get_merged_middle_translation()
        str(trans)
        self.assertTrue(True)

    def test_translate_interval_nontrivial_unequal_starts_and_ends(self):
        # Special case where start of interval translates to
        # more/less intervals than end of interval
        graph = Graph(
            {
                1: Block(3),
                2: Block(3)
            },
            {
                1: [2]
            }
        )

        graph2 = Graph(
            {
                3: Block(3),
                4: Block(3),
                5: Block(3)
            },
            {
                3: [5],
                4: [5]
            }
        )
        ab = {
                1: [Interval(0, 3, [3], graph2), Interval(0, 3, [4], graph2)],
                2: [Interval(0, 3, [5], graph2)]
            }

        ba = {
                3: [Interval(0, 3, [1], graph)],
                4: [Interval(0, 3, [1], graph)],
                5: [Interval(0, 3, [2], graph)],
            }

        trans = Translation(ab, ba, graph)

        self.assertEqual(len(trans.translate_position(Position(1, 0))), 2)

        interval = Interval(0, 3, [1, 2], graph)
        translated = trans.translate_interval(interval).get_single_path_intervals()
        self.assertEqual(len(translated), 2)

        correct_translations = [
                Interval(0, 3, [3, 5]),
                Interval(0, 3, [4, 5])
        ]
        self.assertTrue(correct_translations[0] in translated)
        self.assertTrue(correct_translations[1] in translated)


        # Inverse case
        trans = Translation(ba, ab, graph2)
        interval = Interval(0, 3, [1, 2], graph)
        translated = trans.translate_interval(interval, True).get_single_path_intervals()
        self.assertEqual(len(translated), 2)


if __name__ == "__main__":
    unittest.main()
