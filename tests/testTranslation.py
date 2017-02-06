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

    def test_empty_start(self):
        pass  # TODO

    def _test_overlapping_merge(self):
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

        print("=== trans: ===")
        print(full_trans)
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
        intervals = [Interval(5, 15, ["A"]),
                     Interval(0, 10, ["B"])]

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
        print("==== Tets get old edges ===")
        old_edges = trans.get_old_edges(graph)
        self.assertTrue(3 in old_edges)
        self.assertTrue(4 in old_edges[3])
        print("== Old edges ==")
        print(old_edges)

        edges = old_edges
        trans.get_external_edges(graph, edges)
        print("=== External ===")
        print(edges)
        self.assertEqual(edges[0], [10  ])
        self.assertEqual(edges[10], [11])
        self.assertEqual(edges[11], [3])
        self.assertEqual(edges[3], [4])






if __name__ == "__main__":
    unittest.main()
