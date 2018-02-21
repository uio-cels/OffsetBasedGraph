import unittest
import dummygraph
from offsetbasedgraph import Graph, Block, Interval, Position, Translation, BlockCollection


DEBUG = False


def simple_graph():
    blocks = {_id: block for _id, block in
              enumerate([10, 20, 30, 40])}
    adj_list = {0: [1], 1: [2], 2: [3], 0: [4], 4: [3]}
    return blocks, adj_list


def disjoint_graph():
    blocks = {_id: Block(block) for _id, block in
              enumerate([30, 40])}
    adj_list = {}
    return Graph(blocks, adj_list)


class TestGraph(unittest.TestCase):

    def test_init(self):
        blocks = {1: Block(5), 2: Block(3)}
        graph = Graph(blocks, {})

    def test_block_in_graph(self):
        blocks = {1: Block(3), 2: Block(3)}
        graph = Graph(blocks, {})

        self.assertTrue(graph.block_in_graph(1))
        self.assertTrue(graph.block_in_graph(2))
        self.assertTrue(not graph.block_in_graph(3))


    def assert_graph_equals(self, graph, blocks, adj_list):
        new_graph = Graph(blocks, adj_list)
        self.assertEqual(new_graph, graph)

    def test_graph_equals(self):
        # Case 1
        graph1 = Graph({1: Block(10)}, {})
        graph2 = Graph({1: Block(10)}, {})
        self.assertEqual(graph1, graph2)

        # Case 2
        graph1 = Graph({1: Block(10), 2: Block(1)}, {})
        graph2 = Graph({1: Block(10)}, {})
        self.assertTrue(graph1 != graph2)

        # Case 3
        graph1 = Graph({1: Block(10), 2: Block(1)}, {1: []})
        graph2 = Graph({1: Block(10), 2: Block(1)}, {})
        self.assertEqual(graph1, graph2)

        # Case 4
        graph1 = Graph({1: Block(10), 2: Block(1)}, {1: [2]})
        graph2 = Graph({1: Block(10), 2: Block(1)}, {1: [2]})
        self.assertEqual(graph1, graph2)

    def test_init(self):
        blocks, adj_list = simple_graph()
        graph = Graph(blocks, adj_list)
        self.assert_graph_equals(graph, blocks, adj_list)

    def test_reverse_edges(self):
        adj_list = {1: [10, 20],
                    2: [100],
                    3: [100]}

        reverse_list = Graph._get_reverse_edges(adj_list)
        fasit = {-10: [-1], -20: [-1], -100: [-2, -3]}
        self.assertEqual(reverse_list, fasit)

    def _setup_merge(self):
        graph = disjoint_graph()
        interval_a = Interval(
            Position(0, 0),
            Position(0, 10), graph=graph)

        interval_b = Interval(
            Position(1, 0),
            Position(1, 10), graph=graph)

        return graph, interval_a, interval_b

    def _test_split(self):
        graph = dummygraph.get_simple_graph()
        splits = [5, 12]
        trans = graph._split_block(2, splits)
        
        # Test graph structure
        new_id = [b for b in graph.adj_list[1] if not b == 3][0]
        new_id2 = [b for b in graph.adj_list[new_id]][0]
        new_id3 = [b for b in graph.adj_list[new_id2]][0]

        new_blocks = {1: Block(10), new_id: Block(5),
                      new_id2: Block(7), new_id3: Block(8),
                      3: Block(10), 4: Block(15)}

        new_adj_list = {1: [3, new_id], 3: [4],
                        new_id: [new_id2],
                        new_id2: [new_id3],
                        new_id3: [4]}

        self.assert_graph_equals(graph, new_blocks, new_adj_list)

        # Test translation
        new_interval = Interval(Position(new_id, 0),
                                Position(new_id3, 8),
                                [new_id, new_id2, new_id3])
        old_intervals = [Interval(Position(2, 0), Position(2, 5)),
                         Interval(Position(2, 5), Position(2, 12)),
                         Interval(Position(2, 12), Position(2, 20))]

        true_trans = Translation(
            {2: [new_interval]},
            {_id: [interval] for _id, interval in
             zip([new_id, new_id2, new_id3], old_intervals)}
            )

        self.assertEqual(trans, true_trans)

    def _test_join(self):
        graph = dummygraph.get_mergable_graph()
        trans = graph._join_blocks([2, 3])

        # Test graph structure
        new_id = graph.adj_list[1][0]
        blocks = {1: Block(10), new_id: Block(20), 4: Block(15)}
        block_edges = {1: [new_id], new_id: [4]}
        self.assert_graph_equals(graph, blocks, block_edges)

        # Test translation
        new_interval = Interval(Position(new_id, 0),
                                Position(new_id, 20))
        old_intervals = [Interval(Position(_id, 0),
                                  Position(_id, 20))
                         for _id in (2, 3)]
        true_translation = Translation(
            {2: [new_interval],
             3: [new_interval]},
            {new_id: old_intervals})

        self.assertEqual(
            true_translation,
            trans)

    def test_insulate_translation(self):
        graph, interval_a, interval_b = self._setup_merge()
        translation, _ = graph._get_inslulate_translation(
            [interval_a, interval_b])
        a_to_b = translation._a_to_b
        self.assertEqual(len(a_to_b), 2)
        i_first = a_to_b[0][0]
        new_ids = i_first.region_paths
        self.assertEqual(i_first, Interval(
            Position(new_ids[0], 0), Position(new_ids[1], 20)))
        i_last = a_to_b[1][0]
        new_ids = i_last.region_paths
        self.assertEqual(i_last, Interval(
            Position(new_ids[0], 0), Position(new_ids[1], 30)))

    def _test_insulated_translation(self):
        graph = dummygraph.get_insulated_graph()
        interval_a = Interval(Position(2, 0), Position(4, 10),
                              [2, 3, 4], graph=graph)
        interval_b = Interval(Position(5, 0), Position(7, 10),
                              [5, 6, 7], graph=graph)
        intervals = [interval_a, interval_b]
        trans, new_graph = graph._get_insulated_merge_transformation(intervals)
        # 
        # graph, interval_a, interval_b = self._setup_merge()
        # translation, _ = graph._get_inslulate_translation(
        #     [interval_a, interval_b])
        # a_to_b = translation._a_to_b
        # self.assertEqual(len(a_to_b), 2)
        # i_first = a_to_b[0][0]
        # new_ids = i_first.region_paths
        # self.assertEqual(i_first, Interval(
        #     Position(new_ids[0], 0), Position(new_ids[1], 20)))
        # i_last = a_to_b[1][0]
        # new_ids = i_last.region_paths
        # self.assertEqual(i_last, Interval(
        #     Position(new_ids[0], 0), Position(new_ids[1], 30)))

    def _test_merge_translation(self):
        # Not in use any more!
        graph, interval_a, interval_b = self._setup_merge()
        translation = graph.merge_intervals(
            interval_a,
            interval_b)
        A = 0
        B = 1
        a_first = translation._a_to_b[A][0].region_paths[0]
        a_last = translation._a_to_b[A][0].region_paths[1]
        b_first = translation._a_to_b[B][0].region_paths[0]
        b_last = translation._a_to_b[B][0].region_paths[1]

        # Sanity check on merge
        self.assertEqual(a_first, b_first)
        self.assertEqual(a_first.length() == interval_a.length())

        # Check that translation object is correct
        a_to_b = {
            A: Interval(Position(a_first, 0),
                        Position(a_last, 20)),
            B: Interval(Position(b_first, 0),
                        Position(b_last, 30))
            }
        b_to_a = {
            a_first: set(interval_a, interval_b),
            a_last: set(
                Interval(
                    Position(A, 10),
                    Position(A, 30)),
                Interval(
                    Position(B, 10),
                    Position(B, 40))
                )
            }
        self.assertEqual(
            Translation(a_to_b, b_to_a),
            translation)

        # Check graph
        blocks = {a_first: Block(10),
                  a_last: Block(20),
                  b_last: Block(30)
                  }
        adj_list = {a_first: [a_last,  b_last]}
        self.assert_graph_equals(
            graph, Graph(blocks, adj_list))

    def test_merge_translations2(self):
        graph = Graph({1: Block(3), 2: Block(2), 3: Block(2),
                       4: Block(3), 5: Block(1)},
                      {1: [2, 5], 3: [4]} )
        #graph = Graph({1: Block(3), 2: Block(2), 3: Block(2),
        #               4: Block(3)},
        #              {1: [2], 3: [4]})
        interval1 = Interval(1, 1, [1, 2], graph)
        interval2 = Interval(1, 2, [3, 4], graph)
        new_graph, trans = graph.merge([interval1, interval2])
        #print(trans)
        #print(new_graph)

        a = trans._a_to_b[1][0].region_paths
        b = trans._a_to_b[3][0].region_paths
        c = trans._a_to_b[2][0].region_paths
        d = trans._a_to_b[4][0].region_paths
        all_blocks = a+b+c+d+[5]
        blocks = {block_id: Block(1) for block_id in all_blocks}
        edges = {a[0]: [a[1]], b[0]: [b[1]], a[1]: [a[2]],
                 a[2]: [c[0], 5], c[0]: [c[1], d[2]]}
        self.assert_graph_equals(new_graph, blocks, edges)

    def test_connect_positions(self):
        graph = dummygraph.get_disjoint_graph()
        pos_a = Position(1, 4)
        pos_b = Position(2, 10)
        n_g, trans = graph.connect_postitions(pos_a, pos_b)
        a, b = trans.translate(Interval(0, 10, [1])).region_paths
        c, d = trans.translate(Interval(0, 20, [2])).region_paths

        blocks = {a: Block(5), b: Block(5), c: Block(10), d: Block(10),
                  3: Block(30)}

        adj_list = {a: [b, d], c: [d]}
        t_g = Graph(blocks, adj_list)
        self.assertEqual(n_g, t_g)

    def _test_get_all_block_borders(self):
        blocks = {1: Block(20), 2: Block(20),
                  11: Block(20), 12: Block(20)}

        adj_list = {1: [2], 11: [12]}
        graph = Graph(blocks, adj_list)

        interval_a = Interval(Position(1, 10), Position(2, 11))
        interval_b = Interval(Position(11, 15), Position(12, 16))
        border_list = graph._get_all_block_borders(interval_a, interval_b)
        self.assertEqual(border_list, [5, 10, 21])

    def _test_split_blocks_at_starts_and_ends(self):
        graph = dummygraph.get_disjoint_graph()
        intervals = [Interval(Position(1, 0), Position(1, 5), graph=graph),
                     Interval(Position(2, 5), Position(2, 10), graph=graph),
                     Interval(Position(3, 25), Position(3, 30), graph=graph)]

        trans = graph._split_blocks_at_starts_and_ends(intervals)

        rps1 = trans.translate_rp(1).region_paths
        rps2 = trans.translate_rp(2).region_paths
        rps3 = trans.translate_rp(3).region_paths

        # Check that starts and ends have not been split
        self.assertEqual(len(rps1), 2)
        self.assertEqual(len(rps2), 3)
        self.assertEqual(len(rps3), 2)

        true_blocks = dict(zip(rps1, [Block(5), Block(5)]))
        true_blocks.update(dict(zip(rps2, [Block(5), Block(5), Block(10)])))
        true_blocks.update(dict(zip(rps3, [Block(25), Block(5)])))

        true_adjs = {rps1[0]: [rps1[1]],
                     rps2[0]: [rps2[1]], rps2[1]: [rps2[2]],
                     rps3[0]: [rps3[1]]}
        self.assert_graph_equals(graph, true_blocks, true_adjs)

    def test_merge_start_block(self):
        """
        Special case where one interval starts at block
        :return:
        """
        # Real case from grch38
        # graph = Graph({0: Block(185285), 1: Block(248956422), 2: Block(182439)}, {})
        # intervals = [Interval(Position(1, 144488705), Position(1, 144555943), [1], graph),
        #             Interval(Position(0, 0), Position(0, 67238), [0], graph)]

        graph = Graph({0: Block(2), 1: Block(4)}, {})
        intervals = [Interval(Position(1, 1), Position(1, 2), [1], graph),
                     Interval(Position(0, 0), Position(0, 1), [0], graph)]

        new_graph, trans = graph.merge(intervals)

    def test_merge_end_block(self):
        graph = Graph({0: Block(2), 1: Block(4)}, {})
        intervals = [Interval(Position(1, 1), Position(1, 2), [1], graph),
                     Interval(Position(0, 1), Position(0, 2), [0], graph)]

        new_graph, trans = graph.merge(intervals)
        # print("test_merge_end_block_graph:")
        # print(new_graph)
        # print(trans)

    def test_merge_two_end_block2(self):
        # print("test end block 2")
        graph = Graph({8: Block(118047), 1: Block(182439), 3: Block(144488705), 12: Block(67238), 6: Block(104400479)},
                      {3: [12], 12: [8, 6]})

        intervals = [Interval(Position(6, 117843), Position(6, 118838), [6], graph),
                     Interval(Position(8, 117052), Position(8, 118047), [8], graph)]
        new_graph, trans = graph.merge(intervals)
        # print(new_graph)

    def test_connect_intervals(self):
        pass

    def test_from_file(self):
        pass

    def test_has_identical_structure(self):
        # Case 1
        g1 = Graph(
            {
                1: Block(1),
                2: Block(10)
            },
            {
                1: [2]
            }
        )

        g2 = Graph(
            {
                5: Block(10),
                2: Block(1)
            },
            {
                5: [2]
            }
        )

        self.assertTrue(g1.has_identical_structure(g2))

        # Case 2
        g1 = Graph(
            {
                1: Block(1),
                2: Block(1),
                3: Block(1),
                4: Block(1)
            },
            {
                1: [2, 3],
                3: [4],
                2: [4]
            }
        )

        g2 = Graph(
            {
                10: Block(2),
                20: Block(2),
                30: Block(2),
                40: Block(2)
            },
            {
                10: [20, 30],
                30: [40],
                20: [40]
            }
        )

        self.assertTrue(g1.has_identical_structure(g2))

    def test_find_critical_blocks(self):
        pass
        #graph = dummygraph.get_realistic_graph()
        #critical_blocks = graph.find_critical_blocks(0)
        #self.assertEqual(critical_blocks, [1, 5, 6])

    def _test_to_from_file(self):
        for graph in [dummygraph.get_simple_graph(),
                      dummygraph.get_disjoint_graph()]:
            graph.to_file("test_graph")

            graph2 = Graph.from_file("test_graph")
            self.assertEqual(graph, graph2)

    def test_get_arbitrary_linear_graph(self):
        initial_graph = Graph(
            {
                1: Block(10),
                2: Block(10),
                3: Block(10),
                4: Block(10)
            },
            {
                1: [2, 4],
                2: [3],
                4: [3]
            }
        )

        graph, trans = initial_graph.get_arbitrary_linear_graph()

        self.assertTrue(len(graph.blocks) == 1)
        self.assertTrue(len(graph.adj_list) == 0)
        self.assertTrue(list(graph.blocks.values())[0].length() == 30)

        # Case 2
        initial_graph = Graph(
            {
                1: Block(10),
                2: Block(10),
                3: Block(10),
                4: Block(10),
                5: Block(10),
                6: Block(5)
            },
            {
                1: [2, 4],
                2: [3],
                4: [3],
                5: [6]
            }
        )


        graph, trans = initial_graph.get_arbitrary_linear_graph()

        self.assertTrue(len(graph.blocks) == 2)
        self.assertTrue(len(graph.adj_list) == 0)
        sum_length = sum([b.length ()for b in graph.blocks.values()])
        self.assertEqual(sum_length, 45)

    def test_block_collection(self):
        blocks = BlockCollection({1: Block(10), 2: Block(12)})

        self.assertTrue(1 in blocks)
        self.assertTrue(2 in blocks)
        self.assertTrue(3 not in blocks)

if __name__ == "__main__":
    unittest.main()
