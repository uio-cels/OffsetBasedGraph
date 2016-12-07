import unittest
import dummygraph
from offsetbasedgraph import Graph, Block, Interval, Position, Translation

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

    def assert_graph_equals(self, graph, blocks, adj_list):
        if DEBUG:
            print(blocks)
            print(graph.blocks)
            print(adj_list)
            print(graph.adj_list)
        self.assertEqual(graph.blocks, blocks)
        self.assertEqual(graph.adj_list, adj_list)
        self.assertEqual(graph.reverse_adj_list,
                         graph._get_reverse_edges(adj_list))

    def test_init(self):
        blocks, adj_list = simple_graph()
        graph = Graph(blocks, adj_list)
        self.assert_graph_equals(graph, blocks, adj_list)

    def test_reverse_edges(self):
        adj_list = {1: [10, 20],
                    2: [100],
                    3: [100]}

        reverse_list = Graph._get_reverse_edges(adj_list)
        fasit = {10: [1], 20: [1], 100: [2, 3]}
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

    def test_split(self):
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

        print(true_trans)
        print(trans)
        self.assertEqual(trans, true_trans)

    def test_join(self):
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

    def _test_merge_translation(self):
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

    def test_get_all_block_borders(self):
        blocks = {1: Block(20), 2: Block(20),
                  11: Block(20), 12: Block(20)}

        adj_list = {1: [2], 11: [12]}
        graph = Graph(blocks, adj_list)

        interval_a = Interval(Position(1, 10), Position(2, 11))
        interval_b = Interval(Position(11, 15), Position(12, 16))
        border_list = graph._get_all_block_borders(interval_a, interval_b)
        self.assertEqual(border_list, [5, 10, 21])

    def test_split_blocks_at_starts_and_ends(self):
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

    def test_connect_intervals(self):
        pass

    def test_from_file(self):
        pass


if __name__ == "__main__":
    unittest.main()
