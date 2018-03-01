import unittest
import dummygraph
from offsetbasedgraph import GeneralMultiPathInterval, Position, Interval, Graph, Block, CriticalPathsMultiPathInterval


class TestInterval(unittest.TestCase):

    def _assertEqualsSinglePaths(self, multipath_interval,
                                 true_intervals):
        single_path_intervals = list(
            multipath_interval.get_single_path_intervals())
        self.assertEqual(len(single_path_intervals), len(true_intervals))
        for interval in single_path_intervals:
            self.assertTrue(interval in true_intervals)

    def _test_get_single_path_intervals_case(self):
        graph = Graph({2: Block(5), 3: Block(5)}, {2: [3]})
        multipath_interval = GeneralMultiPathInterval(
            [Position(2, 3)],
            [Position(3, 3)],
            [2, 3]
        )


    def test_get_single_path_intervals(self):
        graph = dummygraph.get_simple_graph()
        start_positions = [Position(1, 2)]
        end_positions = [Position(4, 2)]
        region_paths = [1, 2, 3, 4]

        # Test different middle
        multipath_interval = GeneralMultiPathInterval(
            start_positions,
            end_positions,
            region_paths,
            graph)
        true_intervals = [
            Interval(start_positions[0], end_positions[0], [1, 2, 4]),
            Interval(start_positions[0], end_positions[0], [1, 3, 4])
            ]
        self._assertEqualsSinglePaths(multipath_interval, true_intervals)

        # Test differnt end
        multipath_interval = GeneralMultiPathInterval(
            [Position(2, 1), Position(3, 1)],
            end_positions,
            [2, 3, 4],
            graph)
        true_intervals = [
            Interval(Position(2, 1), end_positions[0], [2, 4]),
            Interval(Position(3, 1), end_positions[0], [3, 4])
            ]
        self._assertEqualsSinglePaths(multipath_interval, true_intervals)

        # Test different start
        multipath_interval = GeneralMultiPathInterval(
            start_positions,
            [Position(2, 1), Position(3, 1)],
            [1, 2, 3],
            graph)
        true_intervals = [
            Interval(start_positions[0], Position(2, 1), [1, 2]),
            Interval(start_positions[0], Position(3, 1), [1, 3])
            ]
        self._assertEqualsSinglePaths(multipath_interval, true_intervals)

    def test_critical_path(self):

        graph = Graph(
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

        start = Position(3, 1)
        end = Position(3, 3)
        critical = [Interval(1, 2, [4], graph)]

        mp1 = CriticalPathsMultiPathInterval(start, end, critical)
        mp2 = CriticalPathsMultiPathInterval(start, end, critical)
        mp3 = CriticalPathsMultiPathInterval(Position(4, 1), end, critical)

        self.assertTrue(mp1.equal_critical_intervals(mp2))
        self.assertEqual(mp1, mp2)
        self.assertTrue(mp1 != mp3)
        self.assertTrue(mp1.equal_critical_intervals(mp3))



if __name__ == "__main__":
    unittest.main()
