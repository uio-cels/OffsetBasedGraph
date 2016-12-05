import unittest
import dummygraph
from offsetbasedgraph import GeneralMultiPathInterval, Position, Interval


class TestInterval(unittest.TestCase):

    def _assertEqualsSinglePaths(self, multipath_interval,
                                 true_intervals):
        single_path_intervals = list(
            multipath_interval.get_single_path_intervals())
        self.assertEqual(len(single_path_intervals), len(true_intervals))
        for interval in single_path_intervals:
            self.assertTrue(interval in true_intervals)

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


if __name__ == "__main__":
    unittest.main()
