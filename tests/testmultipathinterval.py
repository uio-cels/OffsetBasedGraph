import unittest
import dummygraph
from offsetbasedgraph import GeneralMultiPathInterval, Position, Interval


class TestInterval(unittest.TestCase):

    def test_get_single_path_intervals(self):
        graph = dummygraph.get_simple_graph()
        start_positions = [Position(1, 2)]
        end_positions = [Position(4, 2)]
        region_paths = [1, 2, 3, 4]
        multipath_interval = GeneralMultiPathInterval(
            start_positions,
            end_positions,
            region_paths,
            graph)

        single_path_intervals = multipath_interval.get_single_path_intervals()
        true_intervals = [
            Interval(start_positions[0], end_positions[0], [1, 2, 3]),
            Interval(start_positions[0], end_positions[0], [1, 2, 4])
            ]
        self.assertEqual(len(single_path_intervals), len(true_intervals))
        for interval in single_path_intervals:
            self.assertTrue(interval in true_intervals)

if __name__ == "__main__":
    unittest.main()
