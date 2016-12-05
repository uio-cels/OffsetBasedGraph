import unittest
import dummygraph
from offsetbasedgraph import Interval, Position


class TestInterval(unittest.TestCase):

    def testSimpelInterval(self):
        region_paths = [1, 3, 4]
        interval = Interval(Position(1, 10), Position(4, 10), region_paths)

        for r in region_paths:
            self.assertTrue(
                r in interval.region_paths,
                "The region path %d is not in interval's region paths")

        self.assertEqual(len(interval.region_paths), 3,
                         "Interval should have 3 region paths")

    def test_interval_length(self):
        graph = dummygraph.get_simple_graph()
        interval = Interval(Position(1, 5),
                            Position(4, 10),
                            [1, 2, 4],
                            graph=graph)
        true_length = 5 + 20 + 10
        self.assertEqual(interval.length(), true_length)

    def test_get_position_from_offset(self):
        graph = dummygraph.get_simple_graph()
        interval = Interval(Position(1, 5),
                            Position(4, 10),
                            [1, 2, 4],
                            graph=graph)
        offsets = [3, 23, 33]
        positions = [Position(1, 8),
                     Position(2, 18),
                     Position(4, 8)]

        for offset, position in zip(offsets, positions):
            self.assertEqual(interval.get_position_from_offset(offset),
                             position)

    def test_split(self):
        graph = dummygraph.get_simple_graph()
        interval = Interval(Position(1, 5),
                            Position(4, 10),
                            [1, 2, 4],
                            graph=graph)
        splits = interval.split([2, 7, 27])
        true_intervals = [
            Interval(Position(1, 5), Position(1, 7)),
            Interval(Position(1, 7), Position(2, 2), [1, 2]),
            Interval(Position(2, 2), Position(4, 2), [2, 4]),
            Interval(Position(4, 2), Position(4, 10))
            ]
        self.assertEqual(splits, true_intervals)

    def test_join(self):
        graph = dummygraph.get_simple_graph()
        interval = Interval(Position(1, 5),
                            Position(4, 10),
                            [1, 2, 4],
                            graph=graph)
        splits = interval.split([7])
        self.assertEqual(splits[0].join(splits[1]),
                         interval)


if __name__ == "__main__":
    unittest.main()
