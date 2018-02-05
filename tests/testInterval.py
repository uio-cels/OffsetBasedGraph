import unittest
import dummygraph
from offsetbasedgraph import Interval, Position, IntervalCollection, Graph, Block, IndexedInterval


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

    def _test_split(self):
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

    def _test_join(self):
        graph = dummygraph.get_simple_graph()
        interval = Interval(Position(1, 5),
                            Position(4, 10),
                            [1, 2, 4],
                            graph=graph)
        splits = interval.split([7])
        self.assertEqual(splits[0].join(splits[1]),
                         interval)

    def test_hash(self):
        interval1 = Interval(5, 10, [1, 2, 3, 4])
        interval2 = Interval(5, 10, [1, 2, 3, 4])
        interval_different = Interval(4, 10, [1, 2, 3, 4])
        interval_different2 = Interval(5, 11, [1, 2, 3, 4])
        interval_different3 = Interval(5, 10, [1, 2, 3, 5])

        interval1_minus = Interval(5, 10, [1, 2, 3, 4], direction=-1)

        self.assertEqual(interval2.hash(), interval2.hash())
        self.assertTrue(interval1.hash() != interval_different.hash())
        self.assertTrue(interval1.hash() != interval_different2.hash())
        self.assertTrue(interval1.hash() != interval_different3.hash())
        self.assertTrue(interval1.hash() != interval1_minus.hash())


    def test_overlap(self):
        graph = Graph(
            {
                1: Block(10),
                2: Block(10),
                3: Block(10)
            },
            {
                1: [2],
                2: [3],
                3: [1]
            }
        )

        interval1 = Interval(0, 10, [1, 2], graph)
        interval2 = Interval(0, 10, [1, 2], graph)
        self.assertEqual(interval1.overlap(interval2), 20)

        interval1 = Interval(5, 10, [1, 2], graph)
        interval2 = Interval(0, 10, [1, 2], graph)
        self.assertEqual(interval1.overlap(interval2), 15)

        interval1 = Interval(5, 7, [1], graph)
        interval2 = Interval(0, 10, [1, 2], graph)
        self.assertEqual(interval1.overlap(interval2), 2)

        interval1 = Interval(5, 7, [1], graph)
        interval2 = Interval(6, 7, [1], graph)
        self.assertEqual(interval1.overlap(interval2), 1)

        interval1 = Interval(8, 2, [1, 2, 3, 1], graph)
        interval2 = Interval(0, 10, [1], graph)
        self.assertEqual(interval1.overlap(interval2), 4)

        interval1 = Interval(8, 2, [1, 2, 3, 1], graph)
        interval2 = Interval(0, 2, [1, 2], graph)
        self.assertEqual(interval1.overlap(interval2), 6)

    def test_contains_correct_order(self):

        interval = Interval(0, 10, [1, 2, 3, 4, 5, 6, 7])

        other = Interval(0, 10, [1, 2, 3])
        self.assertTrue(interval.contains_in_correct_order(other))

        other = Interval(0, 10, [1, 2, 4])
        self.assertFalse(interval.contains_in_correct_order(other))

        other = Interval(0, 10, [1, 2, 5])
        self.assertFalse(interval.contains_in_correct_order(other))

        other = Interval(0, 10, [4, 5, 6])
        self.assertTrue(interval.contains_in_correct_order(other))


class TestIndexedInterval(unittest.TestCase):
    def test_position_at_offset(self):
        graph = Graph(
            {
                1: Block(10),
                2: Block(10),
                3: Block(10)
            },
            {
                1: [2],
                2: [3]
            }
        )
        interval = IndexedInterval(4, 6, [1, 2, 3], graph)

        self.assertEqual(interval.position_at_offset(0), Position(1, 4))
        self.assertEqual(interval.position_at_offset(1), Position(1, 5))
        self.assertEqual(interval.position_at_offset(2), Position(1, 6))
        self.assertEqual(interval.position_at_offset(5), Position(1, 9))
        self.assertEqual(interval.position_at_offset(6), Position(2, 0))
        self.assertEqual(interval.position_at_offset(7), Position(2, 1))
        self.assertEqual(interval.position_at_offset(16), Position(3, 0))
        self.assertEqual(interval.position_at_offset(21), Position(3, 5))

    def test_get_offset_at_position(self):
        graph = Graph(
            {
                1: Block(10),
                2: Block(10),
                3: Block(10)
            },
            {
                1: [2],
                2: [3]
            }
        )
        interval = IndexedInterval(4, 6, [1, 2, 3], graph)
        self.assertEqual(interval.get_offset_at_position(Position(1, 4)), 0)
        self.assertEqual(interval.get_offset_at_position(Position(2, 0)), 6)
        self.assertEqual(interval.get_offset_at_position(Position(2, 1)), 7)
        self.assertEqual(interval.get_offset_at_position(Position(2, 9)), 15)
        self.assertEqual(interval.get_offset_at_position(Position(3, 0)), 16)

        self.assertEqual(interval.get_reverse_offset_at_position(Position(-1, 6)), 0)
        self.assertEqual(interval.get_reverse_offset_at_position(Position(-1, 5), is_end_pos=False), 0)
        self.assertEqual(interval.get_reverse_offset_at_position(Position(-1, 5)), 1)
        self.assertEqual(interval.get_reverse_offset_at_position(Position(-3, 5)), 21)
        self.assertEqual(interval.get_reverse_offset_at_position(Position(-2, 0)), 16)

    def test_get_subinterval(self):
        graph = Graph(
            {
                1: Block(10),
                2: Block(10),
                3: Block(10)
            },
            {
                1: [2],
                2: [3]
            }
        )
        interval = IndexedInterval(4, 6, [1, 2, 3], graph)
        self.assertEqual(interval.get_subinterval(0, 1), Interval(4, 5, [1]))
        self.assertEqual(interval.get_subinterval(0, 10), Interval(4, 4, [1, 2]))
        self.assertEqual(interval.get_subinterval(0, 11), Interval(4, 5, [1, 2]))
        self.assertEqual(interval.get_subinterval(1, 20), Interval(5, 4, [1, 2, 3]))
        self.assertEqual(interval.get_subinterval(10, 20), Interval(4, 4, [2, 3]))
        self.assertEqual(interval.get_subinterval(10, 16), Interval(4, 10, [2]))



class TestIntervalCollection(unittest.TestCase):
    def test_to_file_from_file(self):
        intervals = (
                        Interval(0, 5, [1]),
                        Interval(0, 3, [2])
                    )
        collection = IntervalCollection(intervals)
        print(collection.intervals)
        collection.to_file("test_intervalcollection.tmp", text_file=True)

        collection2 = IntervalCollection.from_file("test_intervalcollection.tmp", text_file=True)

        #self.assertTrue(collection.intervals == collection2.intervals)
        for i, interval in enumerate(collection2):
            self.assertEqual(interval, intervals[i])

    def test_gzip(self):
        intervals = (
                        Interval(0, 5, [1]),
                        Interval(0, 1, [1]),
                        Interval(1, 4, [1,2]),
                        Interval(1, 4, [1,2]),
                        Interval(1, 4, [1,2]),
                        Interval(1, 4, [1,2]),
                        Interval(0, 3, [2])
                    )
        collection = IntervalCollection(intervals)
        print(collection.intervals)
        collection.to_gzip("test_intervalcollection.gzip")

        collection2 = IntervalCollection.from_gzip("test_intervalcollection.gzip")

        #self.assertTrue(collection.intervals == collection2.intervals)
        for i, interval in enumerate(collection2):
            self.assertEqual(interval, intervals[i])


if __name__ == "__main__":
    unittest.main()
