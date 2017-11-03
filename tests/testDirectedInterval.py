import unittest
from offsetbasedgraph import Graph, Block, Interval, DirectedInterval, Position, GraphWithReversals

class TestDirectedInterval(unittest.TestCase):

    def test_correct_subclassing(self):
        directed_interval = DirectedInterval(1, 5, [1, 2, 3])
        interval = Interval(1, 5, [1, 2, 3])
        self.assertEqual(interval, directed_interval)

    def test_contains_rp(self):
        rps = [1, 2, -3, 4, -5, 6, -7]
        interval = DirectedInterval(1, 2, rps)
        for rp in rps:
            self.assertTrue(interval.contains_rp(rp))
            self.assertTrue(interval.contains_rp(-rp))

    def test_contains_position(self):

        # Case 1
        graph = GraphWithReversals({1: Block(10)}, {})
        interval = DirectedInterval(1, 5, [1], graph)
        interval_identical_reversed = DirectedInterval(6, 10, [-1], graph)

        for i in range(1, 5):
            self.assertTrue(interval.contains_position(Position(1, i)))
            self.assertTrue(interval_identical_reversed.contains_position(Position(1, i)))
            self.assertFalse(interval.contains_position(Position(-1, i)), \
                             "Position %d,%d should not be in interval" % (-1, i))
            self.assertFalse(interval_identical_reversed.contains_position(Position(-1, i)), \
                             "Position %d,%d should not be in interval" % (-1, i))

        self.assertFalse(interval.contains_position(Position(1, 5)))

        # Case 2
        graph = GraphWithReversals({1: Block(10), 2: Block(10), 3: Block(10)}, {1: [2], 2: [3]})
        interval = DirectedInterval(Position(1, 5), Position(-3, 5), [1, 2, -3], graph)
        self.assertTrue(interval.contains_position(Position(1, 5)))
        self.assertFalse(interval.contains_position(Position(-1, 9)))
        self.assertTrue(interval.contains_position(Position(2, 5)))
        self.assertTrue(interval.contains_position(Position(3, 9)))
        self.assertTrue(interval.contains_position(Position(3, 8)))
        self.assertTrue(interval.contains_position(Position(3, 7)))
        self.assertTrue(interval.contains_position(Position(3, 6)))
        self.assertFalse(interval.contains_position(Position(3, 5)))
        self.assertFalse(interval.contains_position(Position(3, 4)))
        self.assertTrue(interval.contains_position(Position(-3, 3)))
        self.assertFalse(interval.contains_position(Position(-3, 5)))

    def test_can_be_on_strand(self):

        graph = GraphWithReversals(
            {1: Block(3),
             2: Block(3),
             3: Block(3)},
            {1: [2],
             2: [3]}
        )

        interval = DirectedInterval(1, 3, [2], graph)
        self.assertTrue(interval.can_be_on_negative_strand())
        self.assertTrue(interval.can_be_on_positive_strand())

        interval = DirectedInterval(1, 1, [1, 2, 3], graph)
        self.assertTrue(interval.can_be_on_positive_strand())
        self.assertFalse(interval.can_be_on_negative_strand())

        interval = DirectedInterval(1, 1, [-3, -2, -1], graph)
        self.assertFalse(interval.can_be_on_positive_strand())
        self.assertTrue(interval.can_be_on_negative_strand())



        graph = GraphWithReversals(
            {1: Block(3),
             2: Block(3)},
            {1: [2],
             -2: [-1]}
        )
        interval = DirectedInterval(1, 1, [1, 2], graph)
        self.assertTrue(interval.can_be_on_positive_strand())
        self.assertTrue(interval.can_be_on_negative_strand())





if __name__ == "__main__":
    unittest.main()
