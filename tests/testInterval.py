import unittest
from . import dummygraph
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
        

if __name__ == "__main__":
    unittest.main()
