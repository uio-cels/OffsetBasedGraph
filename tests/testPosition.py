import unittest
from offsetbasedgraph import Interval, Position, Graph, Translation, Block

class testPosition(unittest.TestCase):

    def testPositionEquality(self):
        pos1 = Position(1, 3)
        pos2 = Position(1, 4)
        pos3 = Position(1, 3)

        self.assertEqual(pos1, pos3)
        self.assertTrue(pos1 != pos2)


if __name__ == "__main__":
    unittest.main()

