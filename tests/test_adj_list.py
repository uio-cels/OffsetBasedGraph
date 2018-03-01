import unittest
from offsetbasedgraph import GraphWithReversals
from offsetbasedgraph.graph import AdjListAsNumpyArrays
import numpy as np

class TestAdjListWithNumpyArrays(unittest.TestCase):
    def test(self):
        edges = {1: [2],
                 2: [3, 4],
                 4: [5],
                 3: [5],
                 10: [5, 3, -1],
                 100: [-20, -4],
                 -49: [4, 2, -1]
                 }

        npedges = AdjListAsNumpyArrays.create_from_edge_dict(edges)
        for node, edges in edges.items():
            self.assertTrue(np.all(npedges[node] == edges))

if __name__ == "__main__":
    unittest.main()