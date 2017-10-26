import unittest
import numpy as np
from offsetbasedgraph.distancematrix import DistanceIndex
from offsetbasedgraph import GraphWithReversals, Block
from itertools import chain

class TestDistanceIndex(unittest.TestCase):
    def setUp(self):
        nodes = {i: Block(10) for i in range(1, 11)}
        edges = {i: [i+1] for i in range(1, 10)}
        edges[5] += [10]
        self.graph = GraphWithReversals(nodes, edges)
        self._n_nodes = len(self.graph.blocks.values())
        self._create_true_covered_list()
        self._create_true_partial_list()
        # self.maxDiff = None

    def _create_true_covered_list(self):
        covered_neighbours = {i: [i] for i in range(1, 11)}
        for i in range(1, 10):
            covered_neighbours[i].append((i+1))
        for i in range(2, 11):
            covered_neighbours[i].append(i-1)
        covered_neighbours[5].extend([9, 10])
        covered_neighbours[6].extend([9, 10])
        covered_neighbours[10].extend([5, 6])
        covered_neighbours[9].extend([5, 6])
        self.covered_neighbours = {
            i: list(sorted(l)) for i, l in covered_neighbours.items()}

    def _create_true_partial_list(self):
        partial_neighbours = {
            i: [] for i in chain(range(1, 11), range(-10, 0))}
        for i in range(1, 9):
            partial_neighbours[i].append((i+2, 5))
        for i in range(3, 11):
            partial_neighbours[-i].append((-i+2, 5))
        partial_neighbours[5].append((-8, 5))
        partial_neighbours[-6].append((-8, 5))

        partial_neighbours[4].extend([(-9, 5), (10, 5)])
        partial_neighbours[-7].extend([(-9, 5), (10, 5)])

        partial_neighbours[8].extend([(-5, 5), (6, 5)])
        partial_neighbours[9].extend([(-4, 5), (7, 5)])
        partial_neighbours[-10].extend([(-4, 5), (7, 5)])
        self.partial_neighbours = {i: list(sorted(l, key=lambda x: x[0]))
                                   for i, l in partial_neighbours.items()}

    def test_create(self):
        distance_index = DistanceIndex(self.graph, 25)
        distance_index.create()
        self.assertEqual(distance_index.covered_neighbours,
                         self.covered_neighbours)
        for node_id in distance_index.partial_neighbours:
            print(node_id,
                  distance_index.partial_neighbours[node_id],
                  self.partial_neighbours[node_id])
        self.assertEqual(distance_index.partial_neighbours,
                         self.partial_neighbours)


if __name__ == "__main__":
    unittest.main()
