import unittest
import numpy as np
from graph.OffsetBasedGraph import OffsetBasedGraph
from segments.LinearInterval import LinearInterval
from segments.Interval import Interval

class TestInterval(unittest.TestCase):

    def setUp(self):
        chrm_sizes = {"chr1": 10000,
                      "chr1_K1v1_alt": 1000,
                      "chr1_K2v1_alt": 2000}

        alt_loci = [{"chrom": "chr1", "chromStart": 1000, "chromEnd": 2000, "name": "chr1_K1v1_alt", "length": 1000},\
                    {"chrom": "chr1", "chromStart": 3000, "chromEnd": 4000, "name": "chr1_K2v1_alt", "length": 1000}]
        self.graph = OffsetBasedGraph.create_graph(chrm_sizes, alt_loci)

    def testCreateFromLinearInterval(self):
        interval = Interval(self.graph)
        interval.create_from_linear_interval("chr1", 500, 1500)
        self.assertEqual(len(interval.block_list), 2, "Block list should have length 2")

        interval = Interval(self.graph)
        interval.create_from_linear_interval("chr1", 500, 2500)
        self.assertEqual(len(interval.block_list), 3, "Block list should have length 3")

        interval = Interval(self.graph)
        interval.create_from_linear_interval("chr1", 500, 4500)
        self.assertEqual(len(interval.block_list), 5, "Block list should have length 5 ")



if __name__ == "__main__":
    unittest.main()

