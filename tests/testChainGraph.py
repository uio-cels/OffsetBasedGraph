from __future__ import print_function
import unittest
import numpy as np
from GraphGenerators import *

from graph.OffsetBasedGraph import OffsetBasedGraph
FP = ""
class TestChainGraph(unittest.TestCase):

    def test_create_hg38(self):
        chrm_sizes = {"chr1": 10000,
                      "chr1_K1v1_alt": 1000,
                      "chr1_K2v1_alt": 2000}

        alt_loci = [{"chrom": "chr1", "chromStart": 1000, "chromEnd": 2000, "name": "chr1_K1v1_alt", "length": 1000},\
                    {"chrom": "chr1", "chromStart": 3000, "chromEnd": 4000, "name": "chr1_K2v1_alt", "length": 1000}]
        graph = OffsetBasedGraph.create_graph(chrm_sizes, alt_loci, True) #create_hg38_graph(chrm_sizes, alt_loci)

        print("Graph blocks")
        print([str(b) for b in graph.blocks])
        #print graph.block_edges
        self.assertIsNotNone(graph)

        self.assertEqual(len(graph.blocks), 7)
        self.assertEqual(sum([len(e) for e in graph.block_edges.values()]), 8)

        first_block = graph.get_block(
            LinearReference("hg38", "chr1", 0, 1, "+"))
        neighbour_ids = graph.block_edges[first_block.id]
        self.assertEqual(len(neighbour_ids), 2)
        end_block = graph.get_block(
            LinearReference("hg38", "chr1", 9998, 9999, "+"))
        m_block = graph.get_block(
            LinearReference("hg38", "chr1", 2500, 2600, "+"))
        for n_id in neighbour_ids:
            self.assertEqual(len(graph.block_edges[n_id]), 1)
            self.assertEqual(graph.block_edges[n_id][0], m_block.id)

    def test_create_hg38_2(self):
        chrm_sizes = {"chr1": 10000,
                      "chr1_K1v1_alt": 1000,
                      "chr1_K2v1_alt": 2000}

        alt_loci = [{"chrom": "chr1", "chromStart": 1000, "chromEnd": 2000, "name": "chr1_K1v1_alt", "length": 1000},\
                    {"chrom": "chr1", "chromStart": 1000, "chromEnd": 2000, "name": "chr1_K2v1_alt", "length": 1000}]
        graph = OffsetBasedGraph.create_graph(chrm_sizes, alt_loci, True) #create_hg38_graph(chrm_sizes, alt_loci)

        print([str(b) for b in graph.blocks])
        print(graph.block_edges)
        self.assertIsNotNone(graph)
        self.assertEqual(len(graph.blocks), 5)
        self.assertEqual(sum([len(e) for e in graph.block_edges.values()]), 6)

    def test_get_main_path_line_references(self):
        chrm_sizes = {"chr1": 1000000000,
                      "chr1_KI270762v1_alt": 1000000}

        alt_loci = [{"chrom": "chr1", "chromStart": 2448810, "chromEnd": 2791270, "name": "chr1_KI270762v1_alt", "length": 1000000}]
        graph = OffsetBasedGraph.create_graph(chrm_sizes, alt_loci) #create_hg38_graph(chrm_sizes, alt_loci)
        print("Graph blocks")
        print(graph.blocks.keys())
        self.assertEqual(graph.blocks["chr1-1"].main_path_linear_reference.chromosome, "chr1")
        self.assertEqual(graph.blocks["chr1-1"].main_path_linear_reference.start, 0)
        self.assertEqual(graph.blocks["chr1-1"].main_path_linear_reference.end, 2448810)

        self.assertEqual(graph.blocks["chr1_KI270762v1_alt-0"].main_path_linear_reference.chromosome, "chr1")
        self.assertEqual(graph.blocks["chr1_KI270762v1_alt-0"].main_path_linear_reference.start, 2448810)
        self.assertEqual(graph.blocks["chr1_KI270762v1_alt-0"].main_path_linear_reference.end, 2791270)


if __name__ == "__main__":
    unittest.main()
