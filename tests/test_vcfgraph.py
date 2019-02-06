from offsetbasedgraph.vcfgraph import construct_graph, VCFGraph, AdjList, SNPs, graph_from_indels
from offsetbasedgraph.vcfmap import get_vcf_entries
import pytest
import numpy as np

path = "/home/knut/Documents/phd/two_step_graph_mapper/benchmarking/mhc_graph_data/1000genomes_variants.vcf"


# construct_graph(get_vcf_entries(path))


def test_simple_insertion():
    graph = graph_from_indels(np.array([[], []]), np.array([[2], [3]]), 10)
    true_graph = VCFGraph([3, 3, 7], AdjList([0, 2, 3, 3], [1, 2, 2]), SNPs())
    assert graph == true_graph

