from collections import defaultdict, namedtuple
from offsetbasedgraph.vcfgraph import *
from offsetbasedgraph.vcfmap import get_vcf_entries
from pyfaidx import Fasta
import pytest
import numpy as np


data_path = "/home/knut/Documents/phd/two_step_graph_mapper/benchmarking/mhc_graph_data/"
path = data_path + "1000genomes_variants.vcf"
fasta_path = data_path + "linear_ref.fa"


# construct_graph(get_vcf_entries(path))

SeqDummy = namedtuple("SeqDummy", ["seq"])


class FastaDummy:
    def __init__(self, seq):
        self._seq = seq

    def __getitem__(self, key):
        return SeqDummy(self._seq[key])


def test_simple_insertion():
    graph, _ = graph_from_snps_and_indels(
        np.array([[], []]), np.array([[2], [3]]), [], 10)
    true_graph = VCFGraph([2, 3, 8], AdjList([0, 2, 3, 3], [1, 2, 2]), SNPs())
    assert graph == true_graph


def test_two_insertion():
    graph, _ = graph_from_snps_and_indels(
        np.array([[], []]), np.array([[2, 4], [3, 2]]), [], 10)
    adj_list = AdjList.from_dict({0: [1, 2], 1: [2], 2: [3, 4], 3: [4]}, 5)
    true_graph = VCFGraph([2, 3, 2, 2, 6], adj_list, SNPs())
    assert graph == true_graph


def test_simple_deletion():
    graph, _ = graph_from_snps_and_indels(
        np.array([[3], [2]]), np.array([[], []]), [],  10)
    adj_list = AdjList.from_dict({0: [1, 2], 1: [2]}, 3)
    true_graph = VCFGraph([3, 2, 5], adj_list, SNPs())
    assert graph == true_graph


def test_simple_snp():
    graph, _ = graph_from_snps_and_indels(
        np.array([[], []]), np.array([[], []]), [2, 3, 5],  10)
    adj_list = AdjList.from_dict({0: []}, 1)
    true_graph = VCFGraph([10], adj_list, SNPs([0], [2, 3, 5]))
    assert graph == true_graph


@pytest.mark.skip
def test_construct_graph():
    fasta = Fasta(fasta_path)["6"]
    graph, ref, _ = construct_graph(get_vcf_entries(path), 4970557, fasta)
    graph.save("../vcfgraph_scripts/test_graph")
    ref.save("../vcfgraph_scripts/test_ref")


@pytest.mark.skip
def test_load_graph():
    VCFGraph.load("test_graph")


def test_build_seq_graph():
    insertion_seqs = ["AA", "TT", "GG"]
    insertion_node_ids = [1, 3, 5]
    reference_path = Path([0, 2, 4, 6], [0, 3, 6, 9, 12])
    fasta = FastaDummy("CCCGGGAAATTT")
    seq_graph = build_seq_graph(insertion_seqs, insertion_node_ids,
                                reference_path, fasta, 7)
    true_seqs = np.array(["CCC", "AA", "GGG", "TT", "AAA", "GG", "TTT"])
    assert np.all(seq_graph == true_seqs)


def test_build_full_graph():
    insertion_seqs = ["AA", "TT", "GG"]
    insertions = [[3, 6, 9], [2, 2, 2]]
    fasta = FastaDummy("CCCGGGAAATTT")
    graph, _ = graph_from_snps_and_indels([[], []], insertions, [], 12, insertion_seqs, [], fasta)
    adj_list = AdjList.from_dict({0: [1, 2], 1: [2], 2: [3, 4],
                                  3: [4], 4: [5, 6], 5: [6]}, 7)
    true_seqs = np.array(["CCC","AA", "GGG", "TT", "AAA", "GG", "TTT"])
    true_graph = VCFGraph([3, 2, 3, 2, 3, 2, 3], adj_list, SNPs(), seqs=true_seqs)
    assert graph == true_graph


    # assert np.all(seq_graph == np.array(["CCC","AA", "GGG", "TT", "AAA", "GG", "TTT"]))



if __name__ == "__main__":
    test_construct_graph()
