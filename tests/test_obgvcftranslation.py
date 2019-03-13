import offsetbasedgraph as obg
from offsetbasedgraph.vcfgraph import *
from offsetbasedgraph.obg_vcf_translation import *
import pytest


data_path = "/home/knut/Documents/phd/two_step_graph_mapper/benchmarking/mhc_graph_data/"
path = data_path + "1000genomes_variants.vcf"
fasta_path = data_path + "linear_ref.fa"


@pytest.fixture
def simple_graph():
    lens = [4, 1, 3, 2, 2, 3, 1]
    nodes = dict(enumerate((obg.Block(l) for l in lens), 1))
    adj_list = {i: [i+1] for i in range(1, 7)}
    graph = obg.Graph(nodes, adj_list)
    graph.convert_to_numpy_backend()
    return graph


@pytest.fixture
def ins_graph():
    adj_list = {1: [2, 3], 2: [3], 3: [4, 5], 4: [5]}
    nodes = {i: obg.Block(3) for i in range(1, 6)}
    graph = obg.Graph(nodes, adj_list)
    graph.convert_to_numpy_backend()
    return graph


@pytest.fixture
def vcf_ins_graph():
    adj_list = {0: [1, 2], 1: [2], 2: [3, 4], 3: [4]}
    adj_list = AdjList(adj_list, 5)
    node_lens = [3]*5
    seqs = ["AAT", "TGG", "CCA", "ATT", "GGC"]
    seqs = [s.lower() for s in seqs]
    return VCFGraph(node_lens, adj_list, SNPs(), seqs)


@pytest.fixture
def simple_vcf_graph():
    adj_list = {i: [i+1] for i in range(3)}
    adj_list = AdjList(adj_list, 4)
    vcf_seq = ["AAAG", "GGCC", "CTTT", "AAAG"]
    vcf_seq = [s.lower() for s in vcf_seq]
    return VCFGraph([4, 4, 4, 4], adj_list, SNPs(), vcf_seq)


def test_identical_graph_translation(ins_graph, vcf_ins_graph):
    obg_seq = ["AAT", "TGG", "CCA", "ATT", "GGC"]
    interval = obg.Interval(0, 3, [1, 3, 4, 5], graph=ins_graph)
    indexed_interval = obg.NumpyIndexedInterval.from_interval(interval)
    sequence_graph = obg.SequenceGraph.create_empty_from_ob_graph(ins_graph)
    [sequence_graph.set_sequence(i, seq) for i, seq in enumerate(obg_seq, 1)]
    path = Path([0, 2, 3, 4], [0, 3, 6, 9, 12])
    translation = get_translation(simple_graph, sequence_graph,
                                  indexed_interval, vcf_ins_graph, path)
    assert translation == {i: (i-1, 0) for i in range(1, 6)}


def test_simple_path_translation(simple_graph, simple_vcf_graph):
    obg_seq = ["AAAG", "G", "GCC", "CT", "TT", "AAA", "G"]
    interval = obg.Interval(0, 1, list(range(1, 8)), graph=simple_graph)
    indexed_interval = obg.NumpyIndexedInterval.from_interval(interval)
    sequence_graph = obg.SequenceGraph.create_empty_from_ob_graph(simple_graph)
    [sequence_graph.set_sequence(i, seq) for i, seq in enumerate(obg_seq, 1)]
    path = Path([0, 1, 2, 3], [0, 4, 8, 12, 16])
    translation = get_translation(simple_graph, sequence_graph,
                                  indexed_interval, simple_vcf_graph, path)
    assert translation == {1: (0, 0),
                           2: (1, 0), 3: (1, 1),
                           4: (2, 0), 5: (2, 2),
                           6: (3, 0), 7: (3, 3)}
