import numpy as np
import offsetbasedgraph as obg
from offsetbasedgraph.fullgraph import FullGraph, FullVCFGraph
from offsetbasedgraph.obg_vcf_translation import TranslationBuilder
from offsetbasedgraph.vcfmap import get_all_vcf_entries
from offsetbasedgraph.vcfgraph import construct_graphs
from pyfaidx import Fasta

gpc_path = "/home/knut/Documents/phd/graph_peak_caller/tests/mhc_test_data/"
data_path = "/home/knut/Documents/phd/two_step_graph_mapper/benchmarking/mhc_graph_data/"
obg_base_name = data_path + "6"
vcf_base_name = "test"

path = data_path + "1000genomes_variants.vcf"
fasta_path = data_path + "linear_ref.fa"

chrom_sizes = {
    1:	249250621,
    2:	243199373,
    3:	198022430,
    4:	191154276,
    5:	180915260,
    6: 4970557,
    # 6:	171115067,
    7:	159138663,
    8:	146364022,
    9:	141213431,
    10:	135534747,
    11:	135006516,
    12:	133851895,
    13:	115169878,
    14:	107349540,
    15:	102531392,
    16:	90354753,
    17:	81195210,
    18:	78077248,
    20:	63025520,
    20:	63025520,
    19:	59128983,
    22:	51304566,
    21:	48129895}


def build_vcf_graphs():
    fasta = Fasta(fasta_path)
    entries = get_all_vcf_entries(path)
    for chrom, graphs in construct_graphs(entries, chrom_sizes, fasta):
        graph, ref, _ = graphs
        graph.save((vcf_base_name + "_graph"))
        ref.save((vcf_base_name + "_ref "))


def build_translation():
    obg_full_graph = FullGraph.from_files(obg_base_name)
    vcf_full_graph = FullVCFGraph.from_files(vcf_base_name)
    t = TranslationBuilder(obg_full_graph, vcf_full_graph)
    translator = t.build()
    translator.save(vcf_base_name)


if __name__ == "__main__":
    # build_vcf_graphs()
    build_translation()
