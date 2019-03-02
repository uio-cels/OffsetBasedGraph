import numpy as np
import offsetbasedgraph as obg
import graph_peak_caller as gpc
from offsetbasedgraph.fullgraph import FullGraph, FullVCFGraph
from offsetbasedgraph.obg_vcf_translation import TranslationBuilder, Translator
from offsetbasedgraph.vcfgraph import construct_graph, construct_graphs
from offsetbasedgraph.vcfmap import get_vcf_entries, get_all_vcf_entries
from pyvg.conversion import vg_json_file_to_interval_collection
from graph_peak_caller.callpeaks import CallPeaks, Configuration
from graph_peak_caller.reporter import Reporter
from graph_peak_caller.intervals import UniqueIntervals
from graph_peak_caller.callpeaks_interface import find_or_create_linear_map
from pyfaidx import Fasta
import logging
logging.basicConfig(filename="logfile.out", level="INFO")

chrom_sizes = {
    1:	249250621,
    2:	243199373,
    3:	198022430,
    4:	191154276,
    5:	180915260,
    6:	171115067,
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
    19:	59128983,
    22:	51304566,
    21:	48129895}

data_path = "/data/bioinf/human_1pc/"
out_path = data_path+"smallgraph/"
vcf_path = data_path + "filtered_%s.vcf"
all_vcf_path = data_path + "filtered.vcf"
fasta_path = data_path + "hg19_chr1-Y.fa"
gpc_path = "/data/bioinf/benchmarking/data/HUMAN_CTCF_ENCSR000DUB/1/filtered_low_qual_reads_removed_%s.json"
# gpc_path = "/home/knut/Documents/phd/graph_peak_caller/tests/mhc_test_data/"
obg_base_name = data_path + "%s"
vcf_base_name = out_path + "%s_test"


def build_vcf_graphs():
    print("Building VCF Graphs")
    fasta = Fasta(fasta_path)
    entries = get_all_vcf_entries(all_vcf_path)
    for chrom, graphs in construct_graphs(entries, chrom_sizes, fasta):
        print("CHROMOSOME %s " % chrom)
        graph, ref, _ = graphs
        graph.save((vcf_base_name + "_graph") % chrom)
        ref.save((vcf_base_name + "_ref ") % chrom)


def build_vcf_graph(i=20):
    fasta = Fasta(fasta_path)[str(i)]
    graph, ref, _ = construct_graph(get_vcf_entries(vcf_path % i), chrom_sizes[i], fasta)
    graph.save("%s_test_graph" % i)
    ref.save("%s_test_ref" % i)


def build_translation(i=20):
    print("Building traslation: %s" % i)
    obg_full_graph = FullGraph.from_files(obg_base_name % i)
    vcf_full_graph = FullVCFGraph.from_files(vcf_base_name % i)
    t = TranslationBuilder(obg_full_graph, vcf_full_graph)
    translator = t.build()
    translator.save(vcf_base_name % i)
    

def translate_interval(interval, translator, graph, ob_graph):
    reverse = interval.region_paths[0] < 0
    if reverse:
        interval = interval.get_reverse()
    vcf_interval = translator.translate_interval(interval, graph)
    new_interval = interval.__class__(vcf_interval.start, vcf_interval.end,
                                      list(vcf_interval.node_ids+1), graph=ob_graph)
    if reverse:
        new_interval = new_interval.get_reverse()
    return new_interval


def translate_graph(vcf_graph):
    blocks = {i+1: obg.Block(int(v)) for i, v in enumerate(vcf_graph._node_lens)}
    assert np.all(vcf_graph._node_lens > 0)
    edges = {i+1: vcf_graph._adj_list[i]+1 for
             i in range(vcf_graph._node_lens.size)}
    graph = obg.Graph(blocks, edges)
    graph.convert_to_numpy_backend("uint32")
    for i, b in blocks.items():
        assert graph.blocks._array[i] == b.length(), (i, b, graph.blocks._array[i])
    return graph


def translate_intervals(i=20):
    print("Translating Intervals")
    obg_full_graph = FullGraph.from_files(obg_base_name % i)
    vcf_full_graph = FullVCFGraph.from_files(vcf_base_name % i)
    translator = Translator.load(vcf_base_name % i)
    interval_collection = vg_json_file_to_interval_collection(gpc_path % i, obg_full_graph.graph)
    intervals = list(interval_collection)
    counter = 0
    obg_graph = translate_graph(vcf_full_graph.graph)
    obg_graph.to_file(out_path + "%s_small.npz" % i )
    new_intervals = [translate_interval(interval, translator, vcf_full_graph.graph, obg_graph)
                     for interval in intervals]
    print("----------------------------------")
    counter = 0
    for i1, i2 in zip(intervals, new_intervals):

        if not i2.length() == i1.length():
            counter += 1
            if abs(i2.length()-i1.length())>5:
                print(i, i2)
    
    obg.IntervalCollection(new_intervals).to_file("%s_test_translated.intervalcollection" % i)
    print(counter)

def run_old_callpeaks():
    graph = FullGraph.from_files(obg_base_name).graph
    interval_collection = vg_json_file_to_interval_collection(gpc_path, graph)
    interval_collection2 = vg_json_file_to_interval_collection(gpc_path, graph)
    sample = UniqueIntervals(interval_collection)
    control = UniqueIntervals(interval_collection2)
    config = Configuration()
    config.read_length = 70
    config.fragment_length = 120
    config.linear_map_name = "lin_map_old.npz"
    linear_map = find_or_create_linear_map(graph, config.linear_map_name)
    callpeaks = CallPeaks(graph, config, Reporter("testold"))
    callpeaks.run(sample, control)
    

def run_callpeaks(i=20):
    print("Running Callpeaks on %s" % i)
    obg_graph = obg.Graph.from_file(out_path + "%s_small.npz" % i)
    new_intervals = obg.IntervalCollection.from_file("%s_test_translated.intervalcollection" % i)
    new_intervals2 = obg.IntervalCollection.from_file("%s_test_translated.intervalcollection" % i)
    sample = UniqueIntervals(new_intervals)
    control = UniqueIntervals(new_intervals2)
    
    config = Configuration()
    config.read_length = 34
    config.fragment_length = 141
    config.linear_map_name = out_path +"lin_map_%s.npz" % i
    linear_map = find_or_create_linear_map(obg_graph, config.linear_map_name)
    callpeaks = CallPeaks(obg_graph, config, Reporter("testrun%s" %i))
    callpeaks.run(sample, control)

def compare_peaks():
    new_peaks = gpc.peakcollection.PeakCollection.from_file("testrunmax_paths.intervalcollection", text_file=True)
    old_peaks = gpc.peakcollection.PeakCollection.from_file("testoldmax_paths.intervalcollection", text_file=True)
    obg_full_graph = FullGraph.from_files(obg_base_name)
    vcf_full_graph = FullVCFGraph.from_files(vcf_base_name)
    old_linear = obg_full_graph.linear_path
    old_intervals = old_peaks.to_approx_linear_peaks(old_linear, "20")
    linear_interval = obg.Interval(0, vcf_full_graph.path._distance_to_node[-1]-vcf_full_graph.path._distance_to_node[-2], [int(n) for n in vcf_full_graph.path._node_ids+1], graph=obg.Graph.from_file("20_small.npz"))
    new_linear = obg.NumpyIndexedInterval.from_interval(linear_interval)
    new_intervals = new_peaks.to_approx_linear_peaks(new_linear, "20")
    new_starts = sorted([n.start for n in new_intervals.peaks])
    old_starts = sorted([n.start for n in old_intervals.peaks])
    diffs = [abs(n-o) for n, o in zip(new_starts, old_starts)]
    print(sum(diffs), sum(diffs)/len(diffs))



if __name__ == "__main__":
    build_vcf_graphs()
    for i in range(1, 22):
        build_translation(i)
        translate_intervals(i)
        run_callpeaks(i)
        # run_old_callpeaks()
        # compare_peaks()
