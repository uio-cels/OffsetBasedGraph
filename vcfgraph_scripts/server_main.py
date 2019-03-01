import numpy as np
import offsetbasedgraph as obg
import graph_peak_caller as gpc
from offsetbasedgraph.fullgraph import FullGraph, FullVCFGraph
from offsetbasedgraph.obg_vcf_translation import TranslationBuilder, Translator
from offsetbasedgraph.vcfgraph import construct_graph
from offsetbasedgraph.vcfmap import get_vcf_entries
from pyvg.conversion import vg_json_file_to_interval_collection
from graph_peak_caller.callpeaks import CallPeaks, Configuration
from graph_peak_caller.reporter import Reporter
from graph_peak_caller.intervals import UniqueIntervals
from graph_peak_caller.callpeaks_interface import find_or_create_linear_map
from itertools import chain
from pyfaidx import Fasta
data_path = "/data/bioinf/human_1pc/"
vcf_path = data_path + "filtered_20.vcf"
fasta_path = data_path + "hg19_chr1-Y.fa"
gpc_path = "/data/bioinf/benchmarking/data/HUMAN_CTCF_ENCSR000DUB/1/filtered_low_qual_reads_removed_20.json"
# gpc_path = "/home/knut/Documents/phd/graph_peak_caller/tests/mhc_test_data/"
obg_base_name = data_path + "20"
vcf_base_name = "20_test"

def build_vcf_graph():
    fasta = Fasta(fasta_path)["20"]
    graph, ref, _ = construct_graph(get_vcf_entries(vcf_path), 63025520, fasta)
    graph.save("20_test_graph")
    ref.save("20_test_ref")

def build_translation():
    obg_full_graph = FullGraph.from_files(obg_base_name)
    print(max(obg_full_graph.graph.blocks._array))
    vcf_full_graph = FullVCFGraph.from_files(vcf_base_name)
    t = TranslationBuilder(obg_full_graph, vcf_full_graph)
    translator = t.build()
    translator.save("20_test")
    

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


def translate_intervals():
    obg_full_graph = FullGraph.from_files(obg_base_name)
    vcf_full_graph = FullVCFGraph.from_files(vcf_base_name)
    translator = Translator.load("20_test")
    interval_collection = vg_json_file_to_interval_collection(gpc_path, obg_full_graph.graph)
    intervals = list(interval_collection)
    counter = 0
    obg_graph = translate_graph(vcf_full_graph.graph)
    obg_graph.to_file("20_small.npz")
    new_intervals = [translate_interval(interval, translator, vcf_full_graph.graph, obg_graph)
                     for interval in intervals]
    print("----------------------------------")
    counter = 0
    for i, i2 in zip(intervals, new_intervals):

        if not i2.length() == i.length():
            counter += 1
            if abs(i2.length()-i.length())>5:
                print(i, i2)
    
    obg.IntervalCollection(new_intervals).to_file("20_test_translated.intervalcollection")

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
    

def run_callpeaks():
    obg_graph = obg.Graph.from_file("20_small.npz")
    new_intervals = obg.IntervalCollection.from_file("20_test_translated.intervalcollection")
    new_intervals2 = obg.IntervalCollection.from_file("20_test_translated.intervalcollection")
    sample = UniqueIntervals(new_intervals)
    control = UniqueIntervals(new_intervals2)
    
    config = Configuration()
    config.read_length = 70
    config.fragment_length = 120
    config.linear_map_name = "lin_map.npz"
    linear_map = find_or_create_linear_map(obg_graph, config.linear_map_name)
    
    callpeaks = CallPeaks(obg_graph, config, Reporter("testrun"))
    
    
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
    # build_vcf_graph()
    # build_translation()
    # translate_intervals()
    # run_callpeaks()
    # run_old_callpeaks()
    compare_peaks()
# obg_full_graph = FullGraph.from_files(obg_base_name)
# vcf_full_graph = FullVCFGraph.from_files(vcf_base_name)
# t = TranslationBuilder(obg_full_graph, vcf_full_graph)
# # print("BUILDING")n
# # translator = t.build()
# # exit()
# 
# # translator.save("test")
# translator = Translator.load("test")
# interval_collection = vg_json_file_to_interval_collection(gpc_path + "chips.json", obg_full_graph.graph)
# intervals = list(interval_collection)
# counter = 0
# obg_graph = translate_graph(vcf_full_graph.graph)
# new_intervals = [translate_interval(interval, translator, vcf_full_graph.graph, obg_graph)
#                  for interval in intervals]
# for i, i2 in zip(intervals, new_intervals):
#     if not i2.length() == i.length():
#         print(i, i2)
# sample = UniqueIntervals(new_intervals)
# control = UniqueIntervals(new_intervals)
# # sample = UniqueIntervals(list(intervals))
# # control = UniqueIntervals(list(intervals))
# 
# config = Configuration()
# config.read_length = 70
# config.fragment_length = 120
# config.linear_map_name = "lin_map.npz"
# # obg_graph = obg_full_graph.graph
# linear_map = find_or_create_linear_map(obg_graph, config.linear_map_name)
# 
# callpeaks = CallPeaks(obg_graph, config, Reporter("testrun"))
# 
# 
# callpeaks.run(sample, control)
# 
# print(counter)
# print("FINISHED")
# 
# 
# 
# 
# 
# 
# 
# graph = obg.Graph.from_file(base_name+".nobg")
# seq_graph = obg.SequenceGraph.from_file(base_name + ".nobg.sequences")
# reference = obg.NumpyIndexedInterval.from_file(
#     base_name + "_linear_pathv2.interval")
# vcf_entries = get_vcf_entries(path)
# outfile = open("6_variant_map.tsv", "w")
# entry_to_edge_tmp = entry_to_edge_func(graph, reference, seq_graph)
# 
# 
# def entry_to_edge(full_entry):
#     if full_entry.alt.startswith("<"):
#         return None
#     try:
#         return entry_to_edge_tmp(full_entry)
#     except AssertionError:
#         return None
# 
# variants = (entry_to_edge(entry) for entry in vcf_entries)
# variants = (v for v in variants if v is not None)
# var_maps = make_var_map(graph, enumerate(variants))
# write_variant_maps(var_maps, "6")
# # 
# # for t in enumerate(variants):
# #     outfile.write("%s\t%s\n" % t)
# # outfile.close()
