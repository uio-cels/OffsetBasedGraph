import numpy as np
import offsetbasedgraph as obg
from offsetbasedgraph.fullgraph import FullGraph, FullVCFGraph
from offsetbasedgraph.obg_vcf_translation import TranslationBuilder, Translator
from pyvg.conversion import vg_json_file_to_interval_collection
from graph_peak_caller.callpeaks import CallPeaks, Configuration
from graph_peak_caller.reporter import Reporter
from graph_peak_caller.intervals import UniqueIntervals
from graph_peak_caller.callpeaks_interface import find_or_create_linear_map
from itertools import chain

gpc_path = "/home/knut/Documents/phd/graph_peak_caller/tests/mhc_test_data/"
data_path = "/home/knut/Documents/phd/two_step_graph_mapper/benchmarking/mhc_graph_data/"
obg_base_name = data_path + "6"
vcf_base_name = "test"


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
    print([blocks[i] for i in range(1, 10)])
    graph = obg.Graph(blocks, edges)
    print([graph.blocks[i] for i in range(1, 10)])
    graph.convert_to_numpy_backend("uint32")
    print(graph.blocks._array[:10])
    for i, b in blocks.items():
        assert graph.blocks._array[i] == b.length(), (i, b, graph.blocks._array[i])
    return graph


obg_full_graph = FullGraph.from_files(obg_base_name)
vcf_full_graph = FullVCFGraph.from_files(vcf_base_name)
# t = TranslationBuilder(obg_full_graph, vcf_full_graph)
# print("BUILDING")
# translator = t.build()
# translator.save("test")
translator = Translator.load("test")
interval_collection = vg_json_file_to_interval_collection(gpc_path + "chips.json", obg_full_graph.graph)
intervals = list(interval_collection)
counter = 0
obg_graph = translate_graph(vcf_full_graph.graph)
new_intervals = [translate_interval(interval, translator, vcf_full_graph.graph, obg_graph)
                 for interval in intervals]

print("--------------------------------------")
for i, i2 in zip(intervals, new_intervals):
    if not i2.length() == i.length():
        print(i, i2)
sample = UniqueIntervals(new_intervals)
control = UniqueIntervals(new_intervals)
# sample = UniqueIntervals(list(intervals))
# control = UniqueIntervals(list(intervals))

config = Configuration()
config.read_length = 70
config.fragment_length = 120
config.linear_map_name = "lin_map.npz"
# obg_graph = obg_full_graph.graph
linear_map = find_or_create_linear_map(obg_graph, config.linear_map_name)

callpeaks = CallPeaks(obg_graph, config, Reporter("testrun"))


callpeaks.run(sample, control)

print(counter)
print("FINISHED")

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
