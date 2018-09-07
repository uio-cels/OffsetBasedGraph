from .vcfmap import *
from .graph import Graph
from .indexedinterval import NumpyIndexedInterval

from collections import Counter
from itertools import chain

AnalysisResults = namedtuple("AnalysisResults", ["A_id", "A_count", "B_id", "B_count", "total"])


def interval_to_variants_func(reference, graph):
    reference_nodes = reference.nodes_in_interval()
    get_next_node = next_node_func(graph, reference)

    def interval_to_variants(interval):
        interval.graph = graph
        if interval.region_paths[0] < 0:
            assert all(rp < 0 for rp in interval.region_paths)
            interval = interval.get_reverse()
        nodes = interval.region_paths
        variant_nodes = [node-graph.min_node for node in nodes if node not in reference_nodes] # Offset nodes to array idx
        variant_edges = [edge for edge in zip(nodes[:-1], nodes[1:])
                         if all(node in reference_nodes for node in edge) and
                         edge[1] != get_next_node(edge[0])]
        return variant_nodes, variant_edges

    return interval_to_variants

def get_variants_from_intervals(reference, graph, intervals):
    interval_to_variants = interval_to_variants_func(reference, graph)
    return (interval_to_variants(interval) for interval in intervals)


def get_haplotypes_func(preceneses):
    def compatible_haplotypes(variants):
        valid_haplotypes = set(range(1137))
        return set.intersection(
            valid_haplotypes,
            *(set(haplotypes) for haplotypes in preceneses[list(variants)]))
    return compatible_haplotypes


def _analyze_variants(haplotypes_list):
    haplotypes_list = list(haplotypes_list)
    N = len(haplotypes_list)
    counter = Counter(chain.from_iterable(haplotypes_list))
    A_id, A_count = counter.most_common(1)[0]
    remaining = [h for h in haplotypes_list if A_id not in h]
    counter2 = Counter(chain.from_iterable(remaining))
    B_id, B_count = counter2.most_common(1)[0]
    return AnalysisResults(A_count, A_id, B_count, B_id, N)


def analyze_interval_set_func(precences, reference, graph, variant_maps):
    variants_to_haplotypes = get_haplotypes_func(precences)
    interval_to_variants = interval_to_variants_func(reference, graph)
    variant_maps = load_variant_maps(chromosome, folder)
    def analyze_intervals(intervals):
        variant_lists = (interval_to_variants(interval) for interval in intervals)
        haplotype_sets = [variants_to_haplotypes(variants) for variants in variant_lists]
        return _analyze_variants(haplotype_sets)

    return analyze_intervals


def get_ids_from_variants(variant_maps, variants_list):
    return [nodes_edges_to_variant_ids(variants[0], variants[1], variant_maps)
            for variants in variants_list]


def pipeline_func_for_chromosome(chromosome, folder="./"):
    variant_maps = load_variant_maps(chromosome, folder)
    precences = load_precences(chromosome, folder)
    reference = NumpyIndexedInterval.from_file(folder+chromosome+"_linear_pathv2.interval")
    graph = Graph.from_file(folder+chromosome+".nobg")
    interval_to_variants = interval_to_variants_func(reference, graph)
    analyze = analyze_interval_set_func(precences, reference, graph)

    def pipeline(intervals):
        variant_lists = (interval_to_variants(interval) for interval in intervals)
        variant_id_sets = get_ids_from_variants(variant_maps, variant_lists)
        return analyze(variant_id_sets)

    return pipeline
