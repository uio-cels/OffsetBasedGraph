from .vcfmap import *
from .graph import Graph
from .indexedinterval import NumpyIndexedInterval

from collections import Counter
from itertools import chain

AnalysisResults = namedtuple(
    "AnalysisResults", ["A_id", "A_count", "B_id", "B_count", "total"])


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
        variants = list(variants)
        if not variants:
            return valid_haplotypes
        return set.intersection(
            valid_haplotypes,
            *(set(haplotypes) for haplotypes in preceneses[list(variants)]))
    return compatible_haplotypes


def _analyze_variants(haplotypes_list):
    haplotypes_list = list(haplotypes_list)
    N = len(haplotypes_list)
    if N == 0:
        return AnalysisResults(-1, 0, -1, 0, 0)
    counter = Counter(chain.from_iterable(haplotypes_list))
    assert counter or min(len(s) for s in haplotypes_list) == 0, haplotypes_list
    A_id, A_count = counter.most_common(1)[0] if counter else (-1, 0)
    if A_count == N:
        return AnalysisResults(A_id, A_count, -1, 0, N)
                       
    remaining = [h for h in haplotypes_list if A_id not in h]
    assert A_count+len(remaining) == N, (A_count, len(remaining), N, remaining)
    assert remaining, (haplotypes_list, A_id, A_count)
    counter2 = Counter(chain.from_iterable(remaining))
    assert counter2 or min(len(s) for s in remaining) == 0, remaining

    B_id, B_count = counter2.most_common(1)[0] if counter2 else (-1, 0)
    assert A_count+B_count <= N, (A_count, B_count, N)
    return AnalysisResults(A_id, A_count, B_id, B_count, N)


def summarize_results(results, nbins=20):
    table = np.array(results)
    haplo_counts = table[:, 1]
    diplo_counts = table[:, 3] + haplo_counts
    type_ratios = np.hstack((haplo_counts[:, None], diplo_counts[:, None]))/table[:, 4][:, None]
    table_sum = np.sum(table, axis=0)
    type_sums = np.sum(type_ratios, axis=0)
    haplo_hist = np.histogram(type_ratios[:, 0], bins=100, range=(0, 1))
    diplo_hist = np.histogram(type_ratios[:, 1], bins=100, range=(0, 1))
    return np.concatenate(([table.shape[0]], table_sum[[1, 3, 4]], type_sums, haplo_hist[0], diplo_hist[0]))

def analyze_interval_set_func(precences, reference, graph, variant_maps, debug_func=None):
    interval_to_variants = interval_to_variants_func(reference, graph)
    variant_to_variant_ids = variants_to_variant_ids_func(variant_maps, debug_func)
    variant_ids_to_haplotypes = get_haplotypes_func(precences)

    def analyze_intervals(intervals):
        variant_lists = (interval_to_variants(interval) for interval in intervals)
        variant_id_sets = [variant_to_variant_ids(*variant) for variant in variant_lists]
        haplotype_sets = [variant_ids_to_haplotypes(variant_ids) for variant_ids in variant_id_sets]
        return _analyze_variants(haplotype_sets)

    return analyze_intervals


def get_ids_from_variants(variant_maps, variants_list):
    return [nodes_edges_to_variant_ids(variants[0], variants[1], variant_maps)
            for variants in variants_list]


def pipeline_func_for_chromosome(chromosome, folder="./"):
    variant_maps = load_variant_maps(chromosome, folder)
    precences = load_precences(chromosome, folder)
    reference = NumpyIndexedInterval.from_file(folder+chromosome+"_linear_pathv2.interval")
    get_seq = obg.SequenceGraph.from_file(folder+chromosome+".nobg.sequences").get_node_sequence
    graph = Graph.from_file(folder+chromosome+".nobg")


    get_prev_node = prev_node_func(graph, reference)
    get_next_node = next_node_func(graph, reference)

    def get_prev_ref(node):
        nodes = [-node]
        path = []
        while len(nodes) == 1:
            path.append(-nodes[0])
            if -nodes[0] in reference.nodes_in_interval():
                return path[::-1]

            nodes = graph.reverse_adj_list[nodes[0]]

        ref_nodes = [(reference.get_offset_at_node(-node), -node) for node in nodes if -node in reference.nodes_in_interval()]
        if not ref_nodes:
            return None
        return [max(ref_nodes)[1]]+path[::-1]

    # return (max(ref_nodes)[1])

    def get_next_ref(node):
        nodes = [node]
        path = []
        while len(nodes) == 1:
            path.append(nodes[0])
            if nodes[0] in reference.nodes_in_interval():
                return path
            nodes = graph.adj_list[nodes[0]]

        ref_nodes = [(reference.get_offset_at_node(node), node) for node in nodes if node in reference.nodes_in_interval()]
        if not ref_nodes:
            return None
        return path+[min(ref_nodes)[1]]


    def debug_func(node):
        node += graph.min_node
        prev_path = get_prev_ref(node)
        next_path = get_next_ref(node)
        prev_ref = prev_path[0]
        next_ref = next_path[-1]
        if prev_ref is None or next_ref is None:
            logging.warning("(%s, %s)", node, get_seq(node))
            return
        paralell_refs = []
        ref_node = prev_ref
        while ref_node != next_ref:
            ref_node = get_next_node(ref_node)
            paralell_refs.append(ref_node)
        alt_path = prev_path[1:-1]+next_path
        ref_seq = "".join(get_seq(n) for n in paralell_refs)
        alt_seq = "".join(get_seq(n) for n in alt_path)
        if not ref_seq == alt_seq:
            logging.warning("%s: (%s, %s), (%s, %s)", node, alt_path, paralell_refs, prev_path, next_path)

    analyze = analyze_interval_set_func(precences, reference, graph, variant_maps, debug_func)

    def pipeline(intervals):
        return analyze(intervals)

    return pipeline
