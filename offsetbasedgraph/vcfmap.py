import numpy as np
import offsetbasedgraph as obg
from itertools import chain
from collections import namedtuple, defaultdict

VCFEntry = namedtuple("VCFEntry", ["pos", "ref", "alt"])
SNP = namedtuple("SNP", ["nodes"])
DEL = namedtuple("DEL", ["edge"])
INS = namedtuple("INS", ["nodes"])


def prune_entry_end(entry):
    for i, (c1, c2) in enumerate(zip(entry.ref[::-1], entry.alt[::-1])):
        if c1 != c2:
            return VCFEntry(entry.pos, entry.ref[:len(entry.ref)-i], entry.alt[:len(entry.alt)-i])
    i = min(len(entry.alt), len(entry.ref))
    return VCFEntry(entry.pos, entry.ref[:len(entry.ref)-i], entry.alt[:len(entry.alt)-i])


def prune_entry(entry):
    for i, (c1, c2) in enumerate(zip(entry.ref, entry.alt)):
        if c1 != c2:
            return VCFEntry(entry.pos+i, entry.ref[i:], entry.alt[i:])
    i = min(len(entry.ref), len(entry.alt))
    return VCFEntry(entry.pos+i, entry.ref[i:], entry.alt[i:])


def get_vcf_entries(filename):
    def get_entries(line):
        parts = line.split("\t")
        pos = int(parts[1])-1
        ref = parts[3]
        alts = parts[4].split(",")
        return (VCFEntry(pos, ref.lower(), alt.lower()) for alt in alts)

    return chain.from_iterable(
        get_entries(line) for line in open(filename)
        if not line.startswith("#"))


def paralell_nodes_func(graph, linear_path):
    path_nodes = linear_path.nodes_in_interval()

    def get_paralell_nodes(node):
        prev_node = max(
            (linear_path.get_offset_at_node(-prev_node), -prev_node)
            for prev_node in graph.reverse_adj_list[-node]
            if -prev_node in path_nodes)[1]
        paralell_nodes = graph.adj_list[prev_node]
        assert node in paralell_nodes
        return paralell_nodes
    return get_paralell_nodes


def prev_node_func(graph, linear_path):
    path_nodes = linear_path.nodes_in_interval()

    def get_prev_nodes(node):
        return max(
            (linear_path.get_offset_at_node(-prev_node), -prev_node)
            for prev_node in graph.reverse_adj_list[-node]
            if -prev_node in path_nodes)[1]
    return get_prev_nodes


def next_node_func(graph, linear_path):
    path_nodes = linear_path.nodes_in_interval()

    def get_next_node(node):
        return min(
            (linear_path.get_offset_at_node(new_node), new_node)
            for new_node in graph.adj_list[node]
            if new_node in path_nodes)[1]
    return get_next_node


def prune_SNP(entry):
    idxs = [i for i, (c1, c2) in enumerate(zip(entry.ref, entry.alt)) if (c1 != c2)]
    assert len(idxs) == 1, entry
    i = idxs[0]
    return VCFEntry(entry.pos+i, entry.ref[i], entry.alt[i])


def prune_deletion(entry):
    pruned = prune_entry(prune_entry_end(entry))
    assert len(pruned.alt) == 0, pruned
    return pruned


def prune_insertion(entry):
    pruned = prune_entry(prune_entry_end(entry))
    assert len(pruned.ref) == 0, pruned
    return pruned


def prune(entry):
    if len(entry.ref) == len(entry.alt):
        return prune_SNP(entry)
    elif len(entry.ref) > len(entry.alt):
        return prune_deletion(entry)
    else:
        return prune_insertion(entry)


class EdgeMap:
    def __init__(self, adj_list):
        self.indices = adj_list.values()
        self._values = np.zeros_like(adj_list.values)
        self.node_id_offset = adj_list.node_id_offset
        
    def __getitem__(self, edge):
        from_node = edge[0]
        # Returns all edges for a nod
        index = item - self.node_id_offset
        if index < 0:
            return []
        if index >= len(self._indices):
            return []

        start = self._indices[index]
        end = self._indices[index] + self._n_edges[index]
        return self._values[start:end]        

def entry_to_edge_func(graph, reference, seq_graph):
    get_paralell_nodes = paralell_nodes_func(graph, reference)
    get_next_node = next_node_func(graph, reference)
    get_prev_node = prev_node_func(graph, reference)

    def get_SNP_nodes(entry, node):
        next_node = get_next_node(node)
        paralell_nodes = [
            par for par in get_paralell_nodes(node)
            if seq_graph.get_sequence_on_directed_node(par)==entry.alt and next_node in list(graph.adj_list[par])]
        return SNP(paralell_nodes)

    def get_insertion_node(entry, node):
        assert len(entry.alt) <= 64, entry  # Current limit
        paralell_nodes = get_paralell_nodes(node)
        seqs = [seq_graph.get_sequence_on_directed_node(par) for par in paralell_nodes]
        nextss = [graph.adj_list[par] for par in paralell_nodes]
        if len(entry.alt) <= 32:
            valid_nodes = [par for (par, seq, nexts) in zip(paralell_nodes, seqs, nextss)
                           if seq == entry.alt and node in nexts]
        else:
            fulls = [p for p, seq in zip(paralell_nodes, seqs)
                     if seq == entry.alt[:32]]
            assert len(fulls) == 1, (entry, seqs, paralell_nodes)
            paralell_nodes = graph.adj_list[fulls[0]]
            seqs = [seq_graph.get_sequence_on_directed_node(par) for par in paralell_nodes]
            nextss = [graph.adj_list[par] for par in paralell_nodes]
            valid_nodes = [par for (par, seq, nexts) in zip(paralell_nodes, seqs, nextss)
                           if seq == entry.alt[32:] and node in nexts]
            assert len(valid_nodes) == 1, (entry, paralell_nodes, seqs, nextss)
            return (fulls[0], valid_nodes[0])
        assert len(valid_nodes) == 1, (entry, paralell_nodes, seqs, nextss)
        return INS((valid_nodes[0],))

    def get_deletion_edge(entry, node):
        assert reference.get_node_offset_at_offset(entry.pos) == 0
        prev_node = get_prev_node(node)
        deletion_len = len(entry.ref)
        paralell_nodes = get_paralell_nodes(node)
        while deletion_len > 0:
            deletion_len -= graph.node_size(node)
            node = get_next_node(node)
        assert deletion_len == 0
        assert node in paralell_nodes
        return DEL((prev_node, node))

    def entry_to_edge(full_entry):
        entry = prune(full_entry)
        node_offset = int(
            reference.get_node_offset_at_offset(entry.pos))
        assert node_offset == 0 or not (entry.alt and entry.ref), entry
        node = int(reference.get_node_at_offset(entry.pos))
        ref_seq = seq_graph.get_sequence_on_directed_node(node)

        if entry.ref and entry.alt:
            assert ref_seq == entry.ref, entry
            return get_SNP_nodes(entry, node)
        if entry.alt:
            return get_insertion_node(entry, node)
        return get_deletion_edge(entry, node)

    return entry_to_edge


def make_var_map(graph, variants):
    snp_map = np.zeros_like(graph.node_indexes)
    insertion_map = np.zeros_like(graph.node_indexes)
    snp_bucket = {}
    for i, var in variants:
        if type(var) == SNP:
            snp_map[var.nodes[0]-graph.min_node] = i
            if len(var.nodes) > 1:
                snp_bucket[i] = var.nodes

        elif type(var) == INS:
            for node in var.nodes:
                assert insertion_map[node-graph.min_node] == 0
                insertion_map[node-graph.min_node] = i


def get_interval_nodes_touching_variants(interval, linear_map):
    if interval.region_paths[0] < 0:
        assert np.all(interval.region_paths < 0)
        interval = interval.get_reverse()

    linear_path_nodes = linear_map.nodes_in_interval()
    variant_nodes = set(interval.region_paths) - linear_path_nodes

    # Find deletions
    prev_linear_node = None  # Previous node that was on linear
    nodes = interval.region_paths
    for i, node in enumerate(nodes):
        if i > 0:
            if nodes[i-1] not in linear_path_nodes and node in linear_path_nodes:
                # Have passed a deletion, add it
                variant_nodes.add((prev_node_on_linear, node))

        if node in linear_path_nodes:
            prev_node_on_linear = node

    return variant_nodes


def parse_variants(filename):
    parts = (line.split("\t") for line in open(filename))
    return ((int(part[0]), eval(part[1])) for part in parts)


if __name__ == "__main__":
    graph = obg.Graph.from_file("1.nobg")
    snp_files = "variant_map.tsv"
    make_var_map(graph, parse_variants(snp_files))
    exit()

    seq_graph = obg.SequenceGraph.from_file("1.nobg.sequences")
    reference = obg.NumpyIndexedInterval.from_file("1_linear_pathv2.interval")
    vcf_entries = get_vcf_entries("1_variants_cut.vcf")

    outfile = open("variant_map.tsv", "w")
    ref_nodes = reference.nodes_in_interval()
    entry_to_edge = entry_to_edge_func(graph, reference, seq_graph)
    counter = 0
    variants = (entry_to_edge(entry) for entry in vcf_entries)
    for t in enumerate(variants):
        outfile.write("%s\t%s\n" % t)
    outfile.close()
