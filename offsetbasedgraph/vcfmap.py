import numpy as np
import offsetbasedgraph as obg
import pickle
from itertools import chain
from collections import namedtuple

VCFEntry = namedtuple("VCFEntry", ["pos", "ref", "alt"])
SNP = namedtuple("SNP", ["nodes"])
DEL = namedtuple("DEL", ["edge"])
INS = namedtuple("INS", ["nodes"])
VariantMap = namedtuple("VariantMap", ["snps", "insertions", "deletions"])


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
    deletion_map = {}
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
        elif type(var) == DEL:
            deletion_map[var.edge] = i
    return VariantMap(snp_map, insertion_map, deletion_map)


def parse_variants(filename):
    parts = (line.split("\t") for line in open(filename))
    return ((int(part[0]), eval(part[1])) for part in parts)


def write_variants(chromosome, folder):
    base_name = folder + "/" + str(chromosome)
    graph = obg.Graph.from_file(base_name+".nobg")
    seq_graph = obg.SequenceGraph.from_file(base_name + ".nobg.sequences")
    reference = obg.NumpyIndexedInterval.from_file(
        base_name + "_linear_pathv2.interval")
    vcf_entries = get_vcf_entries(base_name + "_variants_cut.vcf")
    outfile = open(base_name + "_variant_map.tsv", "w")
    entry_to_edge = entry_to_edge_func(graph, reference, seq_graph)
    variants = (entry_to_edge(entry) for entry in vcf_entries)
    for t in enumerate(variants):
        outfile.write("%s\t%s\n" % t)
    outfile.close()


def load_variant_maps(variant_maps, base_name):
    snps = np.load(base_name+"_snp_map.npz")
    insertions = np.load(base_name+"_ins_map.npz")
    deletions = pickle.load(open(base_name+"_del_map.pickle", "rb"))
    return VariantMap(snps=snps, insertions=insertions, deletions=deletions)


def write_variant_maps(variant_maps, base_name):
    np.save(base_name+"_snp_map.npz", variant_maps.snps)
    np.save(base_name+"_ins_map.npz", variant_maps.insertions)
    with open(base_name + "_del_map.pickle", "wb") as f:
        pickle.dump(variant_maps.deletions, f)

    np.save(base_name+"ins_map.npz", variant_maps.insertions)


def get_variant_precences(vcf_file_name):
    def precence_from_line(haplo_types, variant_number):
        precence = np.array([variant_number in [int(a), int(b)] for a, b in haplo_types])
        haplotypes = np.flatnonzero(precence)
        return haplotypes

    def get_precences_from_line(line):
        parts = line.split("\t", 9)
        n_variants = len(parts[4].split(","))
        precence = parts[-1].replace(".", "0").replace("/", "|").split("\t")
        haplo_types = [part.split(":", 2)[0].split("|") for part in precence]
        return (precence_from_line(haplo_types, i+1) for i in range(n_variants))

    return chain.from_iterable(
        get_precences_from_line(line) for line in open(vcf_file_name)
        if not line.startswith("#"))


def get_variant_maps(chromosome, folder):
    base_name = folder + "/" + str(chromosome)
    graph = obg.Graph.from_file(base_name+ ".nobg")
    snp_files = base_name + "_variant_map.tsv"
    make_var_map(graph, parse_variants(snp_files))


def get_precences(chromosome, folder="./"):
    base_name = folder + "/" + str(chromosome)
    with open(base_name+"_precences.npy", "rb") as in_file:
        while True:
            try:
                yield np.load(in_file)
            except OSError:
                break


def write_precences(chromosome, folder="./"):
    base_name = folder + "/" + str(chromosome)
    precences = get_variant_precences(base_name + "_variants_small.vcf")
    with open(base_name+"_precences.npy", "wb") as out_file:
        for precence in precences:
            np.save(out_file, precence)


def nodes_edges_to_variant_ids(nodes, edges, variant_maps):
    snps = variant_maps.snps[nodes]
    insertions = variant_maps.insertions[nodes]
    variants = np.where(insertions > 0, insertions, snps)
    assert np.all(variants)
    deletion_ids = [variant_maps.deletions[edge] for edge in edges]
    return set(chain(variants, deletion_ids))


def simplify_vcf(chromosome, folder="./"):
    write_variants(chromosome, folder)
    write_precences(chromosome, folder)
