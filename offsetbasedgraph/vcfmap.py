import numpy as np
import offsetbasedgraph as obg
import logging
import pickle
from itertools import chain
from collections import namedtuple, defaultdict

VCFEntry = namedtuple("VCFEntry", ["pos", "ref", "alt"])
SNP = namedtuple("SNP", ["nodes", "trail"])
TRAIL = namedtuple("TRAIL", ["nodes"])
DEL = namedtuple("DEL", ["from_nodes", "to_nodes"])
INS = namedtuple("INS", ["nodes"])
VariantMap = namedtuple("VariantMap",
                        ["snps", "insertions", "deletions", "trails"])

DELMap = namedtuple("DELMap", ["from_ids", "to_ids"])


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
        parts = line.split("\t", 5)
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
    if i != 0:
        logging.debug("%s, %s", i, entry)
    if len(entry.ref) > 1:
        logging.debug(entry)
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


def basename_to_entry_func(chromosome, folder):
    base_name = folder + "/" + str(chromosome)
    graph = obg.Graph.from_file(base_name+".nobg")
    seq_graph = obg.SequenceGraph.from_file(base_name + ".nobg.sequences")
    reference = obg.NumpyIndexedInterval.from_file(
        base_name + "_linear_pathv2.interval")
    return entry_to_edge_func(graph, reference, seq_graph)


def entry_to_edge_func(graph, reference, seq_graph):
    get_paralell_nodes = paralell_nodes_func(graph, reference)
    get_next_node = next_node_func(graph, reference)
    get_prev_node = prev_node_func(graph, reference)

    def traverse_n(node, n):
        while n > 0:
            n -= graph.node_size(node)
            node = get_next_node(node)
        if n == 0:
            return node
        return None

    def get_SNP_nodes(entry, node, full_entry):
        next_node = get_next_node(node)
        paralell_nodes = [
            par for par in get_paralell_nodes(node)
            if seq_graph.get_node_sequence(par)==entry.alt and next_node in graph.adj_list[par]]
        trails = find_trails(node, paralell_nodes, entry, full_entry)
        return SNP(paralell_nodes, trails)

    def get_insertion_node(entry, node):
        assert len(entry.alt) <= 64, entry  # Current limit
        paralell_nodes = get_paralell_nodes(node)
        seqs = [seq_graph.get_node_sequence(par) for par in paralell_nodes]
        nextss = [graph.adj_list[par] for par in paralell_nodes]
        if len(entry.alt) <= 32:
            valid_nodes = [par for (par, seq, nexts) in zip(paralell_nodes, seqs, nextss)
                           if seq == entry.alt and node in nexts]
        else:
            fulls = [p for p, seq in zip(paralell_nodes, seqs)
                     if seq == entry.alt[:32]]
            assert len(fulls) == 1, (entry, seqs, paralell_nodes)
            paralell_nodes = graph.adj_list[fulls[0]]
            seqs = [seq_graph.get_node_sequence(par) for par in paralell_nodes]
            nextss = [graph.adj_list[par] for par in paralell_nodes]
            valid_nodes = [par for (par, seq, nexts) in zip(paralell_nodes, seqs, nextss)
                           if seq == entry.alt[32:] and node in nexts]
            assert len(valid_nodes) == 1, (entry, paralell_nodes, seqs, nextss)
            return (fulls[0], valid_nodes[0])
        assert len(valid_nodes) == 1, (entry, paralell_nodes, seqs, nextss)
        return INS((valid_nodes[0],))

    def get_deletion_edge(entry, node):
        assert reference.get_node_offset_at_offset(entry.pos) == 0
        prev_nodes = [-prev for prev in graph.reverse_adj_list[-node]]
        deletion_len = len(entry.ref)
        paralell_nodes = get_paralell_nodes(node)
        while deletion_len > 0:
            deletion_len -= graph.node_size(node)
            node = get_next_node(node)
        after_node = node
        last_node = get_prev_node(after_node)
        after_nodes = graph.adj_list[last_node]
        assert deletion_len == 0
        if not all(after_node in paralell_nodes for after_node in after_nodes):
            logging.warning("%s, %s: %s, %s", entry, node, list(after_nodes), list(paralell_nodes))
        return DEL(list(prev_nodes), list(after_nodes))

    def find_trails(ref_node, snp_nodes, entry, full_entry):
        if entry == full_entry:
            return []
        after_trail = traverse_n(ref_node, len(full_entry.ref))
        if after_trail is None:
            return []
        trail_seq = full_entry.ref[len(entry.ref):]
        trail_nodes = [node for node in chain.from_iterable(graph.adj_list[snp_node] for snp_node in snp_nodes)
                       if seq_graph.get_node_sequence(node) == trail_seq and
                       node not in reference.nodes_in_interval() and after_trail in graph.adj_list[node]]
        assert all(trail_node in graph.adj_list[ref_node] for trail_node in trail_nodes)
        if trail_nodes:
            logging.debug("--->%s, %s", snp_nodes, trail_nodes)
        return trail_nodes

    def entry_to_edge(full_entry):
        entry = prune(full_entry)
        node_offset = int(
            reference.get_node_offset_at_offset(entry.pos))
        assert node_offset == 0 or not (entry.alt and entry.ref), entry
        node = int(reference.get_node_at_offset(entry.pos))
        ref_seq = seq_graph.get_node_sequence(node)

        if entry.ref and entry.alt:
            assert ref_seq == entry.ref, (entry, ref_seq)
            return get_SNP_nodes(entry, node, full_entry)
        if entry.alt:
            return get_insertion_node(entry, node)
        return get_deletion_edge(entry, node)

    return entry_to_edge


def make_var_map(graph, variants):
    snp_map = np.zeros_like(graph.node_indexes)
    insertion_map = np.zeros_like(graph.node_indexes)
    trail_map = np.zeros_like(graph.node_indexes)
    deletion_map = DELMap(defaultdict(set), defaultdict(set))
    for i, var in variants:
        if type(var) == SNP:
            nodes = [node-graph.min_node for node in var.nodes]
            snp_map[nodes] = i
            trail_nodes = [node-graph.min_node for node in var.trail]
            trail_map[trail_nodes] = i

        elif type(var) == INS:
            for node in var.nodes:
                if not insertion_map[node-graph.min_node] == 0:
                    logging.warning(
                        "Double insertion at node %s: %s+%s",
                        node-graph.min_node, insertion_map[node-graph.min_node], i)
                insertion_map[node-graph.min_node] = i
        elif type(var) == DEL:
            for node in var.from_nodes:
                deletion_map.from_ids[node].add(i)
            for node in var.to_nodes:
                deletion_map.to_ids[node].add(i)

    return VariantMap(snp_map, insertion_map, deletion_map, trail_map)


def parse_variants(filename):
    parts = (line.split("\t") for line in open(filename))
    return ((int(part[0]), eval(part[1])) for part in parts)


def get_variants(vcf_entries, entry_to_edge):
    return (entry_to_edge(entry) for entry in vcf_entries)


def write_variants(chromosome, folder):
    base_name = folder + "/" + str(chromosome)
    graph = obg.Graph.from_file(base_name+".nobg")
    seq_graph = obg.SequenceGraph.from_file(base_name + ".nobg.sequences")
    reference = obg.NumpyIndexedInterval.from_file(
        base_name + "_linear_pathv2.interval")
    vcf_entries = get_vcf_entries(base_name + "_variants.vcf")
    outfile = open(base_name + "_variant_map.tsv", "w")
    entry_to_edge = entry_to_edge_func(graph, reference, seq_graph)
    variants = (entry_to_edge(entry) for entry in vcf_entries)
    for t in enumerate(variants):
        outfile.write("%s\t%s\n" % t)
    outfile.close()


def load_variant_maps(chromosome, folder="./"):
    base_name = folder+chromosome
    snps = np.load(base_name+"_snp_map.npy")
    insertions = np.load(base_name+"_ins_map.npy")
    trails = np.load(base_name+"_trail_map.npy")
    deletions = pickle.load(open(base_name+"_del_map.pickle", "rb"))
    return VariantMap(snps=snps, insertions=insertions, deletions=deletions, trails=trails)


def write_variant_maps(variant_maps, base_name):
    np.save(base_name+"_snp_map.npy", variant_maps.snps)
    np.save(base_name+"_ins_map.npy", variant_maps.insertions)
    np.save(base_name+"_trail_map.npy", variant_maps.trails)
    with open(base_name + "_del_map.pickle", "wb") as f:
        pickle.dump(variant_maps.deletions, f)


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
    var_maps = make_var_map(graph, parse_variants(snp_files))
    write_variant_maps(var_maps, base_name)


def load_precences(chromosome, folder="./"):
    base_name = folder + "/" + str(chromosome)
    return np.load(base_name+"_precences.npy")


def write_precences(chromosome, folder="./"):
    base_name = folder + "/" + str(chromosome)
    precences = np.array(list(get_variant_precences(base_name + "_variants.vcf")))
    np.save(base_name+"_precences.npy", precences)


def variants_to_variant_ids_func(variant_maps, debug_func=None):
    node_variants = np.where(variant_maps.insertions > 0, variant_maps.insertions, variant_maps.snps)

    def nodes_edges_to_variant_ids(nodes, edges):
        variants = node_variants[nodes]
        trails = variant_maps.trails[nodes]
        if not np.all((variants > 0) | (trails > 0)):
            [debug_func(node) for i, node in enumerate(nodes) if variants[i] == 0 and trails[i] == 0]
        deletion_ids = []
        for from_node, to_node in edges:
            pos_from = variant_maps.deletions.from_ids[from_node]
            pos_to = variant_maps.deletions.to_ids[to_node]
            del_ids = pos_from & pos_to
            if not del_ids:
                logging.warning("-DEL: %s (%s->%s)",
                                (from_node, to_node), pos_from, pos_to)
            else:
                deletion_ids.extend(del_ids)
        return set(chain(variants, deletion_ids))

    return nodes_edges_to_variant_ids


def simplify_vcf(chromosome, folder="./"):
    write_variants(chromosome, folder)
    write_precences(chromosome, folder)
    get_variant_maps(chromosome, folder)
