from .vcfmap import *
from collections import defaultdict
import logging


class SNPs:
    def __init__(self, node_index=[0], snps=[]):
        self._node_index = np.asanyarray(node_index)
        self._snps = np.asanyarray(snps)

    def __eq__(self, other):
        return True
        if not np.all(self._node_index == other._node_index):
            return False
        return np.all(self._snps == other._snps)


class AdjList:
    def __init__(self, node_index=[0], to_nodes=[]):
        self._node_index = np.asanyarray(node_index)
        self._to_nodes = np.asanyarray(to_nodes)

    def __eq__(self, other):
        print("CHekcking adjlist")
        if not np.all(self._node_index == other._node_index):
            logging.debug("Different node index")
            return False
        if not np.all(self._to_nodes == other._to_nodes):
            logging.debug("Different to nodes")
            return False
        logging.debug("Same adj list")
        return True

    @classmethod
    def from_dict(cls, adj_list, n_nodes):
        node_index = np.empty(n_nodes+1, dtype="int")
        lens = [len(adj_list[from_node]) for from_node in range(n_nodes)]
        node_index[1:] = np.cumsum(lens)
        node_index[0] = 0
        to_nodes = [to_node for from_node in range(n_nodes) for to_node in sorted(adj_list[from_node])]
        return cls(node_index, to_nodes)

    def __repr__(self):
        return "AdjList(%s, %s)" % (self._node_index, self._to_nodes)


class VCFGraph:
    def __init__(self, node_lens, adj_list, snps):
        self._node_lens = np.asanyarray(node_lens)
        self._adj_list = adj_list
        self._snps = snps

    def __eq__(self, other):
        if not np.all(self._node_lens == other._node_lens):
            logging.debug("Different lens")
            return False
        if not self._adj_list == other._adj_list:
            print(self._adj_list, other._adj_list)
            logging.debug("Different adj list")
            return False
        return np.all(self._snps == other._snps)

    def __repr__(self):
        return "VCFGraph(%s, %s)" % (self._node_lens, self._adj_list)

    @classmethod
    def from_vcf(filename):
        with open(filename) as f:
            pass


def classify_vcf_entry(vcf_entry):
    if len(vcf_entry.alt) > len(vcf_entry.ref):
        return INS
    elif len(vcf_entry.alt) < len(vcf_entry.ref):
        return DEL
    return SNP


def graph_from_indels(deletions, insertions, reference_length):
    break_points = np.concatenate((deletions[0]+1, deletions[0]+deletions[1]+1, insertions[0]+1))
    args = np.argsort(break_points, kind="mergesort")
    sorted_break_points = break_points[args]
    print(sorted_break_points)
    diffs = np.diff(sorted_break_points)
    node_starts = np.r_[0, sorted_break_points[diffs>0], sorted_break_points[-1]]
    all_node_starts = np.concatenate((insertions[0]+1, node_starts))
    print(all_node_starts)
    code_args = np.argsort(all_node_starts, kind="mergesort")
    reference_node_ids = code_args[-node_starts.size:]
    insertion_node_ids = code_args[:-node_starts.size]
    adj_list = defaultdict(list)
    for from_node, to_node in zip(reference_node_ids[:-1], reference_node_ids[1:]):
        adj_list[from_node].append(to_node)
    node_start_idxs = np.searchsorted(node_starts, insertions[0]+1)
    from_nodes = reference_node_ids[node_start_idxs-1]
    to_nodes = reference_node_ids[node_start_idxs]
    for from_node, to_node, node_id in zip(from_nodes, to_nodes, insertion_node_ids):
        adj_list[from_node].append(node_id)
        adj_list[node_id].append(to_node)
    reference_node_lens = np.diff(np.r_[node_starts, reference_length])
    print(reference_node_lens)
    insertion_node_lens = insertions[1]
    node_lens = np.empty(reference_node_lens.size+insertion_node_lens.size)
    node_lens[code_args[-node_starts.size:]] = reference_node_lens
    node_lens[code_args[:-node_starts.size]] = insertion_node_lens
    print(node_lens)
    print(adj_list)
    return VCFGraph(node_lens, AdjList.from_dict(adj_list, node_lens.size), SNPs)

def construct_graph(vcf_entries, reference_length):
    cur_node = 0
    insertion_positions = []
    insertion_lens = []
    deletion_starts = []
    deletion_ends = []
    for entry in vcf_entries:
        if entry.alt.startswith("<"):
            continue

        if len(entry.ref) == len(entry.alt):
            continue
        var_type = classify_vcf_entry(entry)
        if var_type == INS:
            try:
                new_entry = prune_insertion(entry)
            except AssertionError:
                print("Invalid insertion")
                continue
            insertion_positions.append(new_entry.pos)
            insertion_lens.append(len(new_entry.alt))
        elif var_type == DEL:
            try:
                new_entry = prune_deletion(entry)
            except AssertionError:
                print("Invalid insertion")
                continue

            deletion_starts.append(new_entry.pos)
            deletion_lens.append(len(new_entry.alt))
        dels = np.array([deletion_starts, deletion_lens])
        insertions = np.array([insertion_positions, insertion_lens])
        graph = graph_from_indels(dels, insertions, reference_length)

    print(len(deletion_starts))
    print(len(insertion_positions))




        
