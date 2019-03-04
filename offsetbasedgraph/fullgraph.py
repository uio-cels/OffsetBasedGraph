import offsetbasedgraph as obg
from offsetbasedgraph.vcfgraph import VCFGraph, Path, IndexedPath,\
    Sequences, AdjList
import numpy as np
from itertools import takewhile, chain
from collections import deque, namedtuple
import logging
Variant = namedtuple("Variant",
                     ["ref_node_ids", "alt_node_ids", "ref_seq",  "alt_seq", "ref_node"])


class FullGraph:
    def __init__(self, graph, seq_graph, linear_path):
        self.graph = graph
        self.seq_graph = seq_graph
        self.linear_path = linear_path
        self.indexed_path = IndexedPath.from_indexed_interval(linear_path)
        self._visited = np.zeros(self.graph.blocks.max_node_id()+1)

    @classmethod
    def from_files(cls, base_name):
        graph = obg.Graph.from_file(base_name+".nobg")
        seq_graph = obg.SequenceGraph.from_file(base_name + ".nobg.sequences")
        reference = obg.NumpyIndexedInterval.from_file(
            base_name + "_linear_pathv2.interval")
        return cls(graph, seq_graph, reference)


class FullVCFGraph:
    def __init__(self, graph, path, node_offset=0):
        self.graph = graph
        self.path = path
        self.node_offset = node_offset

    def get_next_node(self, node_id):
        return min(
            (self.path.distance_to_node_id(new_node), new_node)
            for new_node in self.graph._adj_list[node_id]
            if self.path.is_in_path(new_node))[1]

    def get_prev_node(self, node):
        return max(
            (self.path.distance_to_node_id(prev_node), prev_node)
            for prev_node in self.graph._rev_adj_list[node]
            if self.path.is_in_path(prev_node))[1]

    def traverse_ref_nodes(self, node=None):
        return self.path.traverse_ref_nodes(node)

    def ref_nodes_between(self, start_node, end_node):
        return self.path.ref_nodes_between(start_node, end_node)

    def next_node(self, node_id):
        return self.path.next_node(node_id)

    def find_insertion_from_node(self, node, seq):
        node = int(node)
        next_nodes = self.graph._adj_list[node]
        next_ref = self.next_node(node)
        for next_node in next_nodes:
            if self.path.is_in_path(next_node):
                continue
            if next_ref not in self.graph._adj_list[next_node]:
                continue

            next_seq = self.graph._seqs[next_node]
            if next_seq == seq:
                return next_node
        logging.warning("Missing insertion from node %s: %s", node, seq)
        return next_ref

    def interval_length(self, interval):
        node_lens = [self.graph._node_lens[node]
                     for node in interval.node_ids[:-1]]
        if node_lens:
            assert interval.start < node_lens[0], (interval, node_lens[0])
        return sum(node_lens) - interval.start + interval.end

    def get_variants_from_node(self, node_id):
        next_nodes = self.graph._adj_list[node_id]
        stack = deque([[node] for node in self.graph._adj_list[node_id]
                       if not self.path.is_in_path(node) and node_id == self.get_prev_node(node)])
        variants = []
        while stack:
            path = stack.popleft()
            assert isinstance(path, list), (path, list(path))
            node = path[-1]
            next_nodes = self.graph._adj_list[node]
            if any(self.path.is_in_path(n) for n in next_nodes):
                end_node = self.get_next_node(node)
                ref_path = list(self.ref_nodes_between(node_id, end_node)[1:])
                alt_seq = "".join(self.graph._seqs[path])
                ref_seq = "".join(self.graph._seqs[ref_path])
                variants.append(Variant(ref_path, path, ref_seq, alt_seq, node_id))
            else:
                assert isinstance(path, list), (path, list(path))
                stack.extend([path+[n] for n in next_nodes])
        return variants

    @classmethod
    def from_files(cls, base_name):
        graph = VCFGraph.load(base_name + "_graph")
        path = IndexedPath.load(base_name + "_ref")
        return cls(graph, path)

    @classmethod
    def from_full_graph(cls, full_graph):
        graph = VCFGraph.from_obg_graph(full_graph.graph, full_graph.seq_graph)
        path = IndexedPath.from_indexed_interval(full_graph.linear_path)
        return cls(graph, path, full_graph.graph.min_node)
