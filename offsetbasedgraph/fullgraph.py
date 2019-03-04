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

    def traverse_ref_nodes(self, node=None):
        return self.indexed_path.traverse_ref_nodes(node)
        # if node is None:
        #     node = self.linear_path.get_node_at_offset(0)
    def ref_nodes_between(self, start_node, end_node):
        return self.indexed_path.ref_nodes_between(start_node, end_node)
        # return list(takewhile(lambda node: node != end_node,
        #                       self.traverse_ref_nodes(start_node)))

    def get_variants_from_node(self, node_id):
        self._visited[node_id] = 1
        ref_nodes = self.linear_path.nodes_in_interval()
        next_nodes = self.graph.adj_list[node_id]
        deletions = {node for node in next_nodes if node in ref_nodes}
        stack = deque([[node] for node in self.graph.adj_list[node_id]
                       if node not in ref_nodes and node_id == self.get_prev_node(node)])
        variants = []
        while stack:
            path = stack.popleft()
            assert isinstance(path, list), (path, list(path))
            node = path[-1]
            # self._visited[node] = 1
            next_nodes = self.graph.adj_list[node]
            if any(n in ref_nodes for n in next_nodes):
                end_node = self.get_next_node(node)
                ref_path = list(self.ref_nodes_between(node_id, end_node)[1:])
                alt_seq = self.seq_graph.get_nodes_sequences(path)
                ref_seq = self.seq_graph.get_nodes_sequences(ref_path)
                self._visited[path] = 1
                variants.append(Variant(ref_path, path, ref_seq, alt_seq, node_id))
            else:
                assert isinstance(path, list), (path, list(path))
                stack.extend([path+[n] for n in next_nodes])
        return variants

    def get_all_variants(self):
        ref_nodes = self.traverse_ref_nodes()
        return chain.from_iterable(self.get_variants_from_node(node) for node in ref_nodes if len(self.graph.adj_list[node]))

    def get_next_node(self, node_id):
        path_nodes = self.linear_path.nodes_in_interval()
        return min(
            (self.linear_path.get_offset_at_node(new_node), new_node)
            for new_node in self.graph.adj_list[node_id]
            if new_node in path_nodes)[1]

    def get_prev_node(self, node):
        path_nodes = self.linear_path.nodes_in_interval()
        return max(
            (self.linear_path.get_offset_at_node(-prev_node), -prev_node)
            for prev_node in self.graph.reverse_adj_list[-node]
            if -prev_node in path_nodes)[1]

    @classmethod
    def from_files(cls, base_name):
        graph = obg.Graph.from_file(base_name+".nobg")
        seq_graph = obg.SequenceGraph.from_file(base_name + ".nobg.sequences")
        reference = obg.NumpyIndexedInterval.from_file(
            base_name + "_linear_pathv2.interval")
        return cls(graph, seq_graph, reference)


class FullVCFGraph:
    def __init__(self, graph, path):
        self.graph = graph
        self.path = path

    def get_next_node(self, node_id):
        path_nodes = self.path._node_ids
        return min(
            (self.path.distance_to_ndoe_id(new_node), new_node)
            for new_node in self.graph._adj_list[node_id]
            if new_node in path_nodes)[1]

    def get_prev_node(self, node):
        path_nodes = self.path._node_ids
        path_nodes = self.linear_path.nodes_in_interval()
        return max(
            (self.linear_path.get_offset_at_node(-prev_node), -prev_node)
            for prev_node in self.graph.reverse_adj_list[-node]
            if -prev_node in path_nodes)[1]

    def traverse_ref_nodes(self, node=None):
        return self.indexed_path.traverse_ref_nodes(node)

    def ref_nodes_between(self, start_node, end_node):
        return self.indexed_path.ref_nodes_between(start_node, end_node)

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
        self._visited[node_id] = 1
        next_nodes = self.graph.adj_list[node_id]
        stack = deque([[node] for node in self.graph.adj_list[node_id]
                       if not self.path.is_in_path(node) and node_id == self.get_prev_node(node)])
        variants = []
        while stack:
            path = stack.popleft()
            assert isinstance(path, list), (path, list(path))
            node = path[-1]
            next_nodes = self.graph.adj_list[node]
            if any(self.path.is_in_path(n) for n in next_nodes):
                end_node = self.get_next_node(node)
                ref_path = list(self.ref_nodes_between(node_id, end_node)[1:])
                alt_seq = "".join(self.graph._seqs[path])
                ref_seq = "".join(self.graph._seqs[path])
                self._visited[path] = 1
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
        graph = VCFGraph.from_obg_graph(full_graph.graph)
        path = IndexedPath.from_indexed_interval(full_graph.linear_path)
        return cls(graph, path)
