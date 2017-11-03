from .graph import Graph, BlockCollection, Block
import numpy as np
from collections import defaultdict


class GraphWithReversals(Graph):

    def __init__(self, blocks, adj_list, create_reverse_adj_list=True,
                 rev_adj_list=None):

        blocks = BlockCollection(blocks)
        m = max(blocks.keys())+1
        self._node_sizes = np.zeros(m, dtype="int32")
        for block_id, block in blocks.items():
            if block_id > 0:
                self._node_sizes[block_id] = block.length()
        super(GraphWithReversals, self).__init__(blocks, adj_list,
                                                 create_reverse_adj_list=create_reverse_adj_list,
                                                 rev_adj_list=rev_adj_list)

    def node_size(self, node_id):
        return self._node_sizes[abs(node_id)]

    def block_in_graph(self, block_id):
        if block_id in self.blocks:
            return True

    @staticmethod
    def _get_reverse_edges(adj_list):
        reverse_edges = defaultdict(list)
        for block, edges in adj_list.items():
            for edge in edges:
                reverse_edges[-edge].append(-block)

        return reverse_edges
