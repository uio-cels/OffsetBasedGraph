from .graph import Graph, BlockCollection, Block
import numpy as np
from collections import defaultdict


class GraphWithReversals(Graph):

    def __init__(self, blocks, adj_list, create_reverse_adj_list=True,
                 rev_adj_list=None):

        blocks = BlockCollection(blocks)

        super(GraphWithReversals, self).__init__(blocks, adj_list,
                                                 create_reverse_adj_list=create_reverse_adj_list,
                                                 rev_adj_list=rev_adj_list)

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
