from .graph import Graph, BlockCollection, Block
import numpy as np
from collections import defaultdict
import logging


class GraphWithReversals(Graph):

    def __init__(self, blocks, adj_list, create_reverse_adj_list=True,
                 rev_adj_list=None):

        blocks = BlockCollection(blocks)
        #logging.info("Finding max block id")
        if len(blocks) > 0:
            m = max(blocks.keys())+1
        else:
            m = 1

        self._node_sizes = np.zeros(m, dtype="int32")

        #logging.info("Setting node sizes")
        for block_id, block in blocks.items():
            if block_id > 0:
                self._node_sizes[block_id] = block.length()
        #logging.info("Creating reverse adj list")
        super(GraphWithReversals, self).__init__(blocks, adj_list,
                                                 create_reverse_adj_list=create_reverse_adj_list,
                                                 rev_adj_list=rev_adj_list)

        #self.assert_correct_edge_dicts()
        #logging.info("Init graph done.")

    def _possible_node_ids(self):
        node_ids = list(self.blocks.keys())
        possible_ids = node_ids + [-n for n in node_ids]
        return possible_ids

    def get_last_blocks(self):
        """
        :param graph: Graph
        :return: Returns a list of all blocks having no incoming edges
        :rtype: list(Graph)
        """
        return [b for b in self._possible_node_ids() if self._is_end_block(b)]

    def _is_start_block(self, node_id):
        has_edges_out = bool(len(self.adj_list[node_id]))
        has_edges_in = bool(len(self.reverse_adj_list[-node_id]))
        return has_edges_out and (not has_edges_in)

    def _is_end_block(self, node_id):
        has_edges_out = bool(len(self.adj_list[node_id]))
        has_edges_in = bool(len(self.reverse_adj_list[-node_id]))
        return (not has_edges_out) and has_edges_in

    def get_first_blocks(self):
        """
        :param graph: Graph
        :return: Return a list of all blocks having no incoming edges
        :rtype: list(Graph)
        """
        print("#######################################")
        return [b for b in self._possible_node_ids()
                if self._is_start_block(b)]

    def _add_edge(self, block_a, block_b):
        """Add edge from a to b, in both adj_list
        and reveres_adj_list

        :param block_a: from block
        :param block_b: to block
        """
        print("!!!!!", block_a, block_b)
        self.adj_list[block_a].append(block_b)
        self.reverse_adj_list[-block_b].append(-block_a)

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

    def assert_correct_edge_dicts(self):
        logging.info("Asserting edges are correct")
        for adjs, other_adjs in [(self.adj_list, self.reverse_adj_list),
                                 (self.reverse_adj_list, self.adj_list)]:
            for node in adjs:
                for edge in adjs[node]:
                    assert -node in other_adjs[-edge]
