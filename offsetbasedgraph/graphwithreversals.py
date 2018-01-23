from .graph import Graph, BlockCollection, Block
import numpy as np
import json
from collections import defaultdict
import logging
import pickle


class BlockArray:
    def __init__(self, array):
        if isinstance(array, dict):
            array = self.from_dict(array)
        assert isinstance(array, np.ndarray), type(array)
        self._array = array

        self.node_id_offset = 0  # Subtracted when indexing

    def save(self, file_name):
        np.save(file_name, self._array)

    @classmethod
    def load(cls, filename):
        return cls(np.load(filename))

    @staticmethod
    def from_dict(node_dict):
        max_key = max(node_dict.keys())
        array = np.zeros(max_key+1, dtype="uint8")
        for key, val in node_dict.items():
            array[key] = val.length()
        return array

    def node_size(self, node_id):
        return self._array[abs(node_id) - self.node_id_offset]

    def __contains__(self, node_id):
        node_id = abs(node_id) - self.node_id_offset
        return node_id > 0 and node_id < len(self._array)

    def __iter__(self):
        return self.keys()

    def keys(self):
        return (i + self.node_id_offset for i, v in enumerate(self._array) if v > 0)

    def values(self):
        return (Block(v) for v in self._array if v > 0)

    def items(self):
        return ((i + self.node_id_offset, Block(v)) for i, v in enumerate(self._array) if v > 0)

    def __getitem__(self, node_id):
        v = self._array[abs(node_id) - self.node_id_offset]
        assert v > 0
        return Block(v)




class GraphWithReversals(Graph):

    def __init__(self, blocks, adj_list,
                 create_reverse_adj_list=True,
                 rev_adj_list=None):

        if isinstance(blocks, np.ndarray):
            blocks = BlockArray(blocks)
        elif isinstance(blocks, BlockArray):
            pass
        else:
            blocks = BlockCollection(blocks)
        super(GraphWithReversals, self).__init__(
            blocks, adj_list,
            create_reverse_adj_list=create_reverse_adj_list,
            rev_adj_list=rev_adj_list)

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
        return self.blocks.node_size(node_id)

    def to_numpy_files(self, base_file_name):
        logging.info("Writing blocks to file")
        self.blocks.save(base_file_name + ".npy")
        logging.info("Writing edges to file")

        if False and isinstance(self.adj_list, AdjListAsMatrix):
            self.adj_list.to_file(base_file_name)
        else:
            with open(base_file_name + "edges.pickle", "wb") as f:
                pickle.dump(self.adj_list, f)

        with open(base_file_name + "rev_edges.pickle", "wb") as f:
            pickle.dump(self.reverse_adj_list, f)
        with open(base_file_name + ".node_id_offset", "w") as f:
            f.write(str(self.blocks.node_id_offset))
            print("Wrote node id offset: %d" % self.blocks.node_id_offset)

        """
        with open(base_file_name + "edges.json", "w") as f:
            f.write(json.dumps(self.adj_list))
        logging.info("Writing reverse edges to file")
        with open(base_file_name + "rev_edges.json", "w") as f:
            f.write(json.dumps(self.reverse_adj_list))
        """
    @classmethod
    def from_numpy_files(cls, base_file_name):
        blocks = BlockArray.load(base_file_name + ".npy")

        #adj_list = AdjListAsMatrix.from_file(base_file_name)
        with open(base_file_name + "edges.pickle", "rb") as f:
            adj_list = pickle.loads(f.read())

        #with open(base_file_name + "edges.json") as f:
        #    adj_list = json.loads(f.read())

        with open(base_file_name + "rev_edges.pickle", "rb") as f:
            rev_adj_list = pickle.loads(f.read())
        with  open(base_file_name + ".node_id_offset") as f:
            node_id_offset = int(f.read())

        #with open(base_file_name + "rev_edges.json") as f:
        #    rev_adj_list = json.loads(f.read())
        graph = cls(blocks, adj_list, rev_adj_list=rev_adj_list,
                   create_reverse_adj_list=False)
        graph.blocks.node_id_offset = node_id_offset
        logging.info("Node id offset: %d" % node_id_offset)
        return graph

    @classmethod
    def from_unknown_file_format(cls, base_file_name):
        try:
            graph = cls.from_file(base_file_name)
        except:
            graph = cls.from_numpy_files(base_file_name)

        assert graph is not None, "Graph %s not found" % base_file_name
        return graph

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
        return
        for adjs, other_adjs in [(self.adj_list, self.reverse_adj_list),
                                 (self.reverse_adj_list, self.adj_list)]:
            for node in adjs:
                for edge in adjs[node]:
                    assert -node in other_adjs[-edge]
