from .graph import Graph, BlockCollection, Block, BlockArray, AdjListAsNumpyArrays
import numpy as np
import json
from collections import defaultdict
import logging
import pickle


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

    def get_sorted_node_ids(self, reverse=False):
        return sorted(self.blocks.keys(), reverse=reverse)

    def __eq__(self, other):
        for node, block in self.blocks.items():
            if node not in other.blocks:
                return False
            if other.blocks[node] != block:
                return False
        for node, block in other.blocks.items():
            if node not in self.blocks:
                return False
            if self.blocks[node] != block:
                return False

        return True

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
        if isinstance(self.blocks, BlockArray):
            logging.info("Using fast way to get first block")
            min_block_id = self.blocks.node_id_offset + 1
            assert min_block_id in self.blocks
            if len(self.reverse_adj_list[min_block_id]) == 0:
                logging.info("Fast way used")
                return [min_block_id]

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

    def uses_numpy_backend(self):
        if isinstance(self.blocks, BlockArray):
            assert isinstance(self.adj_list, AdjListAsNumpyArrays), \
                "If blocks is numpy, edges should also be"
            return True
        return False

    def convert_to_numpy_backend(self):
        if self.uses_numpy_backend():
            logging.warning("Trying to convert to numpy backend, but is already on numpy backend")
            return

        logging.info("Converting to numpy backend...")
        self.blocks = BlockArray.from_dict(self.blocks)
        self.adj_list = AdjListAsNumpyArrays.create_from_edge_dict(self.adj_list)
        self.reverse_adj_list = AdjListAsNumpyArrays.create_from_edge_dict(self.reverse_adj_list)
        logging.info("Conversion finished")

    def convert_to_dict_backend(self):
        new_blocks = {}
        new_adj_list = defaultdict(list)
        new_reverse_adj_list = defaultdict(list)

        i = 0
        for node_id, block in self.blocks.items():
            if i % 100000 == 0:
                logging.info("%d nodes converted" % i)
            i += 1

            new_blocks[node_id] = block
            edges = self.adj_list[node_id]
            if len(edges) > 0:
                new_adj_list[node_id].extend(list(edges))

            edges = self.reverse_adj_list[-node_id]
            if len(edges) > 0:
                new_reverse_adj_list[-node_id].extend(list(edges))

            edges = self.reverse_adj_list[node_id]
            if len(edges) > 0:
                print(" Adding %s to %d" % (list(edges), node_id))
                new_reverse_adj_list[node_id].extend(list(edges))

        self.blocks = new_blocks
        self.adj_list = new_adj_list
        self.reverse_adj_list = new_reverse_adj_list

    def to_numpy_file(self, file_name):
        assert isinstance(self.blocks, BlockArray), "Blocks must be represented as BlockArray"
        assert isinstance(self.adj_list, AdjListAsNumpyArrays), "Edges must be on numpy format"
        assert isinstance(self.reverse_adj_list, AdjListAsNumpyArrays), "Reverse edges must be on numpy format"

        logging.info("Saving to numpy format")
        file = open(file_name, "wb")
        np.savez_compressed(file,
                 blocks=self.blocks._array,
                 node_id_offset = self.blocks.node_id_offset,
                 adj_list_indices=self.adj_list._indices,
                 adj_list_values=self.adj_list._values,
                 adj_list_n_edges=self.adj_list._n_edges,
                 reverse_adj_list_indices=self.reverse_adj_list._indices,
                 reverse_adj_list_values=self.reverse_adj_list._values,
                 reverse_adj_list_n_edges=self.reverse_adj_list._n_edges,
                 reverse_adj_list_node_id_offset=self.reverse_adj_list.node_id_offset
                 )
        file.close()
        logging.info("Graph saved to %s" % file_name)

    @classmethod
    def from_numpy_file(cls, file_name):
        logging.info("Reading from numpy file %s" % file_name)

        try:
            file = open(file_name, "rb")
        except FileNotFoundError:
            try:
                file = open(file_name + ".obg", "rb")
            except FileNotFoundError:
                file = open(file_name + ".nobg", "rb")
        

        data = np.load(file)

        node_id_offset = data["node_id_offset"]

        adj_list = AdjListAsNumpyArrays(
            data["adj_list_indices"],
            data["adj_list_values"],
            data["adj_list_n_edges"],
            node_id_offset+1  # Always one more for edges than for blockarray
        )
        rev_adj_list = AdjListAsNumpyArrays(
            data["reverse_adj_list_indices"],
            data["reverse_adj_list_values"],
            data["reverse_adj_list_n_edges"],
            data["reverse_adj_list_node_id_offset"]  # Always one more for edges than for blockarray
        )

        blocks = BlockArray(data["blocks"])
        blocks.node_id_offset = node_id_offset

        graph = cls(blocks,
                    adj_list=adj_list,
                    rev_adj_list=rev_adj_list)
        file.close()
        logging.info("Done reading from numpy file")
        return graph


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
        logging.info("Reading nodes")
        blocks = BlockArray.load(base_file_name + ".npy")
        logging.info("Reading edges")
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
        logging.info("Initing graph")
        graph = cls(blocks, adj_list, rev_adj_list=rev_adj_list,
                   create_reverse_adj_list=False)
        graph.blocks.node_id_offset = node_id_offset
        return graph

    @classmethod
    def from_unknown_file_format(cls, base_file_name):
        try:
            try:
                graph = cls.from_numpy_file(base_file_name + ".nobg")
                return graph
            except:
                graph = cls.from_numpy_files(base_file_name)
                return graph
        except:
            print("Found no numpy graph. Trying pickle.")

        graph = cls.from_file(base_file_name  + ".obg")
        if graph is None:
            graph = cls.from_file(base_file_name)
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
