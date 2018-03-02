from collections import defaultdict
from .interval import Interval, Position
import pickle
import os
import numpy as np
import logging
from scipy.sparse import dok_matrix, csr_matrix, lil_matrix
import scipy.io


class Block(object):
    def __init__(self, length):
        assert length > 0
        self._length = length

    def length(self):
        return self._length

    def __eq__(self, other):
        return self.length() == other.length()

    def __str__(self):
        return "Block(%d)" % self.length()

    def __repr__(self):
        return self.__str__()


class BlockArray:
    def __init__(self, array):
        if isinstance(array, dict):
            array = self.from_dict(array)
        assert isinstance(array, np.ndarray), type(array)
        self._array = array

        self.node_id_offset = 0  # Subtracted when indexing

    def __len__(self):
        return len(self._array)

    def save(self, file_name):
        np.save(file_name, self._array)

    def max_node_id(self):
        return len(self._array) - 1 + self.node_id_offset

    @classmethod
    def load(cls, filename):
        return cls(np.load(filename))

    @staticmethod
    def from_dict(node_dict):
        max_key = max(node_dict.keys())
        min_key = min(node_dict.keys())
        node_id_offset = min_key

        array = np.zeros((max_key-min_key)+1, dtype="uint8")
        for key, val in node_dict.items():
            array[key-node_id_offset] = val.length()

        block_array = BlockArray(array)
        block_array.node_id_offset = node_id_offset

        return block_array

    def node_size(self, node_id):
        return self._array[abs(node_id) - self.node_id_offset]

    def __contains__(self, node_id):
        node_id = abs(node_id) - self.node_id_offset
        return node_id >= 0 and node_id < len(self._array)

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

    def __str__(self):
        return str(list(self.items()))

    def __repr__(self):
        return self.__str__()


class BlockCollection(dict):
    def __init__(self, *args, **kwargs):
        super(BlockCollection, self).__init__(*args, **kwargs)

    def __contains__(self, name):
        assert isinstance(name, (int, np.integer)) , \
            "Block collection can only deal with numeric block IDs"

        return super().__contains__(abs(name))

    def __getitem__(self, key):
        return super(BlockCollection, self).__getitem__(abs(key))

    def node_size(self, node_id):
        return self[node_id].length()


class AdjListAsNumpyArrays:
    def __init__(self, indices, values, n_edges, node_id_offset=0):
        self._indices = indices  # Mapping from node id to position in values where edge list starts
        self._values = values  # to-edges
        self._n_edges = n_edges  # Number of edges for each node
        self.node_id_offset = node_id_offset

    def __setitem__(self, key, value):
        raise NotImplementedError()

    def __getitem__(self, item):
        # Returns all edges for a nod
        index = item - self.node_id_offset
        if index < 0:
            return []
        if index >= len(self._indices):
            return []

        start = self._indices[index]
        end = self._indices[index] + self._n_edges[index]
        return self._values[start:end]

    @classmethod
    def create_from_edge_dict(cls, edge_dict):
        nodes = edge_dict.keys()
        n_edges = sum([len(edges) for edges in edge_dict.values()])

        i = 0
        sorted_nodes = sorted(nodes)
        min_node_id = sorted_nodes[0]
        max_node_id = sorted_nodes[-1]
        node_span = max_node_id - min_node_id

        indices = np.zeros(node_span+1, dtype=np.int32)
        values = np.zeros(n_edges, dtype=np.int32)
        lengths = np.zeros(node_span+1, dtype=np.int32)

        for node in sorted_nodes:
            index = node - min_node_id
            n_edges_out = len(edge_dict[node])
            lengths[index] = n_edges_out
            indices[index] = i
            #print("Range %d:%d" % (i, i+n_edges_out))
            #print(len(values))
            #print("%s" % edge_dict[node])
            values[i:i+n_edges_out] = edge_dict[node]

            i += n_edges_out

        return cls(indices, values, lengths, min_node_id)


class BaseGraph(object):
    """
    Class for holding an Offset-Based Graph
    and performing simple operations on this graph.

    Does this by storing a dict of Blocks and
    a dict of adjency lists (one list for each block
    that has adjencies).

    >>> blocks = {1: 100000, 2: 50000, 3: 40000}
    >>> adj_list = {1: [2, 3]}
    >>> graph = Graph(blocks, adj_list)
    """

    # Graph alterations
    def __init__(self, blocks, adj_list, create_reverse_adj_list=True,
                 rev_adj_list=None, do_not_create_node_indices=False):
        """
        Inits the graph with a list of blocks and an adjency list
        :param blocks: dict{block_id: block_length}
        :param adj_list: dict{block_id: [neighbour_ids...]}
        """
        if isinstance(blocks, np.ndarray):
            blocks = BlockArray(blocks)
        elif isinstance(blocks, dict):
            blocks = BlockCollection(blocks)
        elif not isinstance(blocks, BlockArray):
            raise ValueError("Blocks must be either a dict, BlockArray or a "
                             "numpy array. Type is %s" % type(blocks))

        self.blocks = blocks

        if not isinstance(adj_list, AdjListAsNumpyArrays):
            self.adj_list = defaultdict(list, adj_list)
        else:
            self.adj_list = adj_list

        if rev_adj_list is not None:
            self.reverse_adj_list = rev_adj_list
        elif create_reverse_adj_list:
            self.reverse_adj_list = self._get_reverse_edges(adj_list)

        if isinstance(self.blocks, BlockArray):
            self._id = self.blocks.max_node_id()
        else:
            self._id = max([b for b in blocks if isinstance(b, int)] + [-1])

        if not do_not_create_node_indices:
            self.node_indexes = self._get_node_indexes()

    def _get_node_indexes(self):
        if isinstance(self.blocks, BlockArray):
            # Quicker way to make node_indexes array
            logging.info("(using cumsum on np block array)")
            node_indexes = np.cumsum(self.blocks._array, dtype=np.uint32)
            logging.info("Node indexes created...")
            self.min_node = (self.blocks.node_id_offset+1)
            return node_indexes

        sorted_nodes = sorted(self.blocks.keys())
        if len(sorted_nodes) == 0:
            return None
        else:
            min_node = sorted_nodes[0]

        self.min_node = min_node
        max_node = sorted_nodes[-1]
        span = max_node-min_node+1
        node_indexes = np.zeros(span+1, dtype=np.uint32)
        offset = 0
        for i, node in enumerate(sorted_nodes):
            index = node - min_node
            node_indexes[index] = offset
            offset += self.node_size(node)
            node_indexes[-1] = offset
        return node_indexes

    def summary(self):
        """Return summary text

        :rtype: str

        """

        n_blocks = len(self.blocks)
        n_edges = sum(len(v) for v in self.adj_list.values())
        return "Graph with %s blocks and %s edges" % (n_blocks, n_edges)

    def copy(self):
        """Make a copy of the graph

        :returns: copy of the graph
        :rtype: Graph

        """
        logging.info("Copying graph")
        new_blocks = {}
        new_adjs = defaultdict(list)
        for b in self.blocks:
            new_blocks[b] = Block(self.blocks[b].length())
        for b in self.adj_list:
            new_adjs[b] = list(self.adj_list[b])

        new_graph = self.__class__(new_blocks, new_adjs, True)
        logging.info("Done copying graph")
        return new_graph

    def to_graph_with_reversals(self):
        from .graphwithreversals import GraphWithReversals
        return GraphWithReversals(self.blocks, self.adj_list)

    @classmethod
    def from_file(cls, file_name):
        """
        Load graph from pickle

        :param file_name: File name
        :rtype: Graph
        """
        if os.path.isfile("%s" % file_name):
            with open("%s" % file_name, "rb") as f:
                obj = pickle.loads(f.read())
                assert isinstance(obj, cls)
                return obj
        else:
            print("Warning: Graph file %s not found" % file_name)
            return None

    def to_file(self, file_name):
        """
        Writes the graph to file so that it later can be
        recreated using the from_file method

        :param file_name: File name
        :return:
        """
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    def remove(self, block_id):
        """Remove a block including edges from the graph

        :param block_id: block id of the block

        """
        del self.blocks[block_id]
        for edge in self.adj_list[block_id]:
            self.reverse_adj_list[edge].remove(block_id)
        del self.adj_list[block_id]

        for edge in self.reverse_adj_list[block_id]:
            self.adj_list[edge].remove(block_id)
        del self.reverse_adj_list[block_id]

    def prev_position(self, pos):
        """
        :param pos: Position
        :return: The previous position in graph
        :rtype: Position
        """
        if pos.offset > 0:
            return Position(pos.region_path_id, pos.offset - 1)
        else:
            # Find previous region path
            prev = None
            for b in self.blocks:
                if pos.region_path_id in self.adj_list[b]:
                    prev = b
                    break
            if prev is None:
                raise Exception("Found no previous pos for position %s" % pos)

            return Position(b, self.blocks[b].length() - 1)

    def next_position(self, pos):
        """
        :param pos: Position
        :return: The the next position in graph
        :rtype: Position
        """
        if pos.offset < self.blocks[pos.region_path_id].length() - 1:
            return Position(pos.region_path_id, pos.offset + 1)
        else:
            # Find next region path
            next = self.adj_list[pos.region_path_id][0]
            return Position(next, 0)

    def assert_position_in_graph(self, position, exclusive=False):
        """Check that a position is in the graph

        :param position: Position
        :param exclusive: Wheter to use exclusive end
        """
        assert position.region_path_id in self.blocks
        block = self.blocks[position.region_path_id]
        assert position.offset < block.length() + int(exclusive)

    def assert_interval_in_graph(self, interval):
        """Check that whole interval is contained in self

        :param interval: Interval
        """
        for rp in interval.region_paths:
            assert rp in self.blocks, "%s not in %s" % (rp, self.blocks.keys())
        self.assert_position_in_graph(interval.start_position)
        self.assert_position_in_graph(interval.end_position, exclusive=True)

    def __str__(self):
        if isinstance(self.adj_list, AdjListAsNumpyArrays):
            edges = str({node: list(self.adj_list[node]) for node in self.blocks
                         if len(self.adj_list[node]) > 0})
        else:
            edges = str(self.adj_list)

        return "Graph: \n Blocks: %s\n Edges: %s" % \
            (self.blocks, edges)

    __repr__ = __str__

    @staticmethod
    def is_main_name(name):
        """Check if name comes from a block on the main path

        :param name: block_id
        :rtype: bool
        """

        if "alt" not in name:
            return True
        if name.count("chr") > 1:
            return True
        return False

    def n_edges_in(self, block):
        """
        Finds and returns the number of edges going in to a block

        :param block:
        :return: Returns the number of edges
        """
        n = 0
        for b in self.blocks:
            if block in self.adj_list[b]:
                n += 1
        return n

    def has_identical_structure(self, other):
        """
        Checks if this graph has identical
        structure (edges and blocks) to other graph.
        Size of region paths is ignores (and can be different).

        :param other: Graph to compare with
        :return: True if identical, otherwise False
        """

        if len(self.blocks) != len(other.blocks):
            print("Different number of blocks")
            return False

        # For every block, check that there exists
        # a block in other graph with the same number of
        # edges in and out
        other_blocks = list(other.blocks.keys()).copy()
        for b in self.blocks:
            match = False
            for ob in other_blocks:
                sim_out = len(self.adj_list[b]) == len(other.adj_list[ob])
                sim_in = self.n_edges_in(b) == other.n_edges_in(ob)
                if sim_out and sim_in:
                    # Remove from list to check, and check next (break)
                    other_blocks.remove(ob)
                    match = True
                    break
            if not match:
                # No match for block b, return False
                return False

        return True

    def max_block_id(self):
        return max([id for id in self.blocks.keys()])

    def _next_id(self):
        """Make a new id and return it

        :returns: new id
        :rtype: int

        """
        self._id += 1
        return self._id

    @staticmethod
    def level_dict(blocks):
        """
        Return dict with block as key and level as value

        :param blocks: block_id
        :rtype: dict

        """
        level_mapping = {"alt": 0,
                         "main": 2,
                         "merged": 1}
        return {b: level_mapping[BaseGraph.block_origin(b)]
                for b in blocks}

    def get_indexed_interval_through_graph(self):
        interval = self.get_arbitrary_interval_through_graph()
        return interval.to_indexed_interval(True)

    def get_arbitrary_interval_through_graph(self):
        logging.info("Getting first blocks")
        start = self.get_first_blocks()
        logging.info("First blocks found")
        assert len(start) == 1, "Only works when graph has one start node"
        nodes = []
        current_block = start[0]
        i = 0
        while True:
            if i % 500000 == 0:
                logging.info("Processing node %d" % i)
            i += 1
            nodes.append(current_block)
            next_blocks = self.adj_list[current_block]

            if len(next_blocks) < 1:
                break

            next_block = next_blocks[0]
            current_block = next_block

        return Interval(0, self.node_size(nodes[-1]), nodes, self)

    def contains_interval(self, interval):
        """
        :param interval:
        :return: Returns true if the interval exists in the graph
        """
        print("Warning not tested. Not correct.")
        for rp in interval.region_paths:
            if rp not in self.blocks:
                return False

        if interval.start_position.offset >= \
                self.blocks[interval.start_position.region_path_id]:
            return False

        if interval.end_position.offset > \
                    self.blocks[interval.end_position.region_path_id]:
            return False

        return True

    def number_of_basepairs(self):
        if isinstance(self.blocks, BlockArray):
            return int(np.sum(self.blocks._array))
        else:
            return sum([b.length() for b in self.blocks.values()])

    def max_node_size(self):
        if isinstance(self.blocks, BlockArray):
            return np.max(self.blocks._array)
        else:
            return max([b.length() for b in self.blocks.values()])



class Graph(BaseGraph):

    def __init__(self, blocks, adj_list,
                 create_reverse_adj_list=True,
                 rev_adj_list=None, do_not_create_node_indices=False):

        super(Graph, self).__init__(
            blocks, adj_list,
            create_reverse_adj_list=create_reverse_adj_list,
            rev_adj_list=rev_adj_list,
            do_not_create_node_indices=do_not_create_node_indices)

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
        if False:
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
