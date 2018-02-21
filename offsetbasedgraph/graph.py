from collections import defaultdict
from .interval import Interval, Position
from .translation import Translation
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
        return super().__getitem__(abs(key))

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
                 rev_adj_list=None):
        """
        Inits the graph with a list of blocks and an adjency list
        :param blocks: dict{block_id: block_length}
        :param adj_list: dict{block_id: [neighbour_ids...]}
        """
        if isinstance(blocks, dict):
            blocks = BlockCollection(blocks)
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

    def connect_postitions(self, position_a, position_b):
        """Connect position_a to position_b with an edge
        and split the region paths if necessary.
        Create translation object and new graph

        :param position_a: Position
        :param position_b: Position
        :returns: New graph and Translation object
        :rtype: (Graph, Translation)

        """
        rp_a = position_a.region_path_id
        rp_b = position_b.region_path_id
        offset_a = position_a.offset
        offset_b = position_b.offset
        l_a = self.blocks[rp_a].length()
        l_b = self.blocks[rp_b].length()
        a_to_b = {}
        b_to_a = {}
        rp_out = rp_a
        rp_in = rp_b
        assert offset_a < l_a
        assert offset_b < l_b, "%s, %s" % (offset_b, l_b)

        if offset_a < l_a - 1:
            # Split first block
            idf = self._next_id()
            idl = self._next_id()
            rp_out = idf
            a_to_b[rp_a] = [Interval(0, l_a-offset_a-1,
                                     [idf, idl])]
            b_to_a[idf] = [Interval(0, offset_a+1, [rp_a])]
            b_to_a[idl] = [Interval(offset_a+1, l_a, [rp_a])]

        if offset_b > 0:
            # Split last block
            idf = self._next_id()
            idl = self._next_id()
            rp_in = idl
            a_to_b[rp_b] = [Interval(0, l_b-offset_b, [idf, idl])]
            b_to_a[idf] = [Interval(0, offset_b, [rp_b])]
            b_to_a[idl] = [Interval(offset_b, l_b, [rp_b])]

        trans = Translation(a_to_b, b_to_a, graph=self)
        n_graph = trans.translate_subgraph(self)
        n_graph.adj_list[rp_out].append(rp_in)
        n_graph.reverse_adj_list[rp_in].append(rp_out)
        trans.graph2 = n_graph
        self._update_a_b_graph(trans._a_to_b, n_graph)
        self._update_a_b_graph(trans._b_to_a, self)
        return n_graph, trans

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

    def merge(self, intervals):
        """Merge the given intervals in the graph, and return
        the resulting graph after merge and a translation object.

        :param intervals: list of intevals to merge
        :returns: translation object, resulting graph
        :rtype: (Graph, Translation)
        """
        block_lengths = self.blocks.copy()
        for b in block_lengths:
            block_lengths[b] = self.blocks[b].length()

        [self.assert_interval_in_graph(i) for i in intervals]
        original_graph = intervals[0].graph
        # Assume same lengths for all intervals
        length = intervals[0].length()
        for interval in intervals[1:]:
            assert interval.length() == length, \
                "All intervals should have the same length (%d, %d)" % \
                (interval.length(), length)

        # 1: Translate graph so that all intervals starts and ends at rps
        trans, graph1 = self._get_inslulate_translation(intervals,
                                                        block_lengths)
        trans.graph2 = graph1
        # correct, a to b interval's graph is wrong for some reason
        self._update_a_b_graph(trans._a_to_b, graph1)

        # Update intervals to be on graph1
        new_intervals = []
        for interval in intervals:
            mp = trans.translate_interval(interval)
            new_interval = mp.get_single_path_intervals()[0]
            new_interval.graph = graph1
            new_intervals.append(new_interval)

        intervals = new_intervals

        # Step 2: Translate each interval to one single block
        # In next graph, each interval is replaced by one single block
        a_b = {}
        b_a = {}
        new_blocks = []
        for interval in intervals:
            # Create one new large block
            large_block = graph1._next_id()
            new_blocks.append(large_block)

            # Back translation from block to interval
            b_a[large_block] = [interval]
            block_lengths[large_block] = interval.length()

            # Find forward translations
            prev_offset = 0
            offset = 0
            for rp in interval.region_paths:
                offset += graph1.blocks[rp].length()
                a_b[rp] = [Interval(prev_offset, offset, [large_block])]
                a_b[rp][0].set_length_cache(offset - prev_offset)
                prev_offset = offset

        # Find new graph
        trans2 = Translation(a_b, b_a, graph1, block_lengths)
        # Update graph in a_b intervals to graph2 (which is now created)

        # Step 3: Merge each block (representing one interval each) to one single block
        new_block = graph1._next_id()
        a_b = {}
        b_a = {new_block: []}

        for block in new_blocks:
            a_b[block] = [Interval(0, length, [new_block])]
            a_b[block][0].set_length_cache(length)
            i = Interval(0, length, [block])
            b_a[new_block].append(i)
            i.set_length_cache(length)

        block_lengths[new_block] = length

        # Translate graph
        trans3 = Translation(a_b, b_a, None, block_lengths)

        # Step 4: Divide large middle block
        starts = []
        for interval in intervals:
            prev_start = 0
            for rp in interval.region_paths:
                prev_start += interval.graph.blocks[rp].length()
                starts.append(prev_start)

        starts = list(set(starts))
        starts.sort()

        # Create translation from big block to all small created from each start
        a_b = {new_block: []}
        b_a = {}
        prev_start = 0
        new_blocks = []
        last_length = 0
        sum_len = 0
        for start in starts:
            new_small = graph1._next_id()
            new_blocks.append(new_small)
            i = Interval(prev_start, start, [new_block])
            i.set_length_cache(start - prev_start)
            b_a[new_small] = [i]
            block_lengths[new_small] = start - prev_start
            last_length = start - prev_start
            sum_len += (start - prev_start)
            prev_start = start

        i2 = Interval(0, last_length, new_blocks)
        i2.set_length_cache(sum_len)
        a_b[new_block].append(i2)

        trans4 = Translation(a_b, b_a, None, block_lengths)

        final_trans = trans + trans2
        final_trans += trans3
        final_trans += trans4
        final_trans.graph1 = original_graph  # Should not be needed ?!

        final_graph = final_trans.translate_subgraph(self)
        final_trans.graph2 = final_graph

        return final_graph, final_trans

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

    @staticmethod
    def block_origin(name):
        """
        Checks if block is merged, alt or main based on the name

        :param name: block name
        :return: merged, alt or main
        """
        name = str(name)
        if name.count("chr") == 1 and "alt" not in name:
            return "main"
        elif name.count("chr") == 2:
            return "merged"
        elif "alt" in name:
            return "alt"
        return "main"

    def find_critical_blocks(self, start_block):
        """Find all critical_blocks starting from
        start block

        :param start_block: block id of start_block
        :returns: List of critical block ids
        :rtype: list(str)

        """
        raise NotImplementedError("Not implemented for graph with reversals")
        cur_block = start_block
        counter = 0
        critical_blocks = []
        visited = []

        while(self.adj_list[cur_block]):
            visited.append((cur_block, counter))
            if counter == 0:
                critical_blocks.append(cur_block)
            nexts = self.adj_list[cur_block]
            counter += (len(nexts) - 1)
            cur_block = nexts[0]
            counter -= (len(self.reverse_adj_list[-cur_block])-1)

        if (counter == 0):
            critical_blocks.append(cur_block)

        return critical_blocks

    def find_parallell_blocks(self, blocks, filter_func):
        """Find all parallell_blocks to block
        that satisfies filter_func

        :param block: region path id
        :param filter_func: predicate function
        :returns: all parallell_blocks
        :rtype: list(str)

        """
        if isinstance(blocks, str):
            blocks = [blocks]
        prev_blocks = [b for b in self.reverse_adj_list[blocks[0]]
                       if filter_func(b)]
        assert len(prev_blocks) == 1, blocks
        start_block = prev_blocks[0]
        next_blocks = [b for b in self.adj_list[blocks[-1]]
                       if filter_func(b)]
        assert len(next_blocks) == 1, blocks
        end_block = next_blocks[0]
        cur_block = start_block
        parallell_blocks = []
        while True:
            next_blocks = [b for b in self.adj_list[cur_block]
                           if filter_func(b)]
            assert len(next_blocks) == 1, str(self.adj_list[cur_block]) + cur_block
            cur_block = next_blocks[0]
            if cur_block == end_block:
                break
            parallell_blocks.append(cur_block)
        return parallell_blocks

    def find_all_critical_blocks(self):
        """Find critical blocks in the graph.
        I.e.  blocks that are traversed by all paths

        :returns: list of critical_blocks ids
        :rtype: list(str)

        """

        start_blocks = [block for block in self.blocks if
                        not self.reverse_adj_list[block]]

        critical_blocks = []
        for start_block in start_blocks:
            critical_blocks.extend(
                self.find_critical_blocks(start_block))

        return critical_blocks

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

    def create_subgraph_from_blocks(self, blocks, alt_locus=None):
        """
        Creates a subgraph using existing edges and only the blocks
        send as argument

        :param blocks: list of block ids
        :param alt_locus: If not None, alt loci blocks not from
            this alt locus will not be added to graph.
            This will typically be parallell alt loci that
            potentially can be added
        :return: Returns a new graph
        """
        blocks = set(blocks)
        # add prev and next critical

        new_edges = defaultdict(list)
        new_blocks = {}
        for b in blocks:
            new_blocks[b] = Block(self.blocks[b].length())

        for b in blocks:
            for e in self.adj_list[b]:
                if alt_locus is not None and\
                   Graph.block_origin(e) == "alt" and alt_locus not in e:
                    continue

                new_edges[b].append(e)
                if e not in new_blocks:
                    new_blocks[e] = Block(self.blocks[e].length())

        # Go through all added blocks, add edges into other added blocks
        for b in new_blocks:
            for edge in self.adj_list[b]:
                if edge in new_blocks and edge not in new_edges[b]:
                    new_edges[b].append(edge)

        subgraph = Graph(new_blocks, new_edges)

        # Add all blocks going into first blocks in graph
        firsts = subgraph.get_first_blocks()
        for f in firsts:
            for before in self.reverse_adj_list[f]:

                if alt_locus is not None \
                   and Graph.block_origin(before) == "alt" \
                   and alt_locus not in before:
                    continue

                new_blocks[before] = Block(self.blocks[before].length())
                new_edges[before].append(f)

        subgraph = Graph(new_blocks, new_edges)

        # If two last rps, they should be going to the same next
        lasts = subgraph.get_last_blocks()
        assert len(lasts) == 1 or len(lasts) == 2, "%s is lasts" % lasts

        if len(lasts) == 2:
            # Connect the two last to next
            if len(self.adj_list[lasts[0]]) > 0 \
               and len(self.adj_list[lasts[1]]) > 0:
                assert self.adj_list[lasts[0]][0] == self.adj_list[lasts[1]][0], \
                    "Not identical next: %s != %s" % (self.adj_list[lasts[0]][0], self.adj_list[lasts[1]][0])

                next = self.adj_list[lasts[0]]
                assert len(next) == 1
                next = next[0]
                last = self.adj_list[lasts[0]][0]
                new_blocks[last] = Block(self.blocks[last].length())
                new_edges[lasts[0]].append(next)
                new_edges[lasts[1]].append(next)

        subgraph2 = Graph(new_blocks, new_edges)
        return subgraph2

    def max_block_id(self):
        return max([id for id in self.blocks.keys()])

    def _next_id(self):
        """Make a new id and return it

        :returns: new id
        :rtype: int

        """
        self._id += 1
        return self._id

    def _get_inslulate_translation(self, intervals, block_lengths={}):
        """Get translation for splitting the region paths
        at the start and end of the intervals such that
        all the intervals span complete region paths

        :param intervals:
        :returns:
        :rtype: Translation

        """
        for i in intervals:
            self.assert_interval_in_graph(i)
        # Split starts
        starts = [i.start_position for i in intervals]
        translation = Translation(graph=self)
        translation.graph2 = self
        cur_graph = self
        i = 0
        id_a = None
        id_b = None

        for start in starts:
            i += 1
            if start.offset == 0:
                continue
            rp = start.region_path_id
            offset = start.offset
            id_a, id_b = (self._next_id(), self._next_id())
            L = self.blocks[rp].length()
            prev_graph = cur_graph.copy()

            trans_dict = {}
            reverse_dict = {}
            trans_dict[rp] = [Interval(
                Position(id_a, 0),
                Position(id_b, L-offset))]

            reverse_dict[id_a] = [Interval(Position(rp, 0),
                                           Position(rp, offset),
                                           graph=prev_graph)]
            block_lengths[id_a] = offset

            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L),
                                           graph=prev_graph)]
            block_lengths[id_b] = L - offset

            tmp_trans = Translation(trans_dict, reverse_dict, graph=cur_graph)
            cur_graph = tmp_trans.translate_subgraph(cur_graph)
            self._update_a_b_graph(trans_dict, cur_graph)

            translation += tmp_trans

        if id_a is not None and id_b is not None:
            assert id_a in cur_graph.blocks
            assert id_b in cur_graph.blocks

        # Translate intervals to new graph
        new_intervals = [translation.translate(interval)
                         for interval in intervals]
        ends = [interval.end_position for interval in new_intervals]

        # Update interval's graph. Is this necessary?
        # Should be done in translate interval
        for interval in new_intervals:
            interval.graph = cur_graph

        back_graph_end = cur_graph
        end_translations = Translation(graph=cur_graph)
        end_translations.graph2 = back_graph_end

        # Split end of intervals in new graph
        for end in ends:
            assert cur_graph is not None
            assert end.region_path_id in cur_graph.blocks
            rp = end.region_path_id
            L = cur_graph.blocks[rp].length()
            offset = end.offset
            if offset == L:
                continue

            id_a, id_b = (self._next_id(), self._next_id())

            trans_dict = {}
            reverse_dict = {}
            trans_dict[rp] = [Interval(
                Position(id_a, 0),
                Position(id_b, L-offset)
                )]

            prev_graph = cur_graph
            reverse_dict[id_a] = [Interval(Position(rp, 0),
                                           Position(rp, offset),
                                           graph=prev_graph)]
            block_lengths[id_a] = offset

            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L),
                                           graph=prev_graph)]
            block_lengths[id_b] = L - offset

            tmp_trans = Translation(trans_dict, reverse_dict, graph=prev_graph)
            cur_graph = tmp_trans.translate_subgraph(cur_graph)

            self._update_a_b_graph(tmp_trans._a_to_b, cur_graph)

            tmp_trans.graph2 = cur_graph

            # Asssert translations now have all intervals needed
            self._assert_translation_intervals_has_graphs(translation)

            end_translations += tmp_trans

        translation += end_translations

        return translation, cur_graph

    def _assert_translation_intervals_has_graphs(self, translation):
        for intervals in translation._a_to_b.values():
            for interval in intervals:
                assert interval.graph is not None

        for intervals in translation._b_to_a.values():
            for interval in intervals:
                assert interval.graph is not None

    def _update_a_b_graph(self, a_b, graph, prev_graph=None):
        for block, intvs in a_b.items():
            for interval in intvs:
                interval.graph = graph
                if prev_graph is not None:
                    interval.set_length_cache(
                        prev_graph.blocks[block].length())

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

    def get_arbitrary_linear_graph(self, return_as_interval=False):
        """
        Returns a linear graph consisiting of an arbitrary path through the original graph,
        and a translation to this graph.
        """

        trans_forward = {}
        trans_back = {}
        new_blocks = {}
        new_adj_list = {}
        visited = {}

        i = 0

        for start_block in self.get_first_blocks():
            chrom_name = "chr_artificial_%d" % i
            linear_chromosome_size = 0
            offset = 0
            linear_blocks = []
            current_block = start_block
            last_block_length = 0
            while True:
                visited[current_block] = True
                linear_blocks.append(current_block)
                #new_blocks[current_block] = self.blocks[current_block]
                block_length = self.blocks[current_block].length()
                trans_forward[current_block] = [Interval(offset, offset + block_length, [chrom_name])]
                offset += block_length
                last_block_length = block_length
                next_blocks = self.adj_list[current_block]

                if len(next_blocks) < 1:
                    break

                #new_adj_list[current_block] = [next_blocks[0]]
                next_block = None
                for block in next_blocks:
                    if block in visited:
                        print("Visited %d" % block)
                        print(visited)
                        continue
                    else:
                        next_block = block
                        break

                if next_block is None:
                    break

                current_block = next_block

            new_blocks[chrom_name] = Block(offset)
            # Forward trans
            trans_back[chrom_name] = [Interval(0, last_block_length, linear_blocks, self)]

            if return_as_interval:
                return Interval(0, last_block_length, linear_blocks, self)

            i += 1

        graph = Graph(new_blocks, new_adj_list)


        for intervals in trans_forward.values():
            for interval in intervals:
                interval.graph = graph

        trans = Translation(trans_forward, trans_back, self)
        return graph,trans

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

    def to_msgpack(self, file_name):
        import msgpack
        from io import BytesIO
        # Write msgpack file
        with open(file_name, 'w') as outfile:
            msgpack.pack(data, outfile)

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
                 rev_adj_list=None):

        if isinstance(blocks, np.ndarray):
            blocks = BlockArray(blocks)
        elif isinstance(blocks, BlockArray):
            pass
        else:
            blocks = BlockCollection(blocks)
        super(Graph, self).__init__(
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
