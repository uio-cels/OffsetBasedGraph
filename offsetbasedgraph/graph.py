from collections import defaultdict
from .interval import Interval, Position
from .util import takes
from .translation import Translation
import sys


class Block(object):
    def __init__(self, length):
        self._length = length

    def length(self):
        return self._length

    def to_interval(self):
        pass

    def __eq__(self, other):
        return self.length() == other.length()

    def __str__(self):

        return "B(%s)" % self._length

    def __repr__(self):
        return self.__str__()


class Graph(object):
    adj_list = defaultdict(list)
    reverse_adj_list = defaultdict(list)

    # Graph alterations
    def __init__(self, blocks, adj_list):
        """
        Inits the graph with a list of blocks and an adjency list
        :param blocks:
        :param adj_list:
        :return:
        """
        self.blocks = blocks
        self.adj_list = defaultdict(list, adj_list)
        self.reverse_adj_list = self._get_reverse_edges(adj_list)
        self._id = max(blocks)

    def copy(self):
        """Make a copy of the graph

        :returns: copy of the graph
        :rtype: Graph

        """
        return Graph(self.blocks, self.adj_list)

    def _next_id(self):
        """Make a new id and return it

        :returns: new id
        :rtype: int

        """
        self._id += 1
        return self._id

        # Debug sanity checking
        for block in adj_list:
            assert block in self.blocks, "Edge found from block that is not in blocks"
            for block2 in adj_list[block]:
                assert block2 in self.blocks, "Edge found going to non-existing block"

    @staticmethod
    def from_file(file_name):
        """
        Reads from file and creates the graph
        :param file_name: File name
        :return:
        """
        pass

    def get_edges_as_list(self):
        """
        :return: Returns edges as tuples (from, to) in a list
        """
        edges = []
        for block in self.adj_list:
            for to in self.adj_list[block]:
                edges.append((block, to))

        return edges

    @staticmethod
    def _generate_translation(interval_a, interval_b,
                              sub_intervals_a, sub_intervals_b, ids):

        reverse_dict = {ids: list(pair)
                        for pair in zip(sub_intervals_a, sub_intervals_b)}

        involved_region_paths = interval_a.region_paths+interval_b.region_paths
        front_dict = {rp: [i for i in sub_intervals_a+sub_intervals_b
                           if rp in i.region_paths]
                      for rp in set(involved_region_paths)}

        return Translation(front_dict, reverse_dict)

    def _get_edges_to_subintervals(self, intervals):
        """Return a list of incoming and outgoing edges at the
        divisions between the sub_intervals

        :param intervals: list of consecutive intervals
        :returns: incoming_edges, outgoing_edges
        :rtype: (list, list)

        """
        pre_edges = []
        post_edges = []
        for i_pre, i_post in zip(intervals[:-1], intervals[1:]):
            if not i_post.starts_at_rp():
                post_edges.append([])
                pre_edges.append([])
                continue
            pre_rps = self.reverse_adj_list[i_post.region_paths[0]]
            pre_edges.append(
                [rp for rp in pre_rps if rp != i_pre.region_paths[-1]]
                )

            post_rps = self.adj_list[i_pre.region_paths[-1]]
            post_edges.append(
                [rp for rp in post_rps if rp != i_post.region_paths[0]])

        return pre_edges, post_edges

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

    def _add_edge(self, block_a, block_b):
        """Add edge from a to b, in both adj_list
        and reveres_adj_list

        :param block_a: from block
        :param block_b: to block
        """

        self.adj_list[block_a].append(block_b)
        self.reverse_adj_list[block_b].append(block_a)

    def _join_blocks(self, block_ids):
        """Merge the given blocks in the graph and
        giving them the combination of all edges.

        :param block_ids: block ids of blocks to merge
        :returns: translation object
        :rtype: Translation

        """
        if not block_ids:
            return Translation()
        lens = [self.blocks[_id].length() for _id in block_ids]
        L = lens[0]
        assert all(l == L for l in lens)

        # Find all incoming and outgoing edges
        in_edges = set()
        out_edges = set()
        for block_id in block_ids:
            in_edges.update(self.reverse_adj_list[block_id])
            out_edges.update(self.adj_list[block_id])

        # Add the new block
        new_id = self._next_id()
        self.blocks[new_id] = Block(L)
        for block_id in out_edges:
            self._add_edge(new_id, block_id)
        for block_id in in_edges:
            self._add_edge(block_id, new_id)

        # Remove old blocks
        for block_id in block_ids:
            self.remove(block_id)

        # Create translation object
        translation = Translation(
            {block_id: [Interval(Position(new_id, 0), Position(new_id, L))]
             for block_id in block_ids},
            {new_id:
             [Interval(Position(block_id, 0), Position(block_id, L))
              for block_id in block_ids]})

        return translation

    def _split_block(self, block_id, offsets):
        """Splits the given block at offsets
        and return translation object

        :param block_id: block to split
        :param offset: where to split
        :returns: translation object
        :rtype: Translation

        """
        l = self.blocks[block_id].length()
        if not offsets:
            return Translation()

        # Add blocks
        blocks = [Block(b-a) for a, b in zip([0]+offsets, offsets+[l])]
        ids = [self._next_id() for _ in blocks]
        self.blocks.update(dict(zip(ids, blocks)))

        # Add edges
        for id1, id2 in zip(ids[:-1], ids[1:]):
            self._add_edge(id1, id2)

        for e in self.adj_list[block_id]:
            self._add_edge(ids[-1], e)
        for e in self.reverse_adj_list[block_id]:
            self._add_edge(e, ids[0])

        self.remove(block_id)

        # Set translation
        return Translation(
            {block_id:
             [Interval(
                 Position(ids[0], 0),
                 Position(ids[-1], blocks[-1].length()),
                 ids)]},
            {_id: [Interval(Position(block_id, offset),
                            Position(block_id, offset+block.length()))]
             for _id, block, offset, in zip(ids, blocks, [0]+offsets)}
        )

    def _merge_clean_intervals(self, intervals):
        """Merge intervals that all start at beginning
        of a region path

        :param intervals: 
        :returns: 
        :rtype: Translation

        """
        first_rps = [interval.region_paths[0] for interval in intervals]
        first_rp_lengths = [self.blocks[rp].lenght()
                            for rp in first_rps]
        min_len = min(first_rp_lengths)
        translations = [self._split_block(rp, min_len) for
                        rp in first_rps]

    def get_split_translation(self, rp, offset):

        id_a, id_b = (self._next_id(), self._next_id())
        L = self.blocks[rp].length()
        trans_dict = {}
        reverse_dict = {}
        trans_dict[rp] = [Interval(
            Position(id_a, 0),
            Position(id_b, L-offset))]

        reverse_dict[id_a] = [Interval(Position(rp, 0),
                                       Position(rp, offset),
                                       graph=self)]

        reverse_dict[id_b] = [Interval(Position(rp, offset),
                                       Position(rp, L),
                                       graph=self)]

        trans = Translation(trans_dict, reverse_dict, graph=self)

        new_graph = trans.translate_subgraph(self)
        trans.graph2 = new_graph
        for i_list in trans_dict.values():
            for i in i_list:
                i.graph = new_graph

        return trans, new_graph

    def connect_postitions(self, position_a, position_b):
        """FIXME! briefly describe function

        :param position_a: 
        :param position_b: 
        :returns: 
        :rtype: 

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

        if offset_a < l_a:
            idf = self._next_id()
            idl = self._next_id()
            rp_in = idl
            a_to_b[rp_a] = [Interval(0, l_a-offset_a-1,
                                     [idf, idl])]
            b_to_a[idf] = [Interval(0, offset_a+1, [rp_a])]
            b_to_a[idl] = [Interval(offset_a+1, l_a, [rp_a])]

        if offset_b > 0:
            idf = self._next_id()
            rp_in = idf
            idl = self._next_id()
            a_to_b = self[rp_b] = [Interval(0, l_b-offset_b,
                                            [idf, idl])]
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
        :return: Returns the previous position in graph
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

            return Position(b, self.graph.blocks[b].length() - 1)

    def next_position(self, pos):
        """
        :param pos: Position
        :return: Returns the previous position in graph
        """
        if pos.offset < self.blocks[pos.region_path_id].length():
            return Position(pos.region_path_id, pos.offset + 1)
        else:
            # Find next region path
            next = self.adj_list[pos.region_path_id]
            return Position(next, 0)

    def _get_insulated_merge_transformation(self, intervals):
        """Merge intervals that all start and end at region path
        boundries

        :param intervals:
        :returns: translation and new graph
        :rtype: (Translation, Graph)

        """
        first_rps = [i.region_paths[0] for i in intervals]
        offset = min(self.blocks[rp].length() for rp in first_rps)
        if offset == intervals[0].length():
            return Translation(graph=self), self

        translation = Translation(graph=self)
        graph = self
        for rp in first_rps:
            if self.blocks[rp].length() == offset:
                continue
            trans, graph = graph.get_split_translation(rp, offset)
            translation += trans

        intervals = [translation.translate(i) for i in intervals]
        first_rps = [i.region_paths[0] for i in intervals]

        merge_id = graph._next_id()
        a_to_b = {rp: [Interval(Position(merge_id, 0),
                                Position(merge_id, offset),
                                )]
                  for rp in first_rps}
        b_to_a = {merge_id: [Interval(Position(rp, 0),
                                      Position(rp, offset),
                                      graph=graph)]}

        merge_trans = Translation(a_to_b, b_to_a, graph=graph)
        intervals = [merge_trans.translate(i) for i in intervals]
        full_trans = translation + merge_trans
        new_graph = full_trans.translate_subgraph(graph)
        trans, ng = new_graph._get_insulated_merge_translation(intervals)
        return full_trans+trans, ng

    def _get_inslulate_translation(self, intervals):
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

            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L),
                                           graph=prev_graph)]

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

        # Update interval's graph. Is this necessary? Should be done in translate interval
        for interval in new_intervals:
            interval.graph = cur_graph

        back_graph_end = cur_graph.copy()
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

            prev_graph = cur_graph.copy()
            reverse_dict[id_a] = [Interval(Position(rp, 0),
                                           Position(rp, offset),
                                           graph=prev_graph)]

            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L),
                                           graph=prev_graph)]

            tmp_trans = Translation(trans_dict, reverse_dict, graph=prev_graph)
            cur_graph = tmp_trans.translate_subgraph(cur_graph)

            self._update_a_b_graph(tmp_trans._a_to_b, cur_graph.copy())

            tmp_trans.graph2 = cur_graph.copy()

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

    def _update_a_b_graph(self, a_b, graph):
        for intvs in a_b.values():
            for interval in intvs:
                interval.graph = graph

    def assert_position_in_graph(self, position, exclusive=False):
        assert position.region_path_id in self.blocks
        block = self.blocks[position.region_path_id]
        assert position.offset < block.length() + int(exclusive)

    def assert_interval_in_graph(self, interval):
        for rp in interval.region_paths:
            assert rp in self.blocks, "%s not in %s" % (rp, self.blocks.keys())
        self.assert_position_in_graph(interval.start_position)
        self.assert_position_in_graph(interval.end_position, exclusive=True)

    def merge(self, intervals):
        """Merges the given intervals in the graph, and returns
        a the resulting graph after merge and a translation object.
        :param intervals: list of intevals to merge
        :returns: translation object, resulting graph
        :rtype: Graph, Translation
        """

        [self.assert_interval_in_graph(i) for i in intervals]
        original_graph = intervals[0].graph
        # Assume same lengths for all intervals
        length = intervals[0].length()
        for interval in intervals[1:]:
            assert interval.length() == length, \
                "All intervals should have the same length (%d, %d)" % (interval.length(), length)

        # Step 1: Translate graph so that all intervals starts and ends at region paths
        trans, graph1 = self._get_inslulate_translation(intervals)
        trans.graph2 = graph1
        self._update_a_b_graph(trans._a_to_b, graph1)  # correct, a to b interval's graph is wrong for some reason

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

            # Find forward translations
            prev_offset = 0
            offset = 0
            for rp in interval.region_paths:
                offset += graph1.blocks[rp].length()
                a_b[rp] = [Interval(prev_offset, offset, [large_block])]
                prev_offset = offset

        # Find new graph
        trans2 = Translation(a_b, b_a, graph1)
        graph2 = trans2.translate_subgraph(graph1)
        trans2.graph2 = graph2

        # Update graph in a_b intervals to graph2 (which is now created)
        self._update_a_b_graph(a_b, graph2)

        # Step 3: Merge each block (representing one interval each) to one single block
        new_block = graph2._next_id()
        a_b = {}
        b_a = {new_block: []}

        for block in new_blocks:
            a_b[block] = [Interval(0, length, [new_block])]
            b_a[new_block].append(Interval(0, length, [block], graph2))

        # Translate graph
        trans3 = Translation(a_b, b_a, graph2)
        graph3 = trans3.translate_subgraph(graph2)
        self._update_a_b_graph(a_b, graph3)
        trans3.graph2 = graph3

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
        for start in starts:
            new_small = graph3._next_id()
            new_blocks.append(new_small)
            b_a[new_small] = [Interval(prev_start, start, [new_block], graph3)]
            last_length = start - prev_start
            prev_start = start

        a_b[new_block].append(Interval(0, last_length, new_blocks))

        trans4 = Translation(a_b, b_a, graph3)
        graph4 = trans4.translate_subgraph(graph3)
        self._update_a_b_graph(a_b, graph4)
        trans4.graph2 = graph4
        final_trans = trans3 + trans4
        final_trans = trans2 + final_trans
        final_trans = trans + final_trans
        final_trans.graph1 = original_graph  # Should not be needed ?!

        final_graph = final_trans.translate_subgraph(self)
        final_trans.graph2 = final_graph

        return final_graph, final_trans


    def _split_blocks_at_starts_and_ends(self, intervals):
        """
        Splits the region paths at the starts and ends of the
        intervals, such that all the intervals starts at the
        beginning of the intervals

        :param intervals: list(intervals)
        :returns: full translation for all block splits
        :rtype: Translation

        """
        full_translation = Translation()

        # Split starts
        for interval in intervals:
            offset = interval.start_position.offset
            if offset == 0:
                continue
            rp = interval.region_paths[0]
            full_translation += self._split_block(rp, [offset])

        for interval in intervals:
            end_offset = interval.end_position.offset
            rp = interval.end_position.region_path_id
            L = self.blocks[rp]
            if end_offset == L:
                continue
            full_translation += self.split_block(rp, [offset])

        return full_translation

    @takes(Interval, Interval)
    def merge_intervals(self, interval_a, interval_b, copy=True):
        """
        Merge the two intervals in the graph and return
        a translation object

        :param interval_a: One interval
        :param interval_b: One interval
        :returns: translation from old to new graph
        :rtype: Translation

        """
        assert interval_a.length() == interval_b.length()
        assert not any(rp_a in interval_b.region_paths
                       for rp_a in interval_a.region_paths)
        break_points = self._get_all_block_borders(
            interval_a,
            interval_b)

        sub_intervals_a = interval_a.split(break_points)
        sub_intervals_b = interval_b.split(break_points)

        ids = (self._next_id() for _ in sub_intervals_a)
        translation = self._generate_translation(
            interval_a, interval_b,
            sub_intervals_a, sub_intervals_b,
            ids)

        graph = self.copy() if copy else self
        # for pair in zip(sub_intervals_a, sub_intervals_b):
        #     start_intervals = [i for i in pair if i.start_position.offset == 0]
        #     in_edges = sum([self.reverse_adj_list[i.region_paths[0]]
        #                     for i in start_intervals], [])
        #     out_intervals = [
        #         i for i in pair if
        #         i.end_position.offset == self.blocks[i.end_position.region_path_id].length()]
        #     out_edges = sum([self.adj_list[i.region_paths[-1]] for i in end_intervals], [])
        # 
        #     graph

        return translation

    def _get_all_block_borders(self, interval_a, interval_b):
        """Return a list of block changes in both a and b

        :param interval_a: Interval
        :param interval_b: Interval
        :returns: list of all break_points
        :rtype: list

        """
        break_points = []
        for interval in [interval_a, interval_b]:
            offset = -interval.start_position.offset
            for region_path in interval.region_paths[:-1]:
                offset += self.blocks[region_path].length()
                break_points.append(offset)
            offset += interval.end_position.offset
            break_points.append(offset)

        unique_points = sorted(set(break_points))
        return unique_points

    @takes(Interval, Interval)
    def connect_intervals(self, interval_a, interval_b):
        pass

    def __str__(self):
        return "Graph: \n Blocks: %s\n Edges: %s" % (self.blocks, self.adj_list)

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def _get_reverse_edges(adj_list):
        reverse_edges = defaultdict(list)
        for block, edges in adj_list.items():
            for edge in edges:
                reverse_edges[edge].append(block)

        return reverse_edges

    def __eq__(self, other):

        for adj in self.adj_list:
            #if adj in other.adj_list and self.adj_list[adj] != other.adj_list[adj]:
            if set(self.adj_list[adj]) != set(other.adj_list[adj]):
                print("Different adj list for key %d" % (adj))
                return False
        for adj in other.adj_list:
            if set(self.adj_list[adj]) != set(other.adj_list[adj]):
            #if adj in self.adj_list and self.adj_list[adj] != other.adj_list[adj]:
                print("Different adj list 2")
                return False


        #if self.adj_list != other.adj_list:
        #    return False

        if self.blocks != other.blocks:
            print("Different blocks")
            return False

        return True

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
            if not match:
                # No match for block b, return False
                return False

        return True





