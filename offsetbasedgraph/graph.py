from collections import defaultdict
from .interval import Interval, Position
from .util import takes
from .translation import Translation


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

        return "Block(%s)" % self._length

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
            print("#", rp)
            if self.blocks[rp].length() == offset:
                continue
            trans, graph = graph.get_split_translation(rp, offset)
            translation += trans

        print(graph)
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
        # Split starts
        trans_dict = {}
        reverse_dict = {}
        starts = [i.start_position for i in intervals]
        translation = Translation(graph=self)
        cur_graph = self
        i = 0

        for start in starts:
            i += 1
            if start.offset == 0:
                continue
            rp = start.region_path_id
            offset = start.offset
            id_a, id_b = (self._next_id(), self._next_id())
            L = self.blocks[rp].length()

            trans_dict[rp] = [Interval(
                Position(id_a, 0),
                Position(id_b, L-offset))]

            reverse_dict[id_a] = [Interval(Position(rp, 0),
                                           Position(rp, offset),
                                           graph=cur_graph)]

            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L),
                                           graph=cur_graph)]

            tmp_trans = Translation(trans_dict, reverse_dict, graph=cur_graph)
            cur_graph = tmp_trans.translate_subgraph(cur_graph)
            for i_list in trans_dict.values():
                for j in i_list:
                    j.graph = cur_graph

            translation += tmp_trans

        new_intervals = [translation.translate(i) for i in intervals]
        ends = [i.end_position for i in new_intervals]

        for end in ends:
            assert cur_graph is not None
            rp = end.region_path_id
            L = self.blocks[rp].length()
            offset = end.offset
            if offset == L:
                continue
            id_a, id_b = (self._next_id(), self._next_id())
            trans_dict[rp] = [Interval(
                Position(id_a, 0),
                Position(id_b, L-offset)
                )]
            reverse_dict[id_a] = [Interval(Position(rp, 0),
                                           Position(rp, offset),
                                           graph=cur_graph)]
            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L),
                                           graph=cur_graph)]
            tmp_trans = Translation(trans_dict, reverse_dict, graph=cur_graph)
            cur_graph = tmp_trans.translate_subgraph(cur_graph)
            for i_list in trans_dict.values():
                for i in i_list:
                    i.graph = cur_graph
            tmp_trans.graph2 = cur_graph
            translation += tmp_trans

        return translation, cur_graph

    def get_merge_translation(self, intervals):
        """Get translation object for merging the
        given intervals in the graph

        :param intervals: list of intevals to merge
        :returns: translation object
        :rtype: Translation

        """
        trans, ng = self._get_inslulate_translation(intervals)
        trans2, ng = ng._get_insulated_merge_translation(
            trans.translate_interval(interval) for interval in intervals)

    def _update_a_b_graph(self, a_b, graph):
        for intvs in a_b.values():
            for interval in intvs:
                interval.graph = graph

    def get_merge_translation2(self, intervals):
        original_graph = intervals[0].graph

        # Assume same lengths for all intervals
        length = intervals[0].length()
        for interval in intervals[1:]:
            assert interval.length() == length, \
                "All intervals should have the same length (%d, %d)" % (interval.length(), length)

        # Step 1: Translate graph so that all intervals starts and ends at region paths
        trans, graph1 = self._get_inslulate_translation(intervals)
        #self._update_a_b_graph(trans._a_to_b, graph1)
        trans.graph2 = graph1
        print("trans1 a to b graph")
        self._update_a_b_graph(trans._a_to_b, graph1)  # correct, a to b interval's graph is wrong for some reason
        print(list(trans._a_to_b.values())[0][0].graph)

        print("=== Graph 1 === ")
        print(graph1)

        print("=== Trans 1 ===")
        print(trans.graph1)
        print(trans)

        # Update intervals to be on graph1
        new_intervals = []
        for interval in intervals:
            mp = trans.translate_interval(interval)
            print("MP" , str(mp))
            new_interval = mp.get_single_path_intervals()[0]
            new_interval.graph = graph1
            new_intervals.append(new_interval)

        intervals = new_intervals
        print("=== Intervals ===")
        print(intervals[0])
        print(intervals[1])

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

        print("a_b", a_b)
        print("b_a", b_a)

        # Find new graph
        trans2 = Translation(a_b, b_a, graph1)
        graph2 = trans2.translate_subgraph(graph1)
        trans2.graph2 = graph2


        # Update graph in a_b intervals to graph2 (which is now created)
        self._update_a_b_graph(a_b, graph2)

        print("=== Graph 2 === ")
        print(graph2)

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

        print("=== Graph 3 === ")
        print(graph3)

        # Step 4: Divide large middle block

        starts = []
        for interval in intervals:
            prev_start = 0
            for rp in interval.region_paths:
                prev_start += interval.graph.blocks[rp].length()
                starts.append(prev_start)

        starts.sort()
        starts = list(set(starts))[:]  # remove last, becauase end
        print(starts)
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
            prev_start = start
            last_length = start - prev_start

        a_b[new_block].append(Interval(0, last_length, new_blocks))

        trans4 = Translation(a_b, b_a, graph3)
        graph4 = trans4.translate_subgraph(graph3)
        self._update_a_b_graph(a_b, graph4)
        trans4.graph2 = graph4

        print("=== Graph 4 === ")
        print(graph4)

        print("=== trans 4 ===")
        print(trans4)
        print("=== trans 3 ===")
        print(trans3)
        final_trans = trans3 + trans4
        print("=== Final trans 3 + 4===")
        print(final_trans)
        final_trans = trans2 + final_trans
        #print(trans2.graph1)
        #print(trans2._b_to_a)
        print("=== Final trans ===")
        print(final_trans)
        print("first trans")
        print(trans)
        final_trans = trans + final_trans

        print("=== Final! trans ===")
        print(final_trans)

        final_trans.graph1 = original_graph  # Should not be needed ?!


        print(list(final_trans._a_to_b.values())[0][0].graph)
        print(list(final_trans._b_to_a.values())[0][0].graph)

        return graph4, final_trans




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
            if self.adj_list[adj] != other.adj_list[adj]:
                return False
        for adj in other.adj_list:
            if self.adj_list[adj] != other.adj_list[adj]:
            #if adj in self.adj_list and self.adj_list[adj] != other.adj_list[adj]:
                return False


        #if self.adj_list != other.adj_list:
        #    return False

        if self.blocks != other.blocks:
            return False

        return True
