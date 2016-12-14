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

        # new_first_intervals = [trans.translate_rp(rp) for 

    def _get_insulated_merge_transformation(self, intervals):
        pass

    def _get_inslulate_transformation(self, intervals):
        # Split starts
        trans_dict = {}
        reverse_dict = {}
        starts = [i.start_position for i in intervals]
        for start in starts:
            if start.offset == 0:
                continue
            rp = start.region_path_id
            offset = start.offset
            id_a, id_b = (self._next_id(), self.next_id())
            L = self.blocks[rp].length()
            trans_dict[rp] = [Interval(
                Position(id_a, 0),
                Position(id_b, L-offset))]
            reverse_dict[id_a] = [Interval(Position(rp, 0),
                                           Position(rp, offset))]
            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L))]
        

    def get_merge_transformation(self, intervals):
        trans = self._get_inslulate_transformation(intervals)
        trans2 = self._get_insulated_merge_transformation(
            [trans.translate(interval) for interval in intervals])
        return trans+trans2

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

        print(interval_a)
        print(break_points)
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
        return "Blocks: %s\Edges: %s" % (self.blocks, self.adj_list)

    @staticmethod
    def _get_reverse_edges(adj_list):
        reverse_edges = defaultdict(list)
        for block, edges in adj_list.items():
            for edge in edges:
                reverse_edges[edge].append(block)

        return reverse_edges

    def __eq__(self, other):

        for adj in self.adj_list:
            if adj in other.adj_list and self.adj_list[adj] != other.adj_list[adj]:
                return False
        for adj in other.adj_list:
            if adj in self.adj_list and self.adj_list[adj] != other.adj_list[adj]:
                return False


        #if self.adj_list != other.adj_list:
        #    return False

        if self.blocks != other.blocks:
            return False

        return True
