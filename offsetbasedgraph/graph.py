from collections import defaultdict
from .interval import Interval, Position
from .util import takes
from .translation import Translation


class Block(object):
    def __init__(self, length):
        self._length = length

    def length(self):
        return self._length

    def __eq__(self, other):
        return self.length() == other.length()


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

    @staticmethod
    def from_file(file_name):
        """
        Reads from file and creates the graph
        :param file_name: File name
        :return:
        """
        pass

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

    def _join_blocks(self, block_ids):
        """Merge the given blocks in the graph and
        giving them the combination of all edges.

        :param block_ids: block ids of blocks to merge
        :returns: translation object
        :rtype: Translation

        """
        if not block_ids:
            return Translation()
        lens = [len(self.blocks[_id] for _id in block_ids)]
        L = lens[0]
        assert all(l == L for l in lens)
        in_edges = set()
        out_edges = set()
        for block_id in block_ids:
            in_edges.update(self.reverse_adj_list[block_id]
                            for block_id in block_ids)
            out_edges.update(self.adj_list[block_id]
                             for block_id in block_ids)
        for block_id in block_ids:
            self.remove(block_id)

        new_id = self._next_id()
        translation = Translation(
            {block_id: Interval(Position(new_id, 0), Position(new_id, L))
             for block_id in block_ids},
            {new_id:
             [Interval(Position(block_id, 0), Position(block_id, L))
              for block_id in block_ids]})
        return translation

    def _split_block(self, block_id, offset):
        """Splits the given block at offset
        and return translation object

        :param block_id: block to split
        :param offset: where to split
        :returns: translation object
        :rtype: Translation

        """
        l = self.blocks[block_id].length()
        if offset == 0 or offset >= l-1:
            return Translation()
        id_a = self._next_id()
        id_b = self._next_id()
        self.blocks[id_a] = Block(offset)
        self.blocks[id_b] = Block(l-offset)

        self.adj_list[id_a] = id_b
        self.reverse_adj_list[id_b] = id_a

        self.adj_list[id_b] = self.adj_list[block_id]
        self.reverse_adj_list[id_a] = self.reveres_adj_list[block_id]

        self.remove(block_id)

        return Translation(
            {block_id:
             Interval(Position(id_a, 0), Position(id_b, l-offset))},
            {id_a: Interval(Position(block_id, 0),
                            Position(block_id, offset)),
             id_b: Interval(Position(block_id, offset),
                            Position(block_id, l))}
            )

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

        for pair in zip(sub_intervals_a, sub_intervals_b):
            start_intervals = [i for i in pair if i.start_position.offset == 0]
            in_edges = sum([self.reverse_adj_list[i.region_paths[0]]
                            for i in start_intervals], [])
            out_intervals = [
                i for i in pair if
                i.end_position.offset == self.blocks[i.end_position.region_path_id].length()]
            out_edges = sum([self.adj_list[i.region_paths[-1]] for i in end_intervals], [])

            graph

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

    @staticmethod
    def _get_reverse_edges(adj_list):
        reverse_edges = defaultdict(list)
        for block, edges in adj_list.items():
            for edge in edges:
                reverse_edges[edge].append(block)

        return reverse_edges

    def __eq__(self, other):
        if self.adj_list != other.adj_list:
            return False

        if self.blocks != other.blocks:
            return False

        return True
