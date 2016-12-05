from collections import defaultdict
from .interval import Interval
from .util import takes


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

    @staticmethod
    def from_file(file_name):
        """
        Reads from file and creates the graph
        :param file_name: File name
        :return:
        """
        pass

    @takes(Interval, Interval)
    def merge_intervals(self, interval_a, interval_b):
        """
        Merge the two intervals in the graph and return
        a translation object

        :param interval_a: One interval
        :param interval_b: One interval
        :returns: translation from old to new graph
        :rtype: Translation

        """
        assert interval_a.length() == interval_b.length()
        # Create translation object
        # # Left side:
        break_points = self._get_all_block_borders(
            interval_a,
            interval_b)

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

        unique_points = list(set(break_points))
        unique_points.sort()
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
