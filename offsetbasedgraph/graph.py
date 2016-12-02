from collections import defaultdict
from offsetbasedgraph import Interval
from util import takes


class Block(object):
    def __init__(self, length):
        self._length = length

    def length(self):
        return self._length

    def __eq__(self, other):
        return self.length == other.length


class Graph(object):
    region_paths = {}
    adj_list = defaultdict(list)
    reverse_adj_list = defaultdict(list)

    # Graph alterations
    def __init__(self, blocks, adj_list):
        self.blocks = blocks
        self.adj_list = adj_list
        self.reverse_adj_list = self._get_reverse_edges(adj_list)

    @takes(Interval, Interval)
    def merge_intervals(self, interval_a, interval_b):
        '''
        Merge the two intervals in the graph and return
        a translation object
        '''
        assert interval_a.length() == interval_b.length()

    @takes(Interval, Interval)
    def connect_intervals(self, interval_a, interval_b):
        pass

    @staticmethod
    def _get_reverse_edges(self, adj_list):
        pass
