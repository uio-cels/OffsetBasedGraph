from collections import defaultdict


class Position(object):
    """
    Represents a position  in the graph
    """
    region_path_id = None
    offset = None


class Interval(object):
    start_position = None
    end_position = None
    region_paths = None

    def length(self):
        pass

    def __init__(self, start_position, end_position, region_paths):
        self.start_position = start_position
        self.end_position = end_position
        self.region_paths = region_paths
        self.rp_lens_tmp = None  # For caching of region path lengths

class Translation(object):
    _a_to_b = {}  # Region paths to interval
    _b_to_a = {}

    def __init__(self, translate_dict, reverse_dict):
        pass

    @takes(Interval)
    def translate_interval(self, interval, inverse=False):
        '''
        Return the interval representation of @interval in
        coordinate system b
        '''
        for region_path in interval.region_paths:
            pass

    @takes(Position)
    def translate_position(self, position, inverse=False):
        pass

    def __add__(self, other):
        """Combine two translations"""


class Block(object):
    def __init__(self, length):
        self._length = length

    def length(self):
        return self._length


class Graph(object):
    blocks = {}
    edge_list = defaultdict(list)
    reverse_edge_list = defaultdict(list)

    # Graph alterations
    def __init__(self, blocks, edge_list):
        self.blocks = blocks
        self.edge_list = edge_list
        self._generate_reverse_edges()

    @classmethod
    def from_file(cls, filename):
        pass

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

