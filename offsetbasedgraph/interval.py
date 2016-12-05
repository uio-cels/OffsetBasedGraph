
class Position(object):
    """
    Represents a position  in the graph
    """
    def __init__(self, region_path_id, offset):
        self.region_path_id = region_path_id
        self.offset = offset

    def __eq__(self, other):
        return self.region_path_id == other.region_path_id \
               and self.offset == other.offset

    def __str__(self):
        return "(%s:%s)" % (self.region_path_id, self.offset)


class Interval(object):
    start_position = None
    end_position = None
    region_paths = None
    graph = None

    def length(self):
        """

        :returns: The length of the interval
        :rtype: int

        """
        r_lengths = [self.graph.blocks[rp].length()
                     for rp in self.region_paths[:-1]]

        r_sum = sum(r_lengths)
        return r_sum-self.start_position.offset+self.end_position.offset

    def split(self, offsets):
        """Split the interval at the given offsets

        :param offsets: list of int
        :returns: list of intervals
        :rtype: list

        """
        split_points = [self.get_position_from_offset(offset) for offset
                        in offsets]
        region_path_idxs = {rp: i for i, rp in enumerate(self.region_paths)}
        split_intervals = []

        starts = [self.start_position] + split_points
        ends = split_points + [self.end_position]

        for start_pos, end_pos in zip(starts, ends):
            start_idx = region_path_idxs[start_pos.region_path_id]
            end_idx = region_path_idxs[end_pos.region_path_id]
            region_paths = self.region_paths[start_idx:end_idx+1]
            split_intervals.append(Interval(start_pos, end_pos, region_paths))

        return split_intervals

    def join(self, other):
        """Create new interval as union of self and other
        other must start at self.end_position

        :param other: Interval
        :returns: joined interval
        :rtype: Interval

        """

        assert self.end_position == other.start_position
        region_paths = self.region_paths[:-1] + other.region_paths
        return Interval(self.start_position, other.end_position, region_paths)

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None):

        self.start_position = start_position
        self.end_position = end_position

        # By default include start and end region path
        if region_paths is None:
            region_paths = [start_position.region_path_id]
            if end_position.region_path_id != start_position.region_path_id:
                region_paths.append(end_position.region_path_id)

        self.region_paths = region_paths
        self.graph = graph

    def __eq__(self, other):
        eq = self.start_position == other.start_position
        eq &= self.end_position == other.end_position
        eq &= self.region_paths == other.region_paths
        return eq

    def __str__(self):
        return "%s, %s, %s" % (self.start_position,
                               self.end_position, self.region_paths)

    def get_position_from_offset(self, offset):
        """Get position of with offset counted from the start of
        the interval

        :param offset:
        :returns: The position in the graph
        :rtype: Position

        """
        assert offset < self.length()
        total_offset = offset + self.start_position.offset
        for region_path in self.region_paths:
            rp = self.graph.blocks[region_path]
            rp_length = rp.length()
            if rp_length > total_offset:
                return Position(region_path, total_offset)
            total_offset -= rp_length

        assert False


def interval_factory(graph):
    """Create function that initializes
    intervals without specifying the graph
    every time

    :param graph: The graph for the intervals
    :returns: initializer for Interval
    :rtype: func

    """
    def get_interval(*args, **kwargs):
        """Initialize an interval with the specified graph

        :returns: 
        :rtype: Interval

        """
        kwargs["graph"] = graph
        return Interval(*args, **kwargs)

    return get_interval
