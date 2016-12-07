
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

    def __repr__(self):
        return self.__str__()

class Interval(object):

    def length(self):
        """

        :returns: The length of the interval
        :rtype: int

        """
        r_lengths = [self.graph.blocks[rp].length()
                     for rp in self.region_paths[:-1]]

        r_sum = sum(r_lengths)
        return r_sum-self.start_position.offset+self.end_position.offset

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None):

        if isinstance(start_position, int):
            self.start_position = Position(region_paths[0], start_position)
            print("Start pos %s" % self.start_position)
        else:
            self.start_position = start_position

        if isinstance(end_position, int):
            self.end_position = Position(region_paths[-1], end_position)
        else:
            self.end_position = end_position

        # By default include start and end region path
        if region_paths is None:
            region_paths = [start_position.region_path_id]
            if end_position.region_path_id != start_position.region_path_id:
                region_paths.append(end_position.region_path_id)

        self.region_paths = region_paths
        self.graph = graph

        # Sanity check interval
        for rp in region_paths:
            assert rp in graph.blocks, "Region path %s not in graph" % rp

        assert self.start_position.region_path_id in self.region_paths
        assert self.end_position.region_path_id in self.region_paths

        # Check offsets
        assert self.end_position.offset <= graph.blocks[self.region_paths[-1]].length(), \
                "Offset %d in block %d with size %d" % (self.end_position.offset,
                                            self.region_paths[-1],
                                            graph.blocks[self.region_paths[-1]].length())

    def __eq__(self, other):
        eq = self.start_position == other.start_position
        eq &= self.end_position == other.end_position
        eq &= self.region_paths == other.region_paths
        return eq

    def __str__(self):
        return "%s, %s, %s" % (self.start_position,
                               self.end_position, self.region_paths)

    def __repr__(self):
        return self.__str__()

    def get_position_from_offset(self, offset):
        """Get position of with offset counted from the start of
        the interval

        :param offset:
        :returns: The position in the graph
        :rtype: Position

        """
        assert offset <= self.length(), \
            "Offset %d is larger than total length of interval %d in interval %s" \
            % (offset, self.length(), self.__str__())
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
