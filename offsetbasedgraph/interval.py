
class Position(object):
    """
    Represents a position  in the graph
    """
    def __init__(self, region_path_id, offset):
        self.region_path_id = region_path_id
        self.offset = offset


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
        r_lengths = [self.graph[rp].length() for rp in self.region_paths[:-1]]
        r_sum = sum(r_lengths)
        return r_sum-self.start_position.offset+self.end_position.offset

    def __init__(self, start_position, end_position, region_paths=None):
        self.start_position = start_position
        self.end_position = end_position

        # By default include start and end region path
        if region_paths is None:
            region_paths = [start_position.region_path_id]
            if end_position.region_path_id != start_position.region_path_id:
                region_paths.append(end_position.region_path_id)

        self.region_paths = region_paths

    def __eq__(self, other):
        eq = self.start_position == other.start_position
        eq &= self.end_position == other.end_position
        eq &= self.region_paths == other.region_paths
        return eq

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
                return Position(rp, total_offset)
            total_offset -= rp_length

        assert False
