
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

    def length(self):
        pass

    def __init__(self, start_position, end_position, region_paths):
        self.start_position = start_position
        self.end_position = end_position
        self.region_paths = region_paths

    def __eq__(self, other):
        eq = self.start_position == other.start_position
        eq &= self.end_position == other.end_position
        eq &= self.region_paths == other.region_paths
        return eq
