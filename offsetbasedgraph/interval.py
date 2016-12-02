
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
        return str(list(self.region_path_id, self.offset))

class Interval(object):
    start_position = None
    end_position = None
    region_paths = None

    def length(self):
        pass

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

    def __str__(self):
        return "%d, %d, %s" % (self.start_position,
                               self.end_position, self.region_paths)

    