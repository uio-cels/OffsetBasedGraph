

class GeneralMultiPathInterval(object):
    """
    Holds a general multipath interval by representing all possible start
    positions, end positions and all possible region paths within the interval.

    A method can return (in some cases) a trivial set of single path
    intervals contained by this multi path interval
    """

    def __init__(self, start_positions, end_position, region_paths, graph):
        self.start_positions = start_positions
        self.end_positions = end_position
        self.region_paths = region_paths
        self.graph = graph

    def is_single_path(self):
        pass

    def get_single_path_intervals(self):
        """
        Returns a list of all single path intervals.
        :return:
        """


