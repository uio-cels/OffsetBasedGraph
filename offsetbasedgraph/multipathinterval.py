from itertools import chain
from .interval import Interval, Position


class MultiPathInterval(object):

    def get_single_path_intervals(self):
        raise NotImplementedError()


class SingleMultiPathInterval(MultiPathInterval):
    def __init__(self, interval, graph=None):
        self.interval = interval.copy()
        if graph is not None:
            self.interval.graph = graph
        self.graph = graph
        rps = interval.region_paths
        if len(rps) > 1 and rps[-1] == rps[0]:
            self.interval.region_paths = rps[:-1]

    def get_single_path_intervals(self):
        return [self.interval]


class CriticalPathsMultiPathInterval(MultiPathInterval):
    def __init__(self, start_pos, end_pos, critical_intervals):
        """
        :param start_pos: start Position
        :param end_pos: end Position
        :param critical_intervals: list of critical
            intervals that the multipath interval needs to go through
        """

        self.start_pos = start_pos
        self.end_pos = end_pos
        self.critical_intervals = critical_intervals

    def equal_critical_intervals(self, other):
        for i in self.critical_intervals:
            if i not in other.critical_intervals:
                return False
        for i in other.critical_intervals:
            if i not in self.critical_intervals:
                return False

        return True

    def faster_equal_critical_intervals(self, other):
        # Assuming all critical intervals are ordered and corresponding
        # in the two mps
        if len(self.critical_intervals) != len(other.critical_intervals):
            return False

        for i in range(0, len(self.critical_intervals)):
            if self.critical_intervals != other.critical_intervals:
                return False

        return True

    def __eq__(self, other):
        """
        NB: Only compare exact representation. Will give wrong result if equal,
        but different critical_intervals.
        TODO: Fix so that it compares the minimum representation
        """
        if other.start_pos != self.start_pos:
            return False
        if other.end_pos != self.end_pos:
            return False
        return self.faster_equal_critical_intervals(other)

    def str_without_critical(self):
        out = "%s, %s" % (str(self.start_pos), str(self.end_pos))
        return out

    def __str__(self):
        out = "Start: " + str(self.start_pos) + "\n"
        out += "Critical paths:\n"
        for c in self.critical_intervals:
            out += "  " + str(c) + "\n"
        out += "End: " + str(self.end_pos) + "\n"
        return out

    def __repr__(self):
        return str(self)


class FuzzyMultipathInterval(MultiPathInterval):
    def __init__(self, start_offset, end_offset, region_paths,
                 threshold, graph=None):
        self.start_position = Position(region_paths[0], start_offset)
        self.end_position = Position(region_paths[-1], end_offset)
        self.region_paths = region_paths
        self.threshold = threshold
        self.graph = graph

    @classmethod
    def from_interval(cls, interval, threshold):
        return cls(interval.start_position.offset,
                   interval.end_position.offset,
                   interval.region_paths, threshold,
                   interval.graph)

    def get_single_path_intervals(self):
        raise NotImplementedError()

    def __eq__(self, other):
        """Check if other deviates from self
        with less than self.threshold bps

        :param other: other interval
        :returns: wheter they are equal
        :rtype: bool

        """
        if not self.start_position == other.start_position:
            return False
        if not self.end_position == other.end_position:
            return False

        comparison_directions = [(self, other), (other, self)]
        for template, test_interval in comparison_directions:
            cur_deviance = 0
            for region_path in test_interval.region_paths:
                if region_path in template.region_paths:
                    cur_deviance = 0
                    continue
                cur_deviance += self.graph.blocks[region_path].length()
                if cur_deviance >= self.threshold:
                    return False

        return True

    def _calculate_length(self):
        if len(self.region_paths) == 1:
            return self.end_position.offset - self.start_position.offset

        if not self.region_paths:
            return 0

        assert self.graph is not None,\
            "Graph is none and length cache is None and\
            interval has more than 1 rp. Cannot compute length"

        rp_lengths = [self.graph.blocks[rp].length()
                      for rp in self.region_paths[:-1]]
        r_sum = sum(rp_lengths)
        length = r_sum-self.start_position.offset+self.end_position.offset

        assert length >= 0,\
            "Length is %d for interval %s. r_lengths: %s. Graph: %s"\
            % (length, self, rp_lengths, self.graph)

        return length

    def length(self):
        """

        :returns: The length of the central path of interval
        :rtype: int

        """

        if hasattr(self, "length_cache") and self.length_cache is not None:
            return self.length_cache
        length = self._calculate_length()
        self.length_cache = length
        return length

    def __str__(self):
        return "%s\t%s\t%s" % (self.start_position.offset,
                               self.end_position.offset,
                               " ".join(str(c) for c in self.region_paths))

    __repr__ = __str__


class SimpleMultipathInterval(object):

    def __init__(self, intervals):
        self._intervals = intervals

    def get_single_path_intervals(self):
        return self._intervals


class GeneralMultiPathInterval(MultiPathInterval):
    """
    Holds a general multipath interval by representing all possible start
    positions, end positions and all possible region paths within the interval.

    A method can return (in some cases) a trivial set of single path
    intervals contained by this multi path interval
    """

    def __init__(self, start_positions, end_positions, region_paths, graph):
        self.start_positions = start_positions
        self.end_positions = end_positions
        self.region_paths = region_paths
        self.end_region_paths = {pos.region_path_id: pos
                                 for pos in self.end_positions}
        self.graph = graph

    def is_single_path(self):
        pass

    def _find_end_region_path(self, start_position, region_paths):
        """Generate all valid intervals starting at start_position
        going through region_paths

        :param start_position: start position of the intervals
        :param region_paths: list of region paths
        :returns: generator of valid intervals
        :rtype: generator

        """
        cur_region_path = region_paths[-1]
        if cur_region_path in self.end_region_paths:
            yield Interval(start_position,
                           self.end_region_paths[cur_region_path],
                           region_paths, self.graph)

        next_level = (
            rp
            for rp in self.graph.adj_list[cur_region_path]
            if rp in self.region_paths)

        for rp in next_level:
            for end in self._find_end_region_path(
                    start_position, region_paths + [rp]):
                yield end

    def _get_single_path_intervals(self, start_position):
        """Generate all valid intervals starting from start position

        :param start_position: Position
        :returns: Interval generator
        :rtype: generator

        """
        for end in self._find_end_region_path(start_position,
                                              [start_position.region_path_id]):
            yield end

    def get_single_path_intervals(self):
        """
        Generate all single path intervals.

        :return: generator
        """
        intervals = []
        for start_position in self.start_positions:
            # yield from self._get_single_path_intervals(start_position)
            intervals.extend(self._get_single_path_intervals(start_position))
        return intervals

    def __str__(self):
        return "Start pos: %s\n end pos: %s\n region paths: %s\n" % \
               (str(self.start_positions), str(self.end_positions), str(self.region_paths))
