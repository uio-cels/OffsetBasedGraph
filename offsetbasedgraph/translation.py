from .util import takes
from .interval import Interval, Position
from .multipathinterval import GeneralMultiPathInterval


class Translation(object):

    def __init__(self, translate_dict={}, reverse_dict={}):
        """
        Init the translation object with two dicts. Each dict has
        region path IDs as keys and a list of intervals as values.
        :param translate_dict: Mapping of region paths in graph1 to intervals
        in graph2
        :param reverse_dict: Mapping of region paths in graph2 to intervals
        in graph 1
        :return:
        """
        self._a_to_b = translate_dict
        self._b_to_a = reverse_dict

        # Store graph references
        self.graph1 = None
        self.graph2 = None
        if len(self._a_to_b) > 0:
            self.graph2 = list(self._a_to_b.values())[0][0].graph  # 1 st translation, 1st interval's graph
        if len(self._b_to_a) > 0:
            self.graph1 = list(self._b_to_a.values())[0][0].graph

        if self.graph1 is None:
            self.graph1 = self.graph2
        if self.graph2 is None:
            self.graph2 = self.graph1

    def _translations(self, region_path, inverse):
        if inverse:
            return self._b_to_a[region_path]
        else:
            return self._a_to_b[region_path]

    def _get_other_graph(self, inverse):
        if inverse:
            return self.graph1
        else:
            return self.graph2

    @takes(Interval)
    def translate_interval(self, interval, inverse=False):
        """
        Translate an interval between the two coordinate systems.
        :param interval:
        :param inverse:
        :return: Returns an interval. If inverse is True, a list of intervals
        is returned.
        """
        """
        Algorithm:
            Translate start to new start
            Translate end to new end
            Find all new region paths
        """

        new_starts = self.translate_position(interval.start_position, inverse)
        new_ends = self.translate_position(interval.end_position, inverse)
        new_region_paths = []

        for region_path in interval.region_paths:
            # Find all region paths that this region path follows
            intervals = self._translations(region_path, inverse)
            new_region_paths.extend([i.region_paths for i in intervals])

        new_interval = GeneralMultiPathInterval(new_starts, new_ends,
                            new_region_paths, self._get_other_graph(inverse))

        return new_interval

    @takes(Position)
    def translate_position(self, position, inverse=False):
        """
        Translates a position
        :param position:
        :param inverse: If True, translate back
        :return: Returns the translated Position
        """

        assert self.graph1 is not None, "Graph1 is none"
        assert self.graph2 is not None, "Graph2 is none"

        # Get interval for region path. Select first region path. Count offset.
        intervals = self._translations(position.region_path_id, inverse)
        positions = [interval.get_position_from_offset(position.offset)
                   for interval in intervals]

        return positions

    def __add__(self, other):
        """
        Combine (add) two translations.
        :param other: Another translation
        :return:
        """

        """
        Algorithm
            FOr every (region path => inter) in a, translate inter to c
            by using other.translate_interval(inter). Replace in a.
            For every (region path => inter) in c, translate recursviely back
            to all intervals using a.translate(inter, inverse). Replace in c
        """

        # Find new forward translation dict
        new_translate_dict = {}
        for t in self._a_to_b:
            new_translate_dict[t] = other.translate_interval(self._a_to_b[t])

        changed = True
        while(changed):
            changed = False
            for t in other._b_to_a:
                # For every block=>interval, map backwards to intervals
                new_intervals = []
                for inter in other._b_to_a[t]:
                    new_intervals.append(other.translate_interval(inter, True))
                new_intervals = set(new_intervals)
                if new_intervals != set(other._b_to_a[t]):
                    changed = True  # If change, translate again
                other._b_to_a[t] = new_intervals

    def __eq__(self, other):
        """Check if two translations are equal"""
        eq = self._a_to_b == other._a_to_b
        eq &= self._b_to_a == other._b_to_a
        return eq

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self._b_to_a)
        return "\n".join(["%s: %s" % (_id, ",".join(intervals))
                          for _id, intervals in
                          self._a_to_b.items()])
