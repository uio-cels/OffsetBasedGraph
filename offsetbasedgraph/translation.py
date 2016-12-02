from .util import takes
from .interval import Interval, Position


class Translation(object):
    _a_to_b = {}  # Region paths to interval
    _b_to_a = {}

    def __init__(self, translate_dict, reverse_dict):
        """
        Init the translation object with two dicts. Each dict has
        region path IDs as keys and a list of intervals as values.
        :param translate_dict: The forward translation dict
        :param reverse_dict: The reverse translation dict
        :return:
        """
        self._a_to_b = translate_dict
        self._b_to_a = reverse_dict

    def _translations(self, region_path, inverse):
        if inverse:
            return self._b_to_a[region_path]
        else:
            return self._a_to_b[region_path]

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
            new_region_paths.append(intervals)

        new_intervals = []
        for start, end, region_paths in zip(new_starts, new_ends, new_region_paths):
            new_interval = Interval(start, end, region_paths)
            new_intervals.append(new_interval)

        return new_intervals

    @takes(Position)
    def translate_position(self, position, inverse=False):
        """
        Translates a position
        :param position:
        :param inverse: If True, translate back
        :return: Returns the translated Position
        """
        # Get interval for region path. Select first region path. Count offset.
        intervals = self._translations(position.region_path_id, inverse)
        positions = [interval[0].get_position_fron_offset(position.offset)
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
