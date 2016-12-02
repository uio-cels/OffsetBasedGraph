from util import takes


class Translation(object):
    _a_to_b = {}  # Region paths to interval
    _b_to_a = {}

    def __init__(self, translate_dict, reverse_dict):
        self._a_to_b = translate_dict
        self._b_to_a = reverse_dict

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

    def __eq__(self, other):
        """Check if two translations are equal"""
        eq = self._a_to_b == other._a_to_b
        eq &= self._b_to_a == other._b_to_a
        return eq
