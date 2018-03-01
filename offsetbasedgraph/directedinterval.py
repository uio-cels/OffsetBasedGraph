from offsetbasedgraph import Interval
import numpy as np


class DirectedInterval(Interval):

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None, direction=None):
        super(DirectedInterval, self).__init__(start_position, end_position,
                                               region_paths, graph)

