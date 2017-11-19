from offsetbasedgraph import Interval
import numpy as np


class DirectedInterval(Interval):

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None, direction=None):
        super(DirectedInterval, self).__init__(start_position, end_position,
                                               region_paths, graph)

    def contains_rp(self, region_path_id):
        return -region_path_id in self.region_paths or \
               region_path_id in self.region_paths

    def contains_position(self, position):
        if not self.contains_rp(position.region_path_id):
            return False

        if -position.region_path_id in self.region_paths[1:-1] \
                or position.region_path_id in self.region_paths[1:-1]:
            return True

        if self.graph is None:
            raise Exception("Interval's graph is None. Graph is needed to check if position is in start or end rp of intnerval.")

        block_length = self.graph.blocks[position.region_path_id].length()
        offset = position.offset

        if abs(position.region_path_id) == abs(self.region_paths[0]):
            if np.sign(position.region_path_id) != np.sign(self.start_position.region_path_id):
                offset = block_length - position.offset

            if offset < self.start_position.offset:
                    return False

        if abs(position.region_path_id) == abs(self.region_paths[-1]):
            if np.sign(position.region_path_id) != np.sign(self.end_position.region_path_id):
                offset = block_length - position.offset

            if offset >= self.end_position.offset:
                return False

        return True

    def can_be_on_negative_strand(self):
       return self._can_be_on_strand(plus_strand=False)

    def can_be_on_positive_strand(self):
        return self._can_be_on_strand(plus_strand=True)

    def _can_be_on_strand(self, plus_strand=True):
        if len(self.region_paths) == 1:
            return True

        if self.graph is None:
            raise Exception("Graph cannot be None when ckecinng 'can be on strand'")

        for i in range(len(self.region_paths) - 1):
            from_rp = self.region_paths[i]
            to_rp = self.region_paths[i+1]
            #print("Checking from %d to %d, positive? %d" % (from_rp, to_rp, plus_strand))

            if to_rp in self.graph.adj_list[from_rp]:
                if plus_strand:
                    return True

            if to_rp in self.graph.reverse_adj_list[from_rp]:
                if not plus_strand:
                    return True

        return False

    def get_subinterval(self, start_offset, end_offset):
        subinterval = super(DirectedInterval, self).get_subinterval(start_offset, end_offset)
        return DirectedInterval(subinterval.start_position,
                                subinterval.end_position,
                                subinterval.region_paths,
                                self.graph,
                                self.direction)

    def get_reverse(self):
        reverse = super(DirectedInterval, self).get_reverse()
        return DirectedInterval(reverse.start_position,
                                reverse.end_position,
                                reverse.region_paths,
                                self.graph,
                                reverse.direction)

    def intersects(self, other):
        raise NotImplementedError()

    def get_position_from_offset(self, offset):
        raise NotImplementedError()

    def get_adj_list(self):
        raise NotImplementedError()

    def copy(self):
        raise NotImplementedError()
