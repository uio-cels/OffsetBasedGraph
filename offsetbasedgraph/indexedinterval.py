from .interval import Interval, Position
from collections import defaultdict
import math

class IndexedInterval(Interval):
    """
    Subclass of interval that indexes distance to nodes, making it faster to find subinterval
    """

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None, direction=1):
        super(IndexedInterval, self).__init__(
            start_position, end_position, region_paths, graph, direction
        )

        self.distance_index = defaultdict(list)  # List of rp number (pos in rp list) where end is at that distance or after (making it safe to start traversing at one of those nodes if finnding subinnterval at bigger thann that offset)
        self.distance_to_node = {}

        self._create_index()
        print("Distance index")
        print(self.distance_index)
        print("Distance to node")
        print(self.distance_to_node)

    def _create_index(self):
        nodes = []
        current_offset = 0
        i = 0
        prev_index_created = 0
        self.distance_index[0].append(0)
        for rp in self.region_paths:
            node_size = self.graph.node_size(rp)
            length_to_end = node_size
            if i == 0:
                length_to_end -= self.start_position.offset
                #current_offset += self.start_position.offset

            index = math.floor((current_offset + length_to_end) / 4) * 4

            for set_index in range(prev_index_created, index + 4, 4):
                self.distance_index[set_index].append(i)

            prev_index_created = index + 4

            self.distance_to_node[rp] = current_offset
            i += 1
            current_offset += length_to_end




    def position_at_offset(self, offset):
        assert offset >= 0, "Offset %d is negative" % offset
        assert offset <= self.length(), "Offset %d > interval length %d" % (offset, self.length())
        assert self.graph is not None

        current_offset = 0
        i = 0

        # Find i using index
        index_pos = math.floor(offset / 4) * 4
        print("Index pos for offset %d: %d" % (offset, index_pos))
        print("Index: %d" % index_pos)
        node_start = self.distance_index[index_pos][-1]  # Safe to choose last node
        i = node_start

        for rp in self.region_paths[i:]:
            current_offset = self.distance_to_node[rp]

            node_size = self.graph.node_size(rp)
            length_to_end = node_size
            if i == 0:
                length_to_end -= self.start_position.offset
            i += 1

            if current_offset + length_to_end > offset:
                return Position(rp, (node_size - length_to_end - current_offset) + offset)

            current_offset += length_to_end

    def _nodes_between_offets(self, start, end):
        nodes = []
        current_offset = 0
        i = 0

        # Find i using index
        index_pos = math.floor(start / 4) * 4
        node_start = self.distance_index[index_pos][-1]  # Safe to choose last node
        i = node_start

        for rp in self.region_paths[i:]:
            current_offset = self.distance_to_node[rp]
            node_size = self.graph.node_size(rp)
            length_to_end = node_size
            if i == 0:
                length_to_end -= self.start_position.offset
                #current_offset += self.start_position.offset
            i += 1
            if start < current_offset + length_to_end and end > current_offset:
                #print("    Rp %d, start: %d, end: %d, current offset: %d" % (rp, start, end, current_offset))
                #print("   adding %d" % rp)
                nodes.append(rp)

            if end < current_offset:
                break
            current_offset += length_to_end

        return nodes


