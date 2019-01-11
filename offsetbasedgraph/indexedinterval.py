from .interval import Interval, Position
from collections import defaultdict
import math
import logging
import pickle
from .graph import BlockArray
import numpy as np
from .directedinterval import DirectedInterval


class IndexedInterval(Interval):
    """
    Subclass of interval that indexes distance to nodes,
    making it faster to find subinterval
    """

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None, direction=1,
                 only_create_distance_to_nodes_index=False):
        assert graph is not None, "Index interval requires graph"

        super(IndexedInterval, self).__init__(
            start_position, end_position, region_paths, graph, direction
        )

        self.distance_index = defaultdict(list)
        # List of rp number (pos in rp list) where end is at that distance
        # or after (making it safe to start traversing at one of those nodes
        # if finnding subinnterval at bigger thann that offset)
        self.distance_to_node = {}
        self._create_index(only_create_distance_to_nodes_index)

    def to_file(self, file_name):
        with open(file_name, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def from_file(cls, file_name):
        f = open(file_name, "rb")
        interval = pickle.load(f)
        f.close()
        return interval

    def _create_index(self, only_create_distance_to_nodes=False):
        logging.info("Creating indexes for indexed interval")
        current_offset = 0
        i = 0
        prev_index_created = 0
        self.distance_index[0].append(0)
        for rp in self.region_paths:
            node_size = self.graph.node_size(rp)
            length_to_end = node_size
            if i == 0:
                length_to_end -= self.start_position.offset

            if not only_create_distance_to_nodes:
                index = math.floor((current_offset + length_to_end) / 4) * 4

                for set_index in range(int(prev_index_created), int(index + 4), 4):
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
        #print("Index pos for offset %d: %d" % (offset, index_pos))
        #print("Index: %d" % index_pos)
        node_start = self.distance_index[index_pos][-1]
        # Safe to choose last node
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

    def get_offset_at_position(self, position, direction="+"):
        # Return number of basepairs to this position
        if direction == "-":
            return self.get_reverse_offset_at_position(position)

        distance = self.distance_to_node[position.region_path_id] + position.offset
        if position.region_path_id == self.region_paths[0]:
            distance -= self.start_position.offset
        return int(distance)

    def get_reverse_offset_at_position(self, position, is_end_pos=True):
        rp = -position.region_path_id
        offset = self.graph.node_size(rp)-position.offset
        distance = self.distance_to_node[rp] + offset
        if rp == self.region_paths[0]:
            distance -= self.start_position.offset

        if not is_end_pos:
            distance -= 1

        return int(distance)


class NodeInterval(Interval):
    """
    Simple interval where start and end position is non-important.
    """
    def __init__(self, region_paths):
        self.region_paths = list(region_paths)
        self.start_position = Position(region_paths[0], 0)
        self.end_position= Position(region_paths[-1], 1)
        self.direction = 1
    #@property
    #def graph(self):
    #    raise Exception("NodeInterval is not connected to a graph. Cannot get graph attribute.")

    """
    @property
    def start_position(self):
        raise Exception("Start position of a NodeInterval is undefined.")
    """

    def length(self):
        raise Exception("Length of NodeInterval not supported as start/end is arbitrary with node")


class NumpyIndexedInterval(IndexedInterval):
    """
    Simple index interval which is not an interval, just an index.
    Used for fast method of finding nodes touched of a subinterval between to offsets.
    Fast writing to and from file.
    """

    def __init__(self, distance_to_node_index, node_to_distance_index, min_node, length):
        """
        Distance to node index is numpy array where element i gives node at offset i
        """
        self._length = length
        self._distance_to_node = distance_to_node_index
        self.min_node = min_node
        self._node_to_distance = node_to_distance_index
        self._nodes_in_interval_cache = None
        self._sorted_nodes_in_interval_cache = None

    def length(self):
        return self._length

    def __str__(self):
        return "Numpy indexed Interval"

    def __repr__(self):
        return self.__str__()

    def nodes_in_interval(self):
        if self._nodes_in_interval_cache is not None:
            return self._nodes_in_interval_cache

        # Hacky method to return nodes in interval
        nodes = set(np.unique(self._distance_to_node).astype(int))  
        self._nodes_in_interval_cache = nodes
        return nodes

    def get_sorted_nodes_in_interval(self):
        if self._sorted_nodes_in_interval_cache is not None:
            return self._sorted_nodes_in_interval_cache

        nodes = self._distance_to_node[np.nonzero(np.ediff1d(self._distance_to_node, to_begin=[-1]))]
        self._sorted_nodes_in_interval_cache = nodes
        return nodes


    @classmethod
    def from_interval(cls, interval, assert_on_sorted_node_ids=True):
        assert interval.start_position.offset == 0, "Currently only simple implementation supporting start offset=0"
        logging.info("Creating indexes for indexed interval")
        assert interval.graph is not None, "graph attribute cannot be None when indexing interval"

        assert isinstance(interval.graph.blocks, BlockArray), "Graph must have numpy backend"
        rps = np.array(interval.region_paths)

        if assert_on_sorted_node_ids:
            #assert np.all(np.diff(interval.region_paths) > 0), "Region paths in linear path is not sorted asc"
            if not np.all(np.diff(interval.region_paths) > 0):
                logging.warning("Region paths in linear path is not sorted asc. Linear path can be buggy.")

        min_node = rps[0]
        max_node = rps[-1]
        max_alternative = np.max(rps)
        min_alternative = np.min(rps)
        if max_alternative > max_node:
            logging.warning("Max node in path is not the last "
                            "region path (last is %d and max is %d). "
                            "Path is not sorted." % (max_alternative, max_node))
            max_node = max_alternative

        if min_alternative < min_node:
            logging.warning("Min node is not the first node. "
                            "First node is %d and min is %d" % (min_node, min_alternative))
            min_node = min_alternative

        node_sizes = interval.graph.blocks._array[rps - interval.graph.blocks.node_id_offset]
        node_to_distance = np.zeros(max_node - min_node + 2)
        length = interval.length()
        index_positions = rps - min_node
        node_to_distance[index_positions[1:]] = np.cumsum(node_sizes)[:-1]

        distance_to_node = np.zeros(length+1)
        index_positions = np.cumsum(node_sizes)[:-1]
        logging.info("Interval length: %d" % interval.length())
        logging.info("Interval end offset: %d" % interval.end_position.offset)
        logging.info("Interval end nnode size: %d" % interval.graph.blocks[interval.end_position.region_path_id].length())
        logging.info("Distance to last node end: %d" % index_positions[-1])
        assert index_positions[-1] < interval.length(), \
            "Interval is length %d, but distance to last node end is %d" % \
            (interval.length(), index_positions[-1])

        distance_to_node[index_positions] = np.diff(rps)
        distance_to_node[0] = rps[0]
        distance_to_node = np.cumsum(distance_to_node, dtype=np.uint32)

        return cls(distance_to_node, node_to_distance, min_node, length)

    def to_file(self, file_name):
        file = open(file_name, "wb")
        np.savez_compressed(file,
                 length=self._length,
                 min_node=self.min_node,
                 node_to_distance=self._node_to_distance,
                 distance_to_node=self._distance_to_node)
        file.close()

    @classmethod
    def from_file(cls, file_name):
        file = open(file_name, "rb")
        data = np.load(file)
        interval = cls(data["distance_to_node"], data["node_to_distance"],
                       data["min_node"], data["length"])
        file.close()
        return interval

    def get_exact_subinterval(self, start, end):
        assert start < end, "Start must be < than end (%d, %d)" % (start, end)
        rps = self.get_nodes_between_offset(start, end)
        start_offset = self.get_node_offset_at_offset(start, rps[0])
        end_offset = self.get_node_offset_at_offset(end-1, rps[-1])
        return DirectedInterval(
            int(start_offset), int(end_offset)+1, list(rps))

    def get_node_at_offset(self, offset):
        return self._distance_to_node[offset]

    def get_node_offset_at_offset(self, offset, true_node=None):
        node = self._distance_to_node[offset]
        if true_node is not None:
            assert node == true_node, (node, true_node)
        return offset - self.get_offset_at_node(node)

    def get_nodes_between_offset(self, start, end):
        assert end <= len(self._distance_to_node)
        assert start >= 0
        values = self._distance_to_node[start:end]
        unique = values[:-1][values[1:]!=values[:-1]]
        return np.concatenate((unique, values[-1:]))

    def position_at_offset(self, offset):
        node_offset = self.get_node_offset_at_offset(offset)
        node = self._distance_to_node[offset]
        return Position(node, node_offset)

    def get_subinterval(self, start, end):
        return NodeInterval(self.get_nodes_between_offset(start, end))

    def get_offset_at_node(self, node):
        return self._node_to_distance[node - self.min_node]

    def get_offset_at_position(self, position, direction="+"):
        if direction != "+":
            raise NotImplementedError("Not implemented")

        return self.get_offset_at_node(position.region_path_id) + position.offset
