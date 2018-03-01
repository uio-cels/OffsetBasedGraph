import json
import hashlib
import gzip
import numpy as np
import logging

class Position(object):
    """ Represents a position  in the graph

    >>> Position("chr1-1", 100)
    """

    def __init__(self, region_path_id, offset):
        """
        :param region_path_id: region path/block identifier
        :param offset: offset from start of region path (int)
        """
        self.region_path_id = region_path_id
        self.offset = offset
        #assert isinstance(offset, int), "Position offset %s is not int, %s" % (offset, type(offset))

    def __eq__(self, other):
        return self.region_path_id == other.region_path_id \
               and self.offset == other.offset

    def __str__(self):
        return "(%s:%s)" % (self.region_path_id, self.offset)

    def __repr__(self):
        return self.__str__()

    def copy(self):
        return Position(self.region_path_id, self.offset)

    def notation(self):
        """
        :return: Returns a numan readable string notation
        """
        return "%s:%d" % (self.region_path_id, self.offset)


class BaseInterval(object):
    def set_length_cache(self, l):
        """Set length cache used by lenght()
        :param l: length (int)
        """
        assert l > 0
        self.length_cache = l

    def _calculate_length(self):
        """Calculate length of interval

        :returns: lenght of interval
        :rtype: int

        """

        if len(self.region_paths) == 1:
            return self.end_position.offset - self.start_position.offset

        if not self.region_paths:
            return 0

        assert self.graph is not None,\
            "Graph is none and length cache is None and\
            interval has more than 1 rp. Cannot compute length"

        rp_lengths = [self.graph.node_size(rp)
                      for rp in self.region_paths[:-1]]
        r_sum = sum(rp_lengths)
        length = r_sum-self.start_position.offset+self.end_position.offset

        assert length >= 0,\
            "Length is %d for interval %s. r_lengths: %s."\
            % (length, self, rp_lengths)

        return length

    def __eq__(self, other):
        """Check that start and end positions and 
        region paths are the same

        :param other: Interval
        :rtype: bool

        """
        eq = self.start_position == other.start_position
        if not eq:
            return False
        eq &= self.end_position == other.end_position
        if not eq:
            return False
        eq &= self.region_paths == other.region_paths
        return eq

    def __repr__(self):
        return self.__str__()

    def length(self):
        """
        :returns: The length of the interval
        :rtype: int

        """

        if hasattr(self, "length_cache") and self.length_cache is not None:
            return self.length_cache
        length = self._calculate_length()
        self.length_cache = length
        return length


class Interval(BaseInterval):
    """
    Represents an interval on a graph

    >>> Interval(Position("chr1-0", 10), Position("chr1-1", 100),
    ["chr1-0", "alt, "chr1-1"]), graph)
    """

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None, direction=1):
        """Initialize interval with either:
        (Position, Position, [list(str)])
        (int, int, list(str)) or

        :param start_position: Start position/offset
        :param end_position: End position/offset
        :param region_paths: list of region_paths
        :param graph: The graph the interval is defined on
        """
        assert region_paths is None or isinstance(region_paths, list),\
            "Region paths must be None or list"

        if np.issubdtype(type(start_position), np.integer):
            start_position = int(start_position)
            end_position = int(end_position)

        if isinstance(start_position, int): # or np.issubdtype(start_position, np.integer):
            self._offset_init(start_position, end_position, region_paths)
        elif isinstance(start_position, Position):
            self._position_init(start_position, end_position, region_paths)
        else:
            raise Exception("Start/end position must be int or of type Position. Now: %s" % type(start_position))


        self.graph = graph
        self.rp_lens_tmp = None

        assert self.end_position.region_path_id == self.region_paths[-1]
        assert self.start_position.region_path_id == self.region_paths[0]

        self.length_cache = None
        self.direction = direction
        self._hash_cashe = None

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

    def contains(self, other, tolerance=0):
        """Check if transcription region and all exons of other
        are contained in self. Accecpt difference in start and end
        coordinates lower than tolerance

        :param other: other gene
        :param tolerance: number of basepair tolerance
        :returns: wheter self conatins other
        :rtype: bool

        """
        if other.region_paths[0] == self.region_paths[0]:
            if other.start_position.offset < self.start_position.offset-tolerance:
                return False

        if other.region_paths[-1] == self.region_paths[-1]:
            if other.end_position.offset > self.end_position.offset + tolerance:
                return False

        if other.region_paths[-1] == self.region_paths[-1] and other.region_paths[0] == self.region_paths[0] \
                and (self.start_position.offset > other.end_position.offset or self.end_position.offset < other.start_position.offset):
            return False

        if not all(rp in self.region_paths for rp in other.region_paths):
            return False

        return True

    def contains_in_order_any_direction(self, other):
        if self.contains_in_correct_order(other):
            return True

        if self.contains_in_correct_order(other.get_reverse()):
            return True

        return False

    def contains_in_correct_order(self, other):
        if not self.contains(other):
            return False

        # Check order
        rps = self.region_paths
        other_rps = other.region_paths
        n = len(other_rps)

        return any((other_rps == rps[i:i+n]) for i in range(len(rps)-n+1))

    def intersects(self, other):
        common_rps = [rp for rp in self.region_paths
                      if rp in other.region_paths]
        if not common_rps:
            return False
        for rp in common_rps:
            if rp == self.region_paths[0] and rp == other.region_paths[-1]:
                if self.start_position.offset >= other.end_position.offset:
                    continue
            if rp == self.region_paths[-1] and rp == other.region_paths[0]:
                if other.start_position.offset >= self.end_position.offset:
                    continue
            return True
        return False

    def get_position_from_offset(self, offset, rp_lens=None):
        """Get position of with offset counted from the start of
        the interval

        :param offset:
        :returns: The position in the graph
        :rtype: Position

        """
        total_offset = offset + self.start_position.offset
        if rp_lens is None:
            rp_lens = [self.graph.blocks[rp].length()
                       for rp in self.region_paths]

        for i, rp_length in enumerate(rp_lens):
            region_path = self.region_paths[i]
            if rp_length > total_offset:
                return Position(region_path, total_offset)
            total_offset -= rp_length

        if total_offset == 0:
            return Position(region_path, rp_lens[-1])

        assert False, "No offset %d from interval %s using rp lens %s" % \
            (offset, self, list(rp_lens))

    def get_adj_list(self):
        """
        :return: Returns every adjency in the interval as a list
         of tuples (from_block_id, to_block_id)
        """
        prev = None
        adjs = []
        for rp in self.region_paths:
            if prev is not None:
                adjs.append((prev, rp))
            prev = rp
        return adjs

    def copy(self):
        c = Interval(self.start_position.copy(), self.end_position.copy(),
                     list(self.region_paths), self.graph,
                     direction=self.direction)
        if self.length_cache is not None:
            c.set_length_cache(self.length_cache)
        return c

    def notation(self):
        """
        :return: Returns a human readable notation of the interval
        """

        out = ""
        out += self.start_position.notation()
        for rp in self.region_paths:
            if rp != self.start_position.region_path_id and \
                rp != self.end_position.region_path_id:
                out += ", %s" % rp
        out += ", %s" % self.end_position.notation()
        return out

    def __str__(self):
        graph = "Graph"
        if self.graph is None:
            graph = "None graph"
        return "Intv(%s, %s, %s, %s, lc=%d)" % (
            self.start_position,
            self.end_position, self.region_paths,
            graph, self.length_cache if self.length_cache is not None else 0)

    def _offset_init(self, start_offset, end_offset, region_paths):
        assert region_paths,\
            "No region paths given and start_position is int %s, %s, %s" % (
                start_offset, end_offset, region_paths)
        self.start_position = Position(region_paths[0], start_offset)
        self.end_position = Position(region_paths[-1], end_offset)
        self.region_paths = region_paths

    def _position_init(self, start_position, end_position, region_paths):
        self.start_position = start_position
        self.end_position = end_position

        if region_paths is None:
            region_paths = [start_position.region_path_id]
            if end_position.region_path_id != start_position.region_path_id:
                region_paths.append(end_position.region_path_id)
        self.region_paths = region_paths

    def to_file_line(self):
        assert isinstance(self.region_paths, list)
        #assert isinstance(self.region_paths[0], int)
        assert isinstance(self.direction, int)

        object = {"start": int(self.start_position.offset),
                  "end": int(self.end_position.offset),
                  "region_paths": [int(rp) for rp in self.region_paths],
                  "direction": self.direction
                 }
        try:
            return json.dumps(object)
        except TypeError:
            print("Json error for interval %s" % self)
            raise

    def get_reverse(self):
        assert self.graph is not None, "Graph cannot be None when reversing interval"
        start_node_size = self.graph.node_size(self.start_position.region_path_id)
        end_node_size = self.graph.node_size(self.end_position.region_path_id)
        start_offset = end_node_size - self.end_position.offset
        end_offset = start_node_size - self.start_position.offset
        reversed_rps = list([-int(rp) for rp in self.region_paths[::-1]])
        return Interval(start_offset, end_offset, reversed_rps, self.graph)

    @classmethod
    def from_file_line(cls, line, graph=None):
        object = json.loads(line)
        return cls(object["start"], object["end"], object["region_paths"], direction=object["direction"], graph=graph)

    def __deepcopy__(self, memo):
        return Interval(self.start_position, self.end_position,
                        self.region_paths, self.graph)

    def hash(self, ignore_direction=False, ignore_end_pos=False):
        """
        :return: Returns a unique hash as int of the path of the interval
        """

        if self._hash_cashe is not None:
            return self._hash_cashe

        if ignore_end_pos:
            string = str(self.start_position)
        elif ignore_direction:
            string = "%s-%s-%s"  % (str(self.start_position),
                                    str(self.end_position),
                                    str(self.region_paths))
        else:
            string = "%s-%s-%s-%d"  % (str(self.start_position),
                                str(self.end_position),
                                str(self.region_paths),
                                self.direction)
        hex_hash = hashlib.md5(string.encode('utf-8')).hexdigest()[0:15]
        h = int(hex_hash, 16)
        self._hash_cashe = h
        return h

    def position_at_offset(self, offset):
        assert offset >= 0, "Offset %d is negative" % offset
        assert offset <= self.length(), "Offset %d > interval length %d" % (offset, self.length())
        assert self.graph is not None

        current_offset = 0
        i = 0
        for rp in self.region_paths:
            node_size = self.graph.node_size(rp)
            length_to_end = node_size
            if i == 0:
                length_to_end -= self.start_position.offset
            i += 1

            if current_offset + length_to_end > offset:
                pos = Position(rp, (node_size - length_to_end - current_offset) + offset)
                return pos

            current_offset += length_to_end

        assert False, "Did not find offset %d for interval %s" % (offset, self)

    def _nodes_between_offets(self, start, end):
        nodes = []
        current_offset = 0
        i = 0
        for rp in self.region_paths:
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

    def get_subinterval(self, start_offset, end_offset):
        start_pos = self.position_at_offset(start_offset)
        end_pos = self.position_at_offset(end_offset - 1)
        #print("Start pos: %s" % start_pos)
        #print("End pos: %s" % end_pos)

        # Hack when end is on end of RP
        #if end_pos.offset == 0:
        #    end_pos = self.position_at_offset(end_offset - 1)
        #    end_pos.offset += 1
        end_pos.offset += 1

        nodes = self._nodes_between_offets(start_offset, end_offset)
        #print("Nodes: %s" % nodes)
        return Interval(start_pos, end_pos, nodes, self.graph)

    def to_areas(self):
        areas = {}
        for rp in self.region_paths:
            areas[rp] = np.zeros(self.graph.node_size(rp))

        if len(self.region_paths) == 1:
            areas[self.start_position.region_path_id][self.start_position.offset:self.end_position.offset] = 1
        else:
            for rp in self.region_paths:
                rp_size = self.graph.node_size(rp)
                if rp == self.region_paths[0]:
                    areas[rp][self.start_position.offset:rp_size] = 1

                if rp == self.region_paths[-1]:
                    areas[rp][0:self.end_position.offset] = 1

                if rp in self.region_paths[1:-1]:
                    areas[rp][0:rp_size] = 1

        return areas

    def overlap(self, other):
        assert self.graph is not None
        assert other.graph is not None
        areas = self.to_areas()
        other_areas = other.to_areas()
        #print(areas)
        #print(other_areas)
        overlap = 0
        for rp in areas:
            if rp in other_areas:
                overlap += np.sum(areas[rp][np.where(areas[rp] == other_areas[rp])])

        return overlap

    def overlaps(self, other, minimum_overlap=1):
        overlap = self.overlap(other)
        if overlap >= minimum_overlap:
            return True
        return False

    def approx_contains(self, other, allowed_mismatches=1):
        overlap = other.overlap(self)
        if  overlap > 0 and overlap >= other.length() - allowed_mismatches:
            return True
        return False

    def is_approx_equal(self, other, allowed_mismatches=1):
        overlap = self.overlap(other)
        if  overlap > 0 and overlap >= self.length() - allowed_mismatches:
            return True
        return False

    def to_indexed_interval(self, only_create_distance_index=False):
        from .indexedinterval import IndexedInterval
        return IndexedInterval(self.start_position,
                             self.end_position,
                             self.region_paths,
                             graph=self.graph,
                             only_create_distance_to_nodes_index=only_create_distance_index)

    def to_numpy_indexed_interval(self, only_create_distance_index=False):
        from .indexedinterval import NumpyIndexedInterval
        return NumpyIndexedInterval.from_interval(self)


class IntervalCollection(object):
    interval_class = Interval

    def __init__(self, intervals):
        self.intervals = intervals

    @classmethod
    def create_generator_from_file(cls, file_name, graph=None):
        f = open(file_name)
        intervals = (cls.interval_class.from_file_line(line, graph=graph) for line in f.readlines())
        f.close()
        return cls(intervals)

    @classmethod
    def create_list_from_file(cls, file_name, graph=None):
        f = open(file_name)
        intervals = [cls.interval_class.from_file_line(line, graph=graph) for line in f.readlines()]
        f.close()
        return cls(intervals)

    def __iter__(self):
        return self.intervals.__iter__()

    def to_text_file(self, file_name):
        logging.info("Writing to text file %s" % file_name)
        f = open(file_name, "w")
        i = 0
        for interval in self.intervals:
            f.writelines(["%s\n" % interval.to_file_line()])
            i += 1
        logging.info("Wrote %d lines" % i)
        f.close()
        return file_name

    def copy(self):
        return IntervalCollection.from_file(self.to_file("copy.tmp"))

    def to_file(self, file_name, text_file=False):
        if text_file:
            return self.to_text_file(file_name)
        else:
            return self.to_gzip(file_name)

    @classmethod
    def from_file(cls, file_name, text_file=False, graph=None):
        if text_file:
            return cls.create_generator_from_file(file_name, graph=graph)
        else:
            return cls.from_gzip(file_name, graph=graph)

    def to_gzip(self, file_name):
        f = gzip.open(file_name, "wb")
        try:
            for interval in self.intervals:
                line = "%s\n" % interval.to_file_line()
                f.write(line.encode())
        finally:
            f.close()

        return file_name

    @classmethod
    def from_gzip(cls, file_name, graph=None):
        import gzip
        import io
        gz = gzip.open(file_name, 'r')
        f = io.BufferedReader(gz)
        intervals = (cls.interval_class.from_file_line(line.decode("utf-8"), graph=graph)
                     for line in f.readlines())
        f.close()
        return cls(intervals)
