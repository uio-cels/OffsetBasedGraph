import json
import hashlib
import gzip
import io
import numpy as np


class IntervalCollection(object):
    def __init__(self, intervals):
        self.intervals = intervals

    @classmethod
    def create_generator_from_file(cls, file_name):
        f = open(file_name)
        intervals = (Interval.from_file_line(line) for line in f.readlines())
        f.close()
        return cls(intervals)

    def __iter__(self):
        return self.intervals.__iter__()

    def to_text_file(self, file_name):
        f = open(file_name, "w")
        for interval in self.intervals:
            f.writelines(["%s\n" % interval.to_file_line()])
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
    def from_file(cls, file_name, text_file=False):
        if text_file:
            return cls.create_generator_from_file(file_name)
        else:
           return cls.from_gzip(file_name)

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
    def from_gzip(cls, file_name):
        import gzip
        import io
        gz = gzip.open(file_name, 'r')
        f = io.BufferedReader(gz)
        intervals = (Interval.from_file_line(line.decode("utf-8")) for line in f.readlines())
        f.close()
        return cls(intervals)


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
            "Length is %d for interval %s. r_lengths: %s. Graph: %s"\
            % (length, self, rp_lengths)

        return length

    def __eq__(self, other):
        """Check that start and end positions and 
        region paths are the same

        :param other: Interval
        :rtype: bool

        """
        eq = self.start_position == other.start_position
        eq &= self.end_position == other.end_position
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

        if isinstance(start_position, int) or isinstance(start_position, np.int64):
            self._offset_init(start_position, end_position, region_paths)
        else:
            self._position_init(start_position, end_position, region_paths)

        self.graph = graph
        self.rp_lens_tmp = None

        assert self.end_position.region_path_id == self.region_paths[-1]
        assert self.start_position.region_path_id == self.region_paths[0]

        self.length_cache = None
        self.direction = direction

    def contains_rp(self, region_path):
        """Check if interval uses region_path

        :param region_path: region path id
        :rtype: bool
        """
        return (region_path in self.region_paths)

    def contains_position(self, position):
        if position.region_path_id not in self.region_paths:
            return False
        if position.region_path_id in self.region_paths[1:-1]:
            return True
        if position.region_path_id == self.region_paths[0]:
            if position.offset < self.start_position.offset:
                return False
        if position.region_path_id == self.region_paths[-1]:
            if position.offset >= self.end_position.offset:
                return False
        return True

    def contains(self, other, tolerance=0):
        """Check if transcription region and all exons of other
        are contained in self. Accecpt difference in start and end
        coordinates lower than tolerance

        :param other: other gene
        :param tolerance: number of basepair tolerance
        :returns: wheter self conatins other
        :rtype: bool

        """
        if not all(rp in self.region_paths for rp in other.region_paths):
            return False
        if other.region_paths[0] == self.region_paths[0]:
            if other.start_position.offset < self.start_position.offset-tolerance:
                return False
        if other.region_paths[-1] == self.region_paths[-1]:
            if other.end_position.offset > self.end_position.offset + tolerance:
                return False

        return True

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
                     list(self.region_paths), self.graph)
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
        object = {"start": self.start_position.offset,
                  "end": self.end_position.offset,
                  "region_paths": self.region_paths,
                  "direction": self.direction
                 }
        return json.dumps(object)

    @classmethod
    def from_file_line(cls, line):
        object = json.loads(line)
        return cls(object["start"], object["end"], object["region_paths"], direction=object["direction"])

    def __deepcopy__(self, memo):
        return Interval(self.start_position, self.end_position,
                        self.region_paths, self.graph)

    def hash(self):
        """
        :return: Returns a unique hash as int of the path of the interval
        """
        string = "%s-%s-%s-%d"  % (str(self.start_position),
                                str(self.end_position),
                                str(self.region_paths),
                                self.direction)
        hex_hash = hashlib.md5(string.encode('utf-8')).hexdigest()[0:15]
        return int(hex_hash, 16)



