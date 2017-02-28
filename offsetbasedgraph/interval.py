
class Position(object):
    """
    Represents a position  in the graph
    """
    def __init__(self, region_path_id, offset):
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


class BaseInterval(object):

    def set_length_cache(self, l):
        assert l > 0
        self.length_cache = l

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

    def __eq__(self, other):
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

    def offset_init(self, start_offset, end_offset, region_paths):
        assert region_paths,\
            "No region paths given and start_position is int %s, %s, %s" % (
                start_offset, end_offset, region_paths)
        self.start_position = Position(region_paths[0], start_offset)
        self.end_position = Position(region_paths[-1], end_offset)
        self.region_paths = region_paths

    def position_init(self, start_position, end_position, region_paths):
        self.start_position = start_position
        self.end_position = end_position

        if region_paths is None:
            region_paths = [start_position.region_path_id]
            if end_position.region_path_id != start_position.region_path_id:
                region_paths.append(end_position.region_path_id)
        self.region_paths = region_paths

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None):
        assert region_paths is None or isinstance(region_paths, list),\
            "Region paths must be None or list"

        if isinstance(start_position, int):
            self.offset_init(start_position, end_position, region_paths)
        else:
            self.position_init(start_position, end_position, region_paths)

        self.graph = graph
        self.rp_lens_tmp = None

        assert self.end_position.region_path_id == self.region_paths[-1]
        assert self.start_position.region_path_id == self.region_paths[0]

        self.length_cache = None

    def contains_rp(self, rp):
        return (rp in self.region_paths)

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

    def __deepcopy__(self, memo):
        return Interval(self.start_position, self.end_position,
                        self.region_paths, self.graph)

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
        return "%s, %s, [%s]" % (
            self.start_position.offset,
            self.end_position.offset,
            ', '.join([str(s) for s in self.region_paths]))

    def __str__(self):
        graph = "Graph"
        if self.graph is None:
            graph = "None graph"
        return "Intv(%s, %s, %s, %s, lc=%d)" % (
            self.start_position,
            self.end_position, self.region_paths,
            graph, self.length_cache if self.length_cache is not None else 0)

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

        assert False, "No offset %d from interval %s using rp lens %s" % (offset, self, list(rp_lens))

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
