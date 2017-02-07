
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

    def __deepcopy__(self, memo):
        return Position(self.region_path_id, self.offset)


class Interval(object):

    def offset_init(self, start_offset, end_offset, region_paths):
        assert region_paths, "No region paths given and start_position is int %s, %s, %s" % (start_position, end_position, region_paths)
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

        #if self.end_position.region_path_id not in self.region_paths:
        #    self.region_paths.append(self.end_position.region_path_id)

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None):
        assert region_paths is None or isinstance(region_paths, list), "Region paths must be None or list"

        if isinstance(start_position, int):
            self.offset_init(start_position, end_position, region_paths)
        else:
            self.position_init(start_position, end_position, region_paths)

        # By default include start and end region path
        self.graph = graph

        #if self.graph is not None:
            #assert self.start_position.offset < self.graph.blocks[self.region_paths[0]].length()

        assert self.end_position.region_path_id == self.region_paths[-1]
        assert self.start_position.region_path_id == self.region_paths[0]

        self.length_cache = None

    def set_length_cache(self, l):
        assert l > 0
        self.length_cache = l

    def length(self):
        """

        :returns: The length of the interval
        :rtype: int

        """

        if self.length_cache is not None:
            return self.length_cache

        if len(self.region_paths) == 1:
            self.length_cache = self.end_position.offset - self.start_position.offset
            return self.length_cache

        if not self.region_paths:
            self.length_cache = 0
            return 0

        assert self.graph is not None, "Graph is none and length cache is None and interval has more than 1 rp. Cannot compute length"

        r_lengths = [self.graph.blocks[rp].length()
                     for rp in self.region_paths[:-1]]
        r_sum = sum(r_lengths)
        length = r_sum-self.start_position.offset+self.end_position.offset
        self.length_cache = length

        assert length >= 0, "Length is %d for interval %s. r_lengths: %s. Graph: %s" % (length, self, r_lengths, self.graph)
        return length

    def starts_at_rp(self):
        return self.start_position.offset == 0

    def ends_at_rp(self):
        rp_end = self.graph.blocks[self.end_position.region_path_id].length()
        self.end_position.offset == rp_end

    def split(self, offsets):
        """Split the interval at the given offsets

        :param offsets: list of int
        :returns: list of intervals
        :rtype: list

        """
        split_points = [self.get_position_from_offset(offset) for offset
                        in offsets]
        region_path_idxs = {rp: i for i, rp in enumerate(self.region_paths)}
        split_intervals = []

        starts = [self.start_position] + split_points
        ends = split_points + [self.end_position]

        for start_pos, end_pos in zip(starts, ends):
            start_idx = region_path_idxs[start_pos.region_path_id]
            end_idx = region_path_idxs[end_pos.region_path_id]
            region_paths = self.region_paths[start_idx:end_idx+1]
            split_intervals.append(Interval(start_pos, end_pos, region_paths))

        return split_intervals

    def join(self, other):
        """Create new interval as union of self and other
        other must start at self.end_position

        :param other: Interval
        :returns: joined interval
        :rtype: Interval

        """

        assert self.end_position == other.start_position
        region_paths = self.region_paths[:-1] + other.region_paths
        return Interval(self.start_position, other.end_position, region_paths)

    def filter_to_main(self):
        filtered_rps = [rp for rp in self.region_paths if
                        self.graph.is_main_name(rp)]
        if not filtered_rps:
            return None
        start_offset = self.start_position.offset if self.graph.is_main_name(self.region_paths[0]) else 0

        end_offset = self.end_position.offset if self.graph.is_main_name(self.region_paths[-1]) else self.graph.blocks[filtered_rps[-1]].length()
        return Interval(start_offset, end_offset, filtered_rps)

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

    def approx_equals(self, other, tolerance=0):
        if not self.region_paths == other.region_paths:
            return False
        if abs(self.start_position.offset-other.start_position.offset) > tolerance:
            return False
        if abs(self.end_position.offset-other.end_position.offset)>tolerance:
            return False
        return True

    def __deepcopy__(self, memo):
        return Interval(self.start_position, self.end_position,
                        self.region_paths, self.graph)

    def copy(self):
        c = Interval(self.start_position, self.end_position,
                        self.region_paths, self.graph)
        if self.length_cache is not None:
            c.set_length_cache(self.length_cache)
        return c

    def __eq__(self, other):
        eq = self.start_position == other.start_position
        eq &= self.end_position == other.end_position
        eq &= self.region_paths == other.region_paths
        return eq

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        graph = "Graph"
        if self.graph is None:
            graph = "None graph"
        return "Intv(%s, %s, %s, %s, lc=%d)" % (
            self.start_position,
            self.end_position, self.region_paths,
            graph, self.length_cache if self.length_cache is not None else 0)

    def diff(self, other):
        codes = []
        for rp in self.region_paths:
            if any(o_rp == rp for o_rp in other.region_paths):
                codes.append("E")
            elif any(self.graph.are_paralell(rp, o_rp) for
                     o_rp in other.region_paths):
                codes.append("P")
            else:
                codes.append("N")
        return codes

        codes = []
        if not any(rp in other.region_paths for rp in self.region_paths):
            return ["D"]
        if len(self.region_paths) != len(other.region_paths):
            codes.append("dl")
        if self.region_paths[0] != other.region_paths[0]:
            codes.append("ds")
        if self.region_paths[-1] != other.region_paths[-1]:
            codes.append("de")
        if len(self.region_paths) > 2 and len(other.region_paths) > 2:
            if self.region_paths[1] != other.region_paths[1]:
                codes.append("dm")
        return codes

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

        for i, region_path in enumerate(self.region_paths):
            rp_length = rp_lens[i]
            if rp_length > total_offset:
                return Position(region_path, total_offset)
            total_offset -= rp_length

        assert False, "No offset %d from interval %s using rp lens %s" % (offset, self, rp_lens)

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

    def partition_region_paths(self):
        rps = self.region_paths
        start = self.start_position.offset
        partitions = []
        for rp in rps:
            if rp == rps[-1]:
                end = self.end_position.offset
            else:
                end = self.graph.blocks[rp].length()
            partitions.append(Interval(start, end, [rp]))
            start = 0
        return partitions


def interval_factory(graph):
    """Create function that initializes
    intervals without specifying the graph
    every time

    :param graph: The graph for the intervals
    :returns: initializer for Interval
    :rtype: func

    """
    def get_interval(*args, **kwargs):
        """Initialize an interval with the specified graph

        :returns: 
        :rtype: Interval

        """
        kwargs["graph"] = graph
        return Interval(*args, **kwargs)

    return get_interval
