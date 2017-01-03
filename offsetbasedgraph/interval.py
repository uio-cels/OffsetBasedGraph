
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

    def length(self):
        """

        :returns: The length of the interval
        :rtype: int

        """
        try:
            r_lengths = [self.graph.blocks[rp].length()
                     for rp in self.region_paths[:-1]]
            r_sum = sum(r_lengths)
            return r_sum-self.start_position.offset+self.end_position.offset

        except KeyError as e:
            print("Error on interval %s on graph %s" % (str(self), str(self.graph)))
            print(str(e))
            raise

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

    def __deepcopy__(self, memo):
        return Interval(self.start_position, self.end_position, self.region_paths, self.graph)

    def copy(self):
        return Interval(self.start_position, self.end_position, self.region_paths, self.graph)

    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None):

        if isinstance(start_position, int):
            self.start_position = Position(region_paths[0], start_position)
        else:
            self.start_position = start_position

        if isinstance(end_position, int):
            self.end_position = Position(region_paths[-1], end_position)
        else:
            self.end_position = end_position

        # By default include start and end region path
        if region_paths is None:
            region_paths = [start_position.region_path_id]
            if end_position.region_path_id != start_position.region_path_id:
                region_paths.append(end_position.region_path_id)

        self.region_paths = region_paths
        self.graph = graph

        # Sanity check interval
        # assert self.start_position.region_path_id in self.region_paths
        # assert self.end_position.region_path_id in self.region_paths

        if self.graph is None:
            return
        # for rp in region_paths:
        # assert rp in graph.blocks, "Region path %s not in graph \n%s" % (rp, graph)

        # Check offsets
#         max_offset = graph.blocks[self.region_paths[-1]].length()
#         msg = "Offset %d in block %d with size %d. Interval: %s" % (
#             self.end_position.offset,
#             self.region_paths[-1],
#             graph.blocks[self.region_paths[-1]].length(), self.__str__())
# 
#         assert self.end_position.offset <= max_offset, msg

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
        return "Intv(%s, %s, %s, %s)" % (self.start_position,
                            self.end_position, self.region_paths,
                            graph)

    def get_position_from_offset(self, offset, rp_lens=None):
        """Get position of with offset counted from the start of
        the interval

        :param offset:
        :returns: The position in the graph
        :rtype: Position

        """

        # assert offset <= self.length(), \
        #    "Offset %d is larger than total length of interval %d in interval %s" \
        #     % (offset, self.length(), self.__str__())

        total_offset = offset + self.start_position.offset
        if rp_lens is None:
            print("Finding rp_lens")
            rp_lens = [self.graph.blocks[rp].length()
                       for rp in self.region_paths]

        #if len(rp_lens) == 1:
        #    return Position()


        for i, region_path in enumerate(self.region_paths):
            #print("Inside")
            #print(self.graph)
            #print(region_path)
            rp_length = rp_lens[i]
            #print(rp_length)
            if rp_length > total_offset:
                return Position(region_path, total_offset)
            total_offset -= rp_length


        print("Tried to get offset %d from interval %s" % (offset, self))
        print(rp_lens)

        assert False

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
