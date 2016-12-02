

class GraphInterval(object):
    """Used to represent intervals on a graph"""

    def __init__(self, block_list, start, end):
        self.block_list = block_list
        self.start = start
        self.end = end

    def __str__(self):
        return "%s: %s :%s" % (
            self.start,
            " ".join([block for block in self.block_list]),
            self.end)

    def __eq__(self, other):
        """ Check if start and end position as well
        the block list are equal"""
        if not self.start == other.start:
            return False
        if not self.end == other.end:
            return False
        return self.block_list == other.block_list

    @staticmethod
    def find_flanks(graph, inteval, region_paths_all):
        pass

    @staticmethod
    def from_linear_interval(graph, interval):
        """ Asssuming graph block indexed on chromosome
        and sorted on start coordinate
        """
        chrom = interval.chromosome
        region_paths_all = graph.index[chrom]
        if not region_paths_all:
            print("!!!Did not find chromosome: %s" % chrom)
            return None
        region_paths = [rp for rp in region_paths_all if
                        rp.linear_references[chrom].intersects(interval)]

        if not region_paths:
            return graph.find_interval_nice(interval)

        assert region_paths, "Region path not found (%s, %s)" % (interval, region_paths_all)

        start = interval.start-region_paths[0].linear_references[chrom].start
        end = interval.end-region_paths[-1].linear_references[chrom].start

        if start < 0 or end >= region_paths[-1].linear_references[chrom].length():
            return graph.find_interval_nice(interval)

        return GraphInterval(region_paths, start, end)
