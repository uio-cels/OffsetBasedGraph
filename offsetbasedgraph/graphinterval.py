
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

    @staticmethod
    def from_linear_interval(graph, interval):
        """ Asssuming graph block indexed on chromosome
        and sorted on start coordinate
        """
        chrom = interval.chromosome
        region_paths = graph.index[chrom]
        region_paths = [rp for rp in region_paths if
                        rp.linear_references[chrom].intersects(interval)]

        start = interval.start-region_paths[0].linear_references[chrom].start
        end = interval.end-region_paths[-1].linear_references[chrom].start

        return GraphInterval(region_paths, start, end)

    @classmethod
    def linear_coordinate_to_graph(cls, graph, chr_id, coordinate):
        """ Maps the linear coordinate to a coordinate in the block graph.
        graph_block_index should be sent as a parameter, and can be obtained
        by calling create_block_index()
        """
        for key, block in graph.blocks.items():
            if block.contains_position((chr_id, coordinate)):
                return (key, coordinate-block.linear_references[chr_id].start)
        else:
            raise Exception("Did not find coordinate (%s, %s)"
                            % (chr_id, coordinate))
        # Our coordinate can be in any of these blocks
        potential_blocks = graph.block_index[chr_id]

        for potential_block in potential_blocks:
            # Get start and end position of this block in linear genome,
            # check whether this is the correct block
            # assume always only one linear refernce, i.e. only one species
            linear_references = potential_block.linear_references.values()
            for lr in linear_references:
                start = lr.start
                end = lr.end
                if start <= coordinate and end >= coordinate:
                    return (potential_block.id, coordinate - start)

        raise Exception("No block found for chr_id %s, coordinate %d"
                        % (chr_id, coordinate))

    @classmethod
    def linear_segment_to_graph(cls, graph,  chr_id, start, end):
        """
        Takes a linear segment on hg38 and returns a ChainSegment object
        """
        block_list = []
        start_pos = 0
        end_pos = 0

        start_block, start_pos = cls.linear_coordinate_to_graph(
            graph, chr_id, start)
        end_block, end_pos = cls.linear_coordinate_to_graph(
            graph, chr_id, end)

        block_list.append(start_block)
        current_block = start_block
        while True:
            if current_block == end_block:
                break
            prev_block = current_block
            if current_block not in graph.block_edges:
                raise Exception("Did not find end block, current blocks = %s" %
                                block_list)

            edges = [e for e in graph.block_edges[current_block]
                     if chr_id in graph.blocks[e].linear_references]
            min_start = 3000000000
            min_edge = None
            for edge in edges:
                start = graph.blocks[edge].linear_references[chr_id].start
                if start < min_start:
                    min_start = start
                    min_edge = edge

            current_block = min_edge

            if current_block == prev_block:
                raise Exception("Error while traversing block. Did not find " +
                                "next block for block %s, %s" %
                                (current_block, edges))

            block_list.append(current_block)

        return cls(block_list, start_pos, end_pos)
