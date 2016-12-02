from __future__ import absolute_import
from __future__ import print_function
from collections import defaultdict
from .graphinterval import GraphInterval
from .regionpath import RegionPath
from .linearinterval import LinearInterval

from .config import DEBUG


class OffsetBasedGraph():
    """
    Class representing an offset based graph,
    by storing  a list of region paths (sometimes referred to as blocks)
    and a list of edges connecting these region paths.

    :Example:

    >>> graph = OffsetBasedGraph()
    >>> graph.create_graph()
    >>> print(graph.blocks)
    >>> print(graph.block_edges)

    """

    def __init__(self, name):
        """
        Initialises the offset based graph with a name. The graph will be empty
        at this point.
        :param name: A textual name of the graph
        :return:
        """
        self.blocks = {}
        self.block_edges = {}
        self.block_edges_back = {}
        self.cur_id = 0
        self.cur_id_counter = {}  # Number of times a base id of block
                                  # has been used

        self.block_index = {}
        self.name = name

    def __str__(self):
        elements = [str(block.id) + ": " + str(block)
                    for block in self.blocks.values()]
        for key, val in self.block_edges.items():
            elements.append("%s: %s" % (str(key), str(val)))
        return "\n".join(elements)

    def deep_copy(self, name):
        """ Copies the graph.

        :param name: Name of the copy
        :return: Returns the copied graph.
        """
        new_graph = OffsetBasedGraph(name)
        new_graph.blocks = self.blocks.copy()
        new_graph.block_edges = self.block_edges.copy()
        new_graph.block_edges_back = self.block_edges_back.copy()
        new_graph.cur_id = self.cur_id
        new_graph.cur_id_counter = self.cur_id_counter.copy()
        new_graph._create_block_index()
        return new_graph

    def find_overlapping_blocks(self, block):
        """
        Finds region paths in the graph that are overlapping with
        the given block (i.e. comming
        from the same source sequences). Uses linear references to determine
        this.

        :param block:  ID of the region path
        :return: Returns a list of region paths (blocks) that are overlapping
        """
        main_path_interval = self.blocks[block].main_path_linear_reference
        return [name for name in self.blocks.keys() if
                self.blocks[name].main_path_linear_reference.intersects(main_path_interval)]

    def get_first_block(self, chromosome_id):
        """
        Get the first region path in a chromosome.
        WARNING: Requires specific naming scheme of region paths,
        where first block always has id chr[number]-1

        :param chromosome_id:
        :return: The first region path of given chromosome
        """
        first_block_id = chromosome_id + "-1"
        return self.blocks[first_block_id]

    def get_last_block(self, chromosome_id):
        """
        Gets the last region path in a chromosome
        WARNING: Requires specific naming scheme (see get_first_block)

        :param chromosome_id:
        :return: Returns the last region path of given chromosome
        """

        # Starts at first block, returns last block that it can
        # find by traversing
        block = self.get_first_block(chromosome_id)
        while True:
            next_blocks = self.block_edges[block.id]
            if not next_blocks:
                return block.id
            block = self.blocks[next_blocks[0]]
        return block.id

    def remove_block(self, name):
        """
        Removes a region path (block).

        :param name: ID of region path / block
        :return:
        """
        del self.blocks[name]

        if name in self.block_edges:
            edges = self.block_edges[name]
            del self.block_edges[name]
        else:
            edges = []

        if name in self.block_edges_back:
            back_edges = self.block_edges_back[name]
            del self.block_edges_back[name]
        else:
            back_edges = []

        for edge in edges:
            self.block_edges_back[edge].remove(name)

        for back_edge in back_edges:
            self.block_edges[back_edge].remove(name)

    def set_back_edges(self):
        """
        Updates the graph with correct back edges from each block.
        This method should always be called after creating the graph
        """
        self.block_edges_back = {}
        for block in self.blocks.values():
            for other in self.blocks.values():

                if other.id not in self.block_edges:
                    continue
                if block.id not in self.block_edges[other.id]:
                    continue
                if block.id not in self.block_edges_back:
                    self.block_edges_back[block.id] = []
                self.block_edges_back[block.id].append(other.id)

    def get_main_path_linear_references(self):
        """
        This method works for a graph created from GRCh38.
        For every region path, find the "parallell" position on the main path,
        i.e. the linear reference on the main path corresponding to the
        alternative loci. Stores this in block.main_path_linear_reference
        """
        # TODO: do this with info in graph
        raise NotImplementedError()

        for block in self.blocks.values():
            if "alt" in block.id:
                # This is an alt, block: Do a db call to find position on main path
                alt_id = list(block.linear_references.values())[0].chromosome
                res = self.db.fetch_all("SELECT chrom, chromStart, chromEnd FROM altLocations where name='%s'" % (alt_id))
                res = res[0]

                chrom = res["chrom"]
                start = int(res["chromStart"])
                end = int(res["chromEnd"])
                linref = LinearInterval("hg38", chrom, start, end)
            else:
                linref = list(block.linear_references.values())[0]

            block.main_path_linear_reference = linref

    def _create_block_index(self):
        """ Creates a dictionary where indices are linear block ids (chromosome ids
        and alt loci ids), and values are list of blocks in the graph that come
        from these blocks.
        This dict is used when mapping linear intervals to graph intervals
        """
        self.index = defaultdict(list)
        for block_id, block in self.blocks.items():
            for linear_reference in block.linear_references.values():
                self.index[linear_reference.chromosome].append(block)

        for k, v in self.index.items():
            v.sort(key=lambda block: block.linear_references[k].start)

    def get_blocks(self, lin_ref):
        """
        Finds all region paths that contains a linear references

        :param lin_ref: The linear reference
        :return: Returns region paths
        """
        return [block for block in
                list(self.blocks.values()) if block.contains(lin_ref)]

    def get_intersecting_blocks(self, lin_ref):
        """
        Finds all region paths that intersects with a linear references

        :param lin_ref: The linear reference
        :return: Returns region paths
        """
        return [block for block in
                list(self.blocks.values()) if block.intersects(lin_ref)]

    def get_block(self, lin_ref):
        filtered = self.get_blocks(lin_ref)
        if not filtered:
            return None
        return filtered[0]

    def set_parallel(self, lin_refs):
        graf_refs = [GraphInterval.from_linear_interval(lin_ref)
                     for lin_ref in lin_refs]

    def divide_block(self, block, offset):
        blocks = block.split([offset])
        self.block_edges[blocks[1].id] = self.block_edges[block.id]
        self.block_edges[block.id] = [block.id]
        self.blocks[blocks[0].id] = blocks[0]
        self.blocks[blocks[1].id] = blocks[1]

    def set_equal(self, lin_refs):
        """ """
        graph_refs = [GraphInterval.from_linear_interval(lin_ref)
                      for lin_ref in lin_refs]

        start_blocks = [graph_ref.block_list[0] for graph_ref in graph_refs]
        end_blocks = [graph_ref.block_list[-1] for graph_ref in graph_refs]

    def find_coordinate_nice(self, chrom, offset):
        """Find the coordinate and allow overflowing
        region paths"""
        region_paths = self.index[chrom]
        assert region_paths, "Did not find region path %s" % chrom
        containing_rps = [rp for rp in region_paths
                          if rp.contains_position(chrom, offset)]
        if containing_rps:
            return containing_rps[0].id, offset-containing_rps[0].linear_references[chrom].start

        ds = [rp.distance_to_position(chrom, offset) for rp in region_paths]
        first_dist = min(ds)
        closest_rp = min(zip(region_paths, ds), key=lambda x: x[1])[0]
        lin_ref = closest_rp.linear_references[chrom]

        traversal_dict = self.block_edges
        offset_calc = lambda block, overflow: overflow

        if offset < lin_ref.start:
            offset_calc = lambda block, overflow: block.length()- overflow
            traversal_dict = self.block_edges_back

        distance_travelled = 0
        cur_block_id = closest_rp.id
        while distance_travelled < first_dist:
            next_ids = traversal_dict[cur_block_id]
            cur_block_id = next_ids[0]
            distance_travelled += self.blocks[cur_block_id].length()

        overflow = distance_travelled - first_dist
        new_offset = offset_calc(self.blocks[cur_block_id], overflow)
        return cur_block_id, new_offset

    def find_interval_nice(self, interval):
        start_block_id, start_offset = self.find_coordinate_nice(
            interval.chromosome,
            interval.start)
        end_block_id, end_offset = self.find_coordinate_nice(
            interval.chromosome,
            interval.end)

        cur_block = start_block_id
        blocks = [cur_block]
        while cur_block != end_block_id:
            next_blocks = self.block_edges[cur_block]
            if interval.chromosome in next_blocks:
                cur_block = interval.chromosome
            else:
                cur_block = next_blocks[0]
            blocks.append(cur_block)

        return GraphInterval(blocks, start_offset, end_offset)


    def merge_lin_refs(self, org_lin_ref, new_lin_ref):
        """
        Merges two parts of the graph into one region path. The two parts to
        merge are the two linear references

        :param org_lin_ref: Linear reference to merge
        :param new_lin_ref: Linear reference to merge
        :return:
        """
        org_block = self.get_block(org_lin_ref)
        if org_block is None:
            return

        pre_ref, post_ref = self.split_block(org_block, org_lin_ref)
        in_edges = self.get_previous_blocks(org_block.id)
        out_edges = []
        if org_block.id in self.block_edges:
            out_edges = self.block_edges[org_block.id]

        lin_refs = [pre_ref, post_ref, new_lin_ref, org_lin_ref]
        blocks = [self._add_block({lr.chromosome: lr}) for lr in lin_refs]

        pre_block = blocks[0]
        post_block = blocks[1]
        new_block = blocks[2]
        org_block_1 = blocks[3]
        if pre_block.is_empty():
            for p_block in in_edges:
                self.add_block_edge(p_block, org_block_1.id)
                self.add_block_edge(p_block, new_block.id)
            del self.blocks[pre_block.id]
        else:
            for p_block in in_edges:
                self.add_block_edge(p_block, pre_block.id)
            self.add_block_edge(pre_block.id,  new_block.id)
            self.add_block_edge(pre_block.id,  org_block_1.id)

        if post_block.is_empty():
            self.block_edges[org_block_1.id] = out_edges[:]
            self.block_edges[new_block.id] = out_edges[:]
            del self.blocks[post_block.id]
        else:
            self.add_block_edge(new_block.id, post_block.id)
            self.add_block_edge(org_block_1.id, post_block.id)
            self.block_edges[post_block.id] = out_edges[:]

        # Remove org_block
        if org_block.id in self.block_edges:
            del self.block_edges[org_block.id]
        del self.blocks[org_block.id]

        for p_block in in_edges:
            self.block_edges[p_block].remove(org_block.id)

    def split_block(self, block, lin_ref):
        """
        Splits a block at the given linear reference. Returns two linear
        references, one before and one after lin_ref
        """
        current_lin_refs = [lr for lr in list(block.linear_references.values())
                            if lr.genome_id == lin_ref.genome_id]
        assert len(current_lin_refs) == 1
        current_lin_ref = current_lin_refs[0]
        splitted = []
        splitted.append(
            LinearInterval(lin_ref.genome_id,
                           lin_ref.chromosome,
                           current_lin_ref.start,
                           lin_ref.start,
                           lin_ref.strand)
        )
        splitted.append(
            LinearInterval(
                lin_ref.genome_id,
                lin_ref.chromosome,
                lin_ref.end,
                current_lin_ref.end,
                lin_ref.strand)
        )
        return splitted

    def get_next_blocks(self, block):
        return [self.blocks[bid] for bid in self.block_edges[block.id]]

    def merge_linear_segments(self, lin_seg_1, lin_seg_2):
        """
        Merges the graph where the two linear segments are
        """
        # First find the two blocks where each linear segment is
        block1 = self.get_block(lin_seg_1)
        block2 = self.get_block(lin_seg_2)
        if (block1 is None or block2 is None):
            return
        # The new block that is a merge of two linear segments
        mid_block = self._add_block({lin_seg_1.chromosome: lin_seg_1,
                                     lin_seg_2.chromosome: lin_seg_2})

        pre1, post1 = self.split_block(block1, lin_seg_1)
        pre2, post2 = self.split_block(block2, lin_seg_2)
        pre_block1 = RegionPath(block1.id, {pre1.chromosome: pre1})
        pre_block2 = RegionPath(block2.id, {pre2.chromosome: pre2})

        self.blocks[pre_block1.id] = pre_block1
        self.blocks[pre_block2.id] = pre_block2
        post_block1 = self._add_block({post1.chromosome: post1})
        post_block2 = self._add_block({post2.chromosome: post2})
        if block1.id in self.block_edges:
            self.block_edges[post_block1.id] = self.block_edges[block1.id][:]

        if block2.id in self.block_edges:
            self.block_edges[post_block2.id] = self.block_edges[block2.id][:]

        self.block_edges[pre_block1.id] = []
        self.block_edges[pre_block2.id] = []

        self.add_block_edge(pre_block1.id, mid_block.id)
        self.add_block_edge(pre_block2.id, mid_block.id)
        self.add_block_edge(mid_block.id, post_block1.id)
        self.add_block_edge(mid_block.id, post_block2.id)

    def add_block_edge(self, block_id, next_block_id):
        """
        Adds an edge between two blocks (region paths)

        :param block_id: The block that the edge goes from
        :param bnext_block_id: The block that the edge goes to

        :Example:

        >>> graph = OffsetBasedGraph()
        >>> graph.add

        """
        if block_id not in self.blocks:
            raise Exception("Edge id (%s) not in blocks: %s"
                            % (block_id, self.blocks))

        if block_id not in self.block_edges:
            self.block_edges[block_id] = [next_block_id]
            return

        self.block_edges[block_id].append(next_block_id)

    def get_subgraph(self, linear_interval, correct_alt):
        """. warnings also:: Modifies current graph."""
        blocks = self.get_intersecting_blocks(linear_interval)
        min_C = 4000000000
        max_C = 0
        min_block = None
        max_block = None
        for block in blocks:
            lin_seg = list(filter(
                lambda li: li.chromosome == linear_interval.chromosome,
                block.linear_references.values()))[0]
            if lin_seg.start < min_C:
                min_block = block
                min_C = lin_seg.start
            if lin_seg.end > max_C:
                max_block = block
                max_C = lin_seg.start

        subgraph = OffsetBasedGraph("pruned graph")

        cur_block = min_block
        subgraph.blocks[cur_block.id] = cur_block
        subgraph.start_block = min_block.id

        list(min_block.linear_references.values())[0].start = linear_interval.start
        list(max_block.linear_references.values())[0].end = linear_interval.end
        if min_block == max_block:
            return subgraph
        while True:
            found = False
            for new_block in self.block_edges[cur_block.id]:
                if any(["alt" in lr.chromosome and not lr.chromosome == correct_alt
                        for lr in self.blocks[new_block].linear_references.values()]):
                    continue
                subgraph.blocks[new_block] = self.blocks[new_block]
                if self.blocks[new_block] == max_block:
                    found = True
                    break
            if found:
                break
            if not self.block_edges[cur_block.id]:
                break
            cur_block = self.blocks[self.block_edges[cur_block.id][0]]

        for block in subgraph.blocks:
            subgraph.block_edges[block] = list(filter(
                lambda b: b in subgraph.blocks,
                self.block_edges[block]))

        return subgraph

    def get_previous_blocks(self, block_id):
        previous_blocks = []
        for p_block, edges in self.block_edges.items():
            if block_id in edges:
                previous_blocks.append(p_block)
        return previous_blocks

    def add_chromosome(self, genome_id, chromosome, chromosome_size):
        """
        Adds a whole chromosome as a linear block
        """
        linear_references = {
            chromosome: LinearInterval(
                genome_id, chromosome, 1, chromosome_size)}
        self._add_block(linear_references)

    def _generate_id_from_linear_references(self, lrs):
        suggestion = ""
        sorted_genome_ids = sorted([lr for lr in lrs])
        if len(sorted_genome_ids) > 1:
            suggestion += "merge-"
        for genome_id in sorted_genome_ids:
            # suggestion += genome_id + "." + lrs[genome_id].chromosome + "-"
            # Skip non-alternative
            if len(sorted_genome_ids) == 1 or "alt" in lrs[genome_id].chromosome:
                suggestion += lrs[genome_id].chromosome + "-"

        base_id = suggestion

        if suggestion in self.cur_id_counter:
            suggestion += str(self.cur_id_counter[base_id] + 1)
            self.cur_id_counter[base_id] += 1
        else:
            suggestion += "0"
            self.cur_id_counter[base_id] = 0

        return suggestion

    def _add_block(self, linear_references):
        """
        Adds block for one genome
        linear_referes is a dict, having keys that are genome id and values
        that are LinearInterval objects
        """
        id = self._generate_id_from_linear_references(linear_references)
        block = RegionPath(id, linear_references)
        self.blocks[id] = block
        return block

    @classmethod
    def create_graph_from_blocks(cls, blocks, block_edges):
        graph = cls("hg38")
        graph.blocks = blocks
        graph.block_edges = block_edges
        graph._create_block_index()
        graph.set_back_edges()
        return graph

    @classmethod
    def create_graph(cls, chrom_sizes, alt_loci_infos, dummy_graph=False):
        """
        Creates and returns a block graph from two dicts
        """

        graph = cls("hg38")

        for chromosome, length in chrom_sizes.items():
            if "_" not in chromosome:
                graph.add_chromosome("hg38", chromosome, length)

        for info in alt_loci_infos:
            main_interval = LinearInterval(
                "hg38", info["chrom"], info["chromStart"], info["chromEnd"])

            alt_interval = LinearInterval(
                "hg38", info["name"], info["altStart"],
                info["altStart"]+info["length"])

            graph.merge_lin_refs(main_interval, alt_interval)

        graph._create_block_index()
        graph.set_back_edges()
        # if not dummy_graph:
        #    graph.get_main_path_linear_references()
        return graph

    def include_alignments(self, alignments):
        for lr1, lr2 in alignments:
            self.merge_linear_segments(lr1, lr2)
