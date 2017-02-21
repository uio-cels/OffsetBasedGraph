from collections import defaultdict
from .interval import Interval, Position
from .util import takes
from .translation import Translation
import pickle
import os

class Block(object):
    def __init__(self, length):

        assert length > 0

        self._length = length

    def length(self):
        return self._length

    def to_interval(self):
        pass

    def __eq__(self, other):
        return self.length() == other.length()

    def __str__(self):

        return "B(%s)" % self._length

    def __repr__(self):
        return self.__str__()


class Graph(object):
    adj_list = defaultdict(list)
    reverse_adj_list = defaultdict(list)

    # Graph alterations
    def __init__(self, blocks, adj_list, create_reverse_adj_list=True, rev_adj_list=None):
        """
        Inits the graph with a list of blocks and an adjency list
        :param blocks:
        :param adj_list:
        :return:
        """
        self.blocks = blocks
        self.adj_list = defaultdict(list, adj_list)
        if rev_adj_list is not None:
            self.reverse_adj_list = rev_adj_list
        elif create_reverse_adj_list:
            self.reverse_adj_list = self._get_reverse_edges(adj_list)

        self._id = max([b for b in blocks if isinstance(b, int)] + [-1])

        if isinstance(self._id, str):
            self._id = 0

    def copy(self):
        """Make a copy of the graph

        :returns: copy of the graph
        :rtype: Graph

        """
        import copy
        #return copy.deepcopy(self)

        new_blocks = {}
        new_adjs = {}
        for b in self.blocks:
            new_blocks[b] = Block(self.blocks[b].length())
        # new_blocks = copy.copy(self.blocks)

        for b in self.adj_list:
            new_adjs[b] = list(self.adj_list[b])
        # new_adjs = self.adj_list.copy()


        new_graph = Graph(new_blocks, new_adjs, True)
        #new_graph.reverse_adj_list = self.reverse_adj_list.copy()
        return new_graph

    def _next_id(self):
        """Make a new id and return it

        :returns: new id
        :rtype: int

        """
        self._id += 1
        return self._id

        # Debug sanity checking
        for block in adj_list:
            assert block in self.blocks, "Edge found from block that is not in blocks"
            for block2 in adj_list[block]:
                assert block2 in self.blocks, "Edge found going to non-existing block"

    @staticmethod
    def from_file(file_name):
        """
        Reads from file and creates the graph
        :param file_name: File name
        :return:
        """
        if os.path.isfile("%s" % file_name):
            with open("%s" % file_name, "rb") as f:
                return pickle.loads(f.read())
        else:
            print("Warning: Graph not found" % file_name)
            return None

    def to_file(self, file_name):
        """
        Writes the graph to file so that it later can be
        recreated using the from_file method
        :param file_name: File name
        :return:
        """
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    def get_edges_as_list(self):
        """
        :return: Returns edges as tuples (from, to) in a list
        """
        edges = []
        for block in self.adj_list:
            for to in self.adj_list[block]:
                edges.append((block, to))

        return edges

    @staticmethod
    def _generate_translation(interval_a, interval_b,
                              sub_intervals_a, sub_intervals_b, ids):

        reverse_dict = {ids: list(pair)
                        for pair in zip(sub_intervals_a, sub_intervals_b)}

        involved_region_paths = interval_a.region_paths+interval_b.region_paths
        front_dict = {rp: [i for i in sub_intervals_a+sub_intervals_b
                           if rp in i.region_paths]
                      for rp in set(involved_region_paths)}

        return Translation(front_dict, reverse_dict)

    def _get_edges_to_subintervals(self, intervals):
        """Return a list of incoming and outgoing edges at the
        divisions between the sub_intervals

        :param intervals: list of consecutive intervals
        :returns: incoming_edges, outgoing_edges
        :rtype: (list, list)

        """
        pre_edges = []
        post_edges = []
        for i_pre, i_post in zip(intervals[:-1], intervals[1:]):
            if not i_post.starts_at_rp():
                post_edges.append([])
                pre_edges.append([])
                continue
            pre_rps = self.reverse_adj_list[i_post.region_paths[0]]
            pre_edges.append(
                [rp for rp in pre_rps if rp != i_pre.region_paths[-1]]
                )

            post_rps = self.adj_list[i_pre.region_paths[-1]]
            post_edges.append(
                [rp for rp in post_rps if rp != i_post.region_paths[0]])

        return pre_edges, post_edges

    def get_first_blocks(self):
        """
        :param graph: Graph
        :return: Returns a list of all blocks having no incoming edges
        """
        firsts = []
        for b in self.blocks:
            n_in = 0
            if len(self.reverse_adj_list[b]) == 0:
                firsts.append(b)

        return firsts

    def get_last_blocks(self):
        """
        :param graph: Graph
        :return: Returns a list of all blocks having no incoming edges
        """
        lasts = []
        for b in self.blocks:
            if len(self.adj_list[b]) == 0:
                lasts.append(b)

        return lasts

    @staticmethod
    def intervals_contains_rp(rp, intervals):
        """
        Returns true if any of the intervals cover rp
        """
        for i in intervals:
            if rp in i.region_paths:
                return True

        return False

    def create_subgraph_from_intervals(self, intervals, padding=10000, alt_locus=None, base_trans=None):
        """
        Creates a subgraph containing all the intervals
        :param intervals: list of intervals. All region paths in the intervals must create a connected subgraph.
        :param padding: number of baseapairs that should be in graph before first and after last intervals
        :return:
        """

        new_first = None
        new_last = None

        blocks = []
        for i in intervals:
            blocks.extend(i.region_paths)

        subgraph = self.create_subgraph_from_blocks(blocks, alt_locus=alt_locus)

        #if base_trans is None:
        trans = Translation({}, {}, graph=subgraph)
        #else:
        #trans = base_trans

        remove = []  # Blocks to remove in the end

        # Prune graph from beginning
        while True:
            intervals = [trans.translate(i) for i in intervals]
            # Find first block, has no genes and only one edge out: Delete
            # If contains genes: Prune
            # IF contains no genes, but multiple edges out: Prune
            first_rps = subgraph.get_first_blocks()
            assert len(first_rps) == 1 , "%s has not length 1" % (first_rps)
            first_rp = first_rps[0]

            start_position = Position(first_rp, 0)

            if not Graph.intervals_contains_rp(first_rp, intervals):
                # If only one edge out, delete this one
                if len(self.adj_list[first_rp]) == 1:
                    subgraph.remove(first_rp)
                    continue


            # Prune this rp
            first = first_rp
            first_length = subgraph.blocks[first].length()

            # Find lowest start in this region path
            first_start = first_length
            for i in intervals:
                if i.region_paths[0] == first:
                    first_start = min(i.start_position.offset, first_start)


            padding_first = max(1, (first_length - first_start) + padding)
            split_trans = Translation({}, {}, graph=subgraph)
            split_trans.graph2 = subgraph
            if first_length > padding:

                # Divide first block
                new_first = subgraph._next_id()
                new_first_length = first_length - padding_first
                new_second = subgraph._next_id()
                trans_first = Translation(
                    {first: [Interval(0, padding_first, [new_first, new_second])]},
                    {new_first: [Interval(0, new_first_length, [first], subgraph)],
                     new_second: [Interval(new_first_length, first_length,
                                           [first], subgraph)]}, graph=subgraph)

                start_position = Position(first, new_first_length)

                #print(trans_first._a_to_b)
                #print(trans_first)
                subgraph = trans_first.translate_subgraph(subgraph)
                split_trans = split_trans + trans_first
                #subgraph.remove(new_first)
                remove.append(new_first)
            subgraph = trans.translate_subgraph(subgraph)
            trans = trans + split_trans
            break


        # Prune from end
        while True:
            intervals = [trans.translate(i) for i in intervals]
            # Find first block, has no genes and only one edge out: Delete
            # If contains genes: Prune
            # IF contains no genes, but multiple edges out: Prune
            last_rps = subgraph.get_last_blocks()
            #assert len(last_rps) == 1 , "%s is more than one last rp" % last_rps
            last_rp = last_rps[0]

            if not Graph.intervals_contains_rp(last_rp, intervals):
                # If only one edge out, delete this one
                if len(self.reverse_adj_list[last_rp]) == 1:
                    subgraph.remove(last_rp)
                    continue


            # Prune this rp
            last = last_rp
            last_length = subgraph.blocks[last].length()

            # Find last end in end rp
            last_end = 0
            for i in intervals:
                if i.region_paths[-1] == last:
                    last_end = max(i.end_position.offset, last_end)

            padding_end = min(last_length - 1, last_end + padding)

            if last_length > padding:
                # Divide last block
                new_last = subgraph._next_id()
                new_last_length  = last_length - padding_end
                new_second = subgraph._next_id()
                trans_last = Translation({last: [Interval(0, padding_end, [new_second, new_last])]},
                                    {new_last: [Interval(padding_end, last_length, [last], subgraph)],
                                     new_second: [Interval(0, padding_end, [last], subgraph)]}, graph=subgraph)
                subgraph = trans_last.translate_subgraph(subgraph)

                trans = trans + trans_last
                #subgraph.remove(new_second)
                remove.append(new_last)
            break

        for r in remove:
            subgraph.remove(r)
        starts = subgraph.get_first_blocks()
        assert len(starts) == 1, " %s has not len 1" % (str(starts))

        #assert False
        #print("<p>Subgraph</p>")
        #print(subgraph)
        return subgraph, trans, start_position

    def create_subgraph_from_blocks(self, blocks, alt_locus=None):
        """
        Creates a subgraph using existing edges and only the blocks send as argument
        :param blocks: list of block ids
        :param alt_locus: If not None, alt loci blocks not from this alt locus will not be added to graph
        This wil typically be parallell alt loci that potentially can be added
        :return: Returns a new graph
        """
        blocks = set(blocks)
        # add prev and next critical

        new_edges = defaultdict(list)
        new_blocks = {}
        for b in blocks:
            new_blocks[b] = Block(self.blocks[b].length())

        for b in blocks:
            for e in self.adj_list[b]:
                #if e in new_blocks:
                if alt_locus is not None and Graph.block_origin(e) == "alt" and alt_locus not in e:
                    continue

                new_edges[b].append(e)
                if e not in new_blocks:
                    #print("<p>adding %s</p>" % e)
                    new_blocks[e] = Block(self.blocks[e].length())

        # Go through all added blocks, add edges into other added blocks
        for b in new_blocks:
            #print("Checking edges for %s" % (b))
            for edge in self.adj_list[b]:
                #print(" Edge %s" % edge)
                if edge in new_blocks and edge not in new_edges[b]:
                    new_edges[b].append(edge)

        subgraph = Graph(new_blocks, new_edges)
        subgraph_without_critical = subgraph.copy()

        # Add all blocks going into first blocks in graph
        firsts = subgraph.get_first_blocks()
        for f in firsts:
            for before in self.reverse_adj_list[f]:

                if alt_locus is not None and Graph.block_origin(before) == "alt" and alt_locus not in before:
                    continue

                new_blocks[before] = Block(self.blocks[before].length())
                #print("Adding before: %s" % before)
                new_edges[before].append(f)
                #print("  adding edge from %s to %s" % (before, f))

        subgraph = Graph(new_blocks, new_edges)

        # If two last rps, they should be going to the same next
        #print(subgraph)
        lasts = subgraph.get_last_blocks()
        assert len(lasts) == 1 or len(lasts) == 2, "%s is lasts" % lasts

        if len(lasts) == 2:
            # Connect the two last to next
            #print("Last")
            #print(lasts)
            if len(self.adj_list[lasts[0]]) > 0 and len(self.adj_list[lasts[1]]) > 0:
                assert self.adj_list[lasts[0]][0] == self.adj_list[lasts[1]][0], "Not identical next: %s != %s" % (self.adj_list[lasts[0]][0], self.adj_list[lasts[1]][0])
                next = self.adj_list[lasts[0]]
                assert len(next) == 1
                next = next[0]
                last = self.adj_list[lasts[0]][0]
                new_blocks[last] = Block(self.blocks[last].length())
                new_edges[lasts[0]].append(next)
                new_edges[lasts[1]].append(next)

        subgraph2 = Graph(new_blocks, new_edges)
        return subgraph2

        # Delete everything below here


        # Append with prev critical and next critical (only if first/last is not critical)

        """
        first_blocks = subgraph.get_first_blocks()
        assert len(first_blocks) > 0, "Subgraph has no first blocks. %s" % (subgraph)
        first = first_blocks[0]
        critical = self.find_critical_blocks(first)
        critical.append(self.get_last_blocks()[0])
        critical = set(critical)

        prev_critical = self.find_previous_critical_block(first, critical)
        new_blocks[prev_critical] = Block(self.blocks[prev_critical].length())
        """
        lasts = subgraph.get_last_blocks()

        assert len(lasts) == 1 or len(lasts) == 2

        if len(lasts) == 2:
            # Connect the two last to next
            assert self.adj_list[lasts[0]] == self.adj_list[lasts[1]], "Not identical next: %s != %s" % (self.adj_list[lasts[0]], self.adj_list[lasts[0]])
            next = self.adj_list[lasts[0]]
            assert len(next) == 1
            next = next[0]
            new_blocks.append(self.adj_list[lasts[0]])
            new_edges[lasts[0]].append(next)
            new_edges[lasts[1]].append(next)

        """
        if len(lasts) > 1:

            last = lasts[0]
            next_critical = self.find_next_critical_block(last, critical)
            print("Next critical block: %s" % next_critical)
            new_blocks[next_critical] = Block(self.blocks[next_critical].length())

            # Create new edges to the new blocks
            for l in subgraph_without_critical.get_last_blocks():
                print("Connecting edges from %s" % l)
                if l != next_critical and next_critical in self.adj_list[l]:
                    new_edges[l] = [next_critical]
                    print(" Connecting %s to %s" % (l, next_critical))


        for f in subgraph_without_critical.get_first_blocks():

            if f != prev_critical:
                new_edges[prev_critical] = [f]
        """
        subgraph2 = Graph(new_blocks, new_edges)

        return subgraph2

    def remove(self, block_id):
        """Remove a block including edges from the graph

        :param block_id: block id of the block

        """
        del self.blocks[block_id]
        for edge in self.adj_list[block_id]:
            self.reverse_adj_list[edge].remove(block_id)
        del self.adj_list[block_id]

        for edge in self.reverse_adj_list[block_id]:
            self.adj_list[edge].remove(block_id)
        del self.reverse_adj_list[block_id]

    def _add_edge(self, block_a, block_b):
        """Add edge from a to b, in both adj_list
        and reveres_adj_list

        :param block_a: from block
        :param block_b: to block
        """

        self.adj_list[block_a].append(block_b)
        self.reverse_adj_list[block_b].append(block_a)

    def _join_blocks(self, block_ids):
        """Merge the given blocks in the graph and
        giving them the combination of all edges.

        :param block_ids: block ids of blocks to merge
        :returns: translation object
        :rtype: Translation

        """
        if not block_ids:
            return Translation()
        lens = [self.blocks[_id].length() for _id in block_ids]
        L = lens[0]
        assert all(l == L for l in lens)

        # Find all incoming and outgoing edges
        in_edges = set()
        out_edges = set()
        for block_id in block_ids:
            in_edges.update(self.reverse_adj_list[block_id])
            out_edges.update(self.adj_list[block_id])

        # Add the new block
        new_id = self._next_id()
        self.blocks[new_id] = Block(L)
        for block_id in out_edges:
            self._add_edge(new_id, block_id)
        for block_id in in_edges:
            self._add_edge(block_id, new_id)

        # Remove old blocks
        for block_id in block_ids:
            self.remove(block_id)

        # Create translation object
        translation = Translation(
            {block_id: [Interval(Position(new_id, 0), Position(new_id, L))]
             for block_id in block_ids},
            {new_id:
             [Interval(Position(block_id, 0), Position(block_id, L))
              for block_id in block_ids]})

        return translation

    def _split_block(self, block_id, offsets):
        """Splits the given block at offsets
        and return translation object

        :param block_id: block to split
        :param offset: where to split
        :returns: translation object
        :rtype: Translation

        """
        l = self.blocks[block_id].length()
        if not offsets:
            return Translation()

        # Add blocks
        blocks = [Block(b-a) for a, b in zip([0]+offsets, offsets+[l])]
        ids = [self._next_id() for _ in blocks]
        self.blocks.update(dict(zip(ids, blocks)))

        # Add edges
        for id1, id2 in zip(ids[:-1], ids[1:]):
            self._add_edge(id1, id2)

        for e in self.adj_list[block_id]:
            self._add_edge(ids[-1], e)
        for e in self.reverse_adj_list[block_id]:
            self._add_edge(e, ids[0])

        self.remove(block_id)

        # Set translation
        return Translation(
            {block_id:
             [Interval(
                 Position(ids[0], 0),
                 Position(ids[-1], blocks[-1].length()),
                 ids)]},
            {_id: [Interval(Position(block_id, offset),
                            Position(block_id, offset+block.length()))]
             for _id, block, offset, in zip(ids, blocks, [0]+offsets)}
        )

    def _merge_clean_intervals(self, intervals):
        """Merge intervals that all start at beginning
        of a region path

        :param intervals: 
        :returns: 
        :rtype: Translation

        """
        first_rps = [interval.region_paths[0] for interval in intervals]
        first_rp_lengths = [self.blocks[rp].lenght()
                            for rp in first_rps]
        min_len = min(first_rp_lengths)
        translations = [self._split_block(rp, min_len) for
                        rp in first_rps]

    def get_split_translation(self, rp, offset):

        id_a, id_b = (self._next_id(), self._next_id())
        L = self.blocks[rp].length()
        trans_dict = {}
        reverse_dict = {}
        trans_dict[rp] = [Interval(
            Position(id_a, 0),
            Position(id_b, L-offset))]

        reverse_dict[id_a] = [Interval(Position(rp, 0),
                                       Position(rp, offset),
                                       graph=self)]

        reverse_dict[id_b] = [Interval(Position(rp, offset),
                                       Position(rp, L),
                                       graph=self)]

        trans = Translation(trans_dict, reverse_dict, graph=self)

        new_graph = trans.translate_subgraph(self)
        trans.graph2 = new_graph
        for i_list in trans_dict.values():
            for i in i_list:
                i.graph = new_graph

        return new_graph, trans

    def connect_postitions(self, position_a, position_b):
        """Connect position_a to position_b with an edge
        and split the region paths if necessary.
        Create translation object and new graph
        :param position_a: Position
        :param position_b: Position
        :returns: New graph and Translation object
        :rtype: Graph, Translation

        """
        rp_a = position_a.region_path_id
        rp_b = position_b.region_path_id
        offset_a = position_a.offset
        offset_b = position_b.offset
        l_a = self.blocks[rp_a].length()
        l_b = self.blocks[rp_b].length()
        a_to_b = {}
        b_to_a = {}
        rp_out = rp_a
        rp_in = rp_b
        assert offset_a < l_a
        assert offset_b < l_b, "%s, %s" % (offset_b, l_b)

        if offset_a < l_a - 1:
            # Split first block
            idf = self._next_id()
            idl = self._next_id()
            rp_out = idf
            a_to_b[rp_a] = [Interval(0, l_a-offset_a-1,
                                     [idf, idl])]
            b_to_a[idf] = [Interval(0, offset_a+1, [rp_a])]
            b_to_a[idl] = [Interval(offset_a+1, l_a, [rp_a])]

        if offset_b > 0:
            # Split last block
            idf = self._next_id()
            idl = self._next_id()
            rp_in = idl
            a_to_b[rp_b] = [Interval(0, l_b-offset_b, [idf, idl])]
            b_to_a[idf] = [Interval(0, offset_b, [rp_b])]
            b_to_a[idl] = [Interval(offset_b, l_b, [rp_b])]

        trans = Translation(a_to_b, b_to_a, graph=self)
        n_graph = trans.translate_subgraph(self)
        n_graph.adj_list[rp_out].append(rp_in)
        n_graph.reverse_adj_list[rp_in].append(rp_out)
        trans.graph2 = n_graph
        self._update_a_b_graph(trans._a_to_b, n_graph)
        self._update_a_b_graph(trans._b_to_a, self)
        return n_graph, trans

    def prev_position(self, pos):
        """
        :param pos: Position
        :return: Returns the previous position in graph
        """
        if pos.offset > 0:
            return Position(pos.region_path_id, pos.offset - 1)
        else:
            # Find previous region path
            prev = None
            for b in self.blocks:
                if pos.region_path_id in self.adj_list[b]:
                    prev = b
                    break
            if prev is None:
                raise Exception("Found no previous pos for position %s" % pos)

            return Position(b, self.blocks[b].length() - 1)

    def next_position(self, pos):
        """
        :param pos: Position
        :return: Returns the previous position in graph
        """
        if pos.offset < self.blocks[pos.region_path_id].length() - 1:
            return Position(pos.region_path_id, pos.offset + 1)
        else:
            # Find next region path
            next = self.adj_list[pos.region_path_id][0]
            return Position(next, 0)

    def _get_insulated_merge_transformation(self, intervals):
        """Merge intervals that all start and end at region path
        boundries

        :param intervals:
        :returns: translation and new graph
        :rtype: (Translation, Graph)

        """
        first_rps = [i.region_paths[0] for i in intervals]
        offset = min(self.blocks[rp].length() for rp in first_rps)
        if offset == intervals[0].length():
            return Translation(graph=self), self

        translation = Translation(graph=self)
        graph = self
        for rp in first_rps:
            if self.blocks[rp].length() == offset:
                continue
            trans, graph = graph.get_split_translation(rp, offset)
            translation += trans

        intervals = [translation.translate(i) for i in intervals]
        first_rps = [i.region_paths[0] for i in intervals]

        merge_id = graph._next_id()
        a_to_b = {rp: [Interval(Position(merge_id, 0),
                                Position(merge_id, offset),
                                )]
                  for rp in first_rps}
        b_to_a = {merge_id: [Interval(Position(rp, 0),
                                      Position(rp, offset),
                                      graph=graph)]}

        merge_trans = Translation(a_to_b, b_to_a, graph=graph)
        intervals = [merge_trans.translate(i) for i in intervals]
        full_trans = translation + merge_trans
        new_graph = full_trans.translate_subgraph(graph)
        trans, ng = new_graph._get_insulated_merge_translation(intervals)
        return full_trans+trans, ng

    def _get_inslulate_translation(self, intervals, block_lengths = {}):
        """Get translation for splitting the region paths
        at the start and end of the intervals such that
        all the intervals span complete region paths

        :param intervals:
        :returns:
        :rtype: Translation

        """
        for i in intervals:
            self.assert_interval_in_graph(i)
        # Split starts
        starts = [i.start_position for i in intervals]
        translation = Translation(graph=self)
        translation.graph2 = self
        cur_graph = self
        i = 0
        id_a = None
        id_b = None

        for start in starts:
            i += 1
            if start.offset == 0:
                continue
            rp = start.region_path_id
            offset = start.offset
            id_a, id_b = (self._next_id(), self._next_id())
            L = self.blocks[rp].length()
            prev_graph = cur_graph.copy()

            trans_dict = {}
            reverse_dict = {}
            trans_dict[rp] = [Interval(
                Position(id_a, 0),
                Position(id_b, L-offset))]

            reverse_dict[id_a] = [Interval(Position(rp, 0),
                                           Position(rp, offset),
                                           graph=prev_graph)]
            block_lengths[id_a] = offset

            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L),
                                           graph=prev_graph)]
            block_lengths[id_b] = L - offset

            tmp_trans = Translation(trans_dict, reverse_dict, graph=cur_graph)
            cur_graph = tmp_trans.translate_subgraph(cur_graph)
            self._update_a_b_graph(trans_dict, cur_graph)

            translation += tmp_trans

        if id_a is not None and id_b is not None:
            assert id_a in cur_graph.blocks
            assert id_b in cur_graph.blocks

        # Translate intervals to new graph
        new_intervals = [translation.translate(interval)
                         for interval in intervals]
        ends = [interval.end_position for interval in new_intervals]

        # Update interval's graph. Is this necessary? Should be done in translate interval
        for interval in new_intervals:
            interval.graph = cur_graph

        back_graph_end = cur_graph
        end_translations = Translation(graph=cur_graph)
        end_translations.graph2 = back_graph_end

        # Split end of intervals in new graph
        for end in ends:
            assert cur_graph is not None
            assert end.region_path_id in cur_graph.blocks
            rp = end.region_path_id
            L = cur_graph.blocks[rp].length()
            offset = end.offset
            if offset == L:
                continue

            id_a, id_b = (self._next_id(), self._next_id())

            trans_dict = {}
            reverse_dict = {}
            trans_dict[rp] = [Interval(
                Position(id_a, 0),
                Position(id_b, L-offset)
                )]

            prev_graph = cur_graph#.copy()
            reverse_dict[id_a] = [Interval(Position(rp, 0),
                                           Position(rp, offset),
                                           graph=prev_graph)]
            block_lengths[id_a] = offset

            reverse_dict[id_b] = [Interval(Position(rp, offset),
                                           Position(rp, L),
                                           graph=prev_graph)]
            block_lengths[id_b] = L - offset

            tmp_trans = Translation(trans_dict, reverse_dict, graph=prev_graph)
            cur_graph = tmp_trans.translate_subgraph(cur_graph)

            self._update_a_b_graph(tmp_trans._a_to_b, cur_graph)

            tmp_trans.graph2 = cur_graph

            # Asssert translations now have all intervals needed
            self._assert_translation_intervals_has_graphs(translation)

            end_translations += tmp_trans


        translation += end_translations

        return translation, cur_graph

    def _assert_translation_intervals_has_graphs(self, translation):
        for intervals in translation._a_to_b.values():
            for interval in intervals:
                assert interval.graph is not None

        for intervals in translation._b_to_a.values():
            for interval in intervals:
                assert interval.graph is not None

    def _update_a_b_graph(self, a_b, graph, prev_graph = None):
        for block, intvs in a_b.items():
            for interval in intvs:
                interval.graph = graph
                if prev_graph != None:
                    interval.set_length_cache(prev_graph.blocks[block].length())

    def assert_position_in_graph(self, position, exclusive=False):
        assert position.region_path_id in self.blocks
        block = self.blocks[position.region_path_id]
        assert position.offset < block.length() + int(exclusive)

    def assert_interval_in_graph(self, interval):
        for rp in interval.region_paths:
            assert rp in self.blocks, "%s not in %s" % (rp, self.blocks.keys())
        self.assert_position_in_graph(interval.start_position)
        self.assert_position_in_graph(interval.end_position, exclusive=True)

    def merge(self, intervals):
        """Merges the given intervals in the graph, and returns
        a the resulting graph after merge and a translation object.
        :param intervals: list of intevals to merge
        :returns: translation object, resulting graph
        :rtype: Graph, Translation
        """

        block_lengths = self.blocks.copy()
        for b in block_lengths:
            block_lengths[b] = self.blocks[b].length()

        [self.assert_interval_in_graph(i) for i in intervals]
        original_graph = intervals[0].graph
        # Assume same lengths for all intervals
        length = intervals[0].length()
        for interval in intervals[1:]:
            assert interval.length() == length, \
                "All intervals should have the same length (%d, %d)" % (interval.length(), length)


        # 1: Translate graph so that all intervals starts and ends at rps
        trans, graph1 = self._get_inslulate_translation(intervals, block_lengths)
        trans.graph2 = graph1
        self._update_a_b_graph(trans._a_to_b, graph1)  # correct, a to b interval's graph is wrong for some reason

        # Update intervals to be on graph1
        new_intervals = []
        for interval in intervals:
            mp = trans.translate_interval(interval)
            new_interval = mp.get_single_path_intervals()[0]
            new_interval.graph = graph1
            new_intervals.append(new_interval)

        intervals = new_intervals

        # Step 2: Translate each interval to one single block
        # In next graph, each interval is replaced by one single block
        a_b = {}
        b_a = {}
        new_blocks = []
        for interval in intervals:
            # Create one new large block
            large_block = graph1._next_id()
            new_blocks.append(large_block)

            # Back translation from block to interval
            b_a[large_block] = [interval]
            block_lengths[large_block] = interval.length()

            # Find forward translations
            prev_offset = 0
            offset = 0
            for rp in interval.region_paths:
                offset += graph1.blocks[rp].length()
                a_b[rp] = [Interval(prev_offset, offset, [large_block])]
                a_b[rp][0].set_length_cache(offset - prev_offset)
                prev_offset = offset

        # Find new graph
        trans2 = Translation(a_b, b_a, graph1, block_lengths)
        #graph2 = trans2.translate_subgraph(graph1)
        #trans2.graph2 = graph2

        # Update graph in a_b intervals to graph2 (which is now created)
        #self._update_a_b_graph(a_b, graph2)

        # Step 3: Merge each block (representing one interval each) to one single block
        new_block = graph1._next_id()
        a_b = {}
        b_a = {new_block: []}

        for block in new_blocks:
            a_b[block] = [Interval(0, length, [new_block])]
            a_b[block][0].set_length_cache(length)
            i = Interval(0, length, [block])
            b_a[new_block].append(i)
            i.set_length_cache(length)

        block_lengths[new_block] = length

        # Translate graph
        trans3 = Translation(a_b, b_a, None, block_lengths)
        #graph3 = trans3.translate_subgraph(graph2)
        #self._update_a_b_graph(a_b, graph3)
        #trans3.graph2 = graph3

        # Step 4: Divide large middle block
        starts = []
        for interval in intervals:
            prev_start = 0
            for rp in interval.region_paths:
                prev_start += interval.graph.blocks[rp].length()
                starts.append(prev_start)

        starts = list(set(starts))
        starts.sort()


        # Create translation from big block to all small created from each start
        a_b = {new_block: []}
        b_a = {}
        prev_start = 0
        new_blocks = []
        last_length = 0
        sum_len = 0
        for start in starts:
            new_small = graph1._next_id()  # graph3._next_id()
            new_blocks.append(new_small)
            i = Interval(prev_start, start, [new_block]) #, graph3)
            i.set_length_cache(start - prev_start)
            b_a[new_small] = [i]
            block_lengths[new_small] = start - prev_start
            last_length = start - prev_start
            sum_len += (start - prev_start)
            prev_start = start

        i2 = Interval(0, last_length, new_blocks)
        i2.set_length_cache(sum_len)
        a_b[new_block].append(i2)


        trans4 = Translation(a_b, b_a, None, block_lengths)
        #print("==== Trans 4 ====")
        #print(trans4)
        #graph4 = trans4.translate_subgraph(graph3)
        #self._update_a_b_graph(a_b, graph4)
        #trans4.graph2 = graph4


        final_trans = trans + trans2
        final_trans += trans3
        final_trans += trans4
        """

        final_trans = trans3 + trans4
        final_trans = trans2 + final_trans
        final_trans = trans + final_trans
        """
        final_trans.graph1 = original_graph  # Should not be needed ?!


        final_graph = final_trans.translate_subgraph(self)
        final_trans.graph2 = final_graph

        return final_graph, final_trans

    def _split_blocks_at_starts_and_ends(self, intervals):
        """
        Splits the region paths at the starts and ends of the
        intervals, such that all the intervals starts at the
        beginning of the intervals

        :param intervals: list(intervals)
        :returns: full translation for all block splits
        :rtype: Translation

        """
        full_translation = Translation()

        # Split starts
        for interval in intervals:
            offset = interval.start_position.offset
            if offset == 0:
                continue
            rp = interval.region_paths[0]
            full_translation += self._split_block(rp, [offset])

        for interval in intervals:
            end_offset = interval.end_position.offset
            rp = interval.end_position.region_path_id
            L = self.blocks[rp]
            if end_offset == L:
                continue
            full_translation += self.split_block(rp, [offset])

        return full_translation

    @takes(Interval, Interval)
    def merge_intervals(self, interval_a, interval_b, copy=True):
        """
        Merge the two intervals in the graph and return
        a translation object

        :param interval_a: One interval
        :param interval_b: One interval
        :returns: translation from old to new graph
        :rtype: Translation

        """
        assert interval_a.length() == interval_b.length()
        assert not any(rp_a in interval_b.region_paths
                       for rp_a in interval_a.region_paths)
        break_points = self._get_all_block_borders(
            interval_a,
            interval_b)

        sub_intervals_a = interval_a.split(break_points)
        sub_intervals_b = interval_b.split(break_points)

        ids = (self._next_id() for _ in sub_intervals_a)
        translation = self._generate_translation(
            interval_a, interval_b,
            sub_intervals_a, sub_intervals_b,
            ids)

        graph = self.copy() if copy else self
        # for pair in zip(sub_intervals_a, sub_intervals_b):
        #     start_intervals = [i for i in pair if i.start_position.offset == 0]
        #     in_edges = sum([self.reverse_adj_list[i.region_paths[0]]
        #                     for i in start_intervals], [])
        #     out_intervals = [
        #         i for i in pair if
        #         i.end_position.offset == self.blocks[i.end_position.region_path_id].length()]
        #     out_edges = sum([self.adj_list[i.region_paths[-1]] for i in end_intervals], [])
        # 
        #     graph

        return translation

    def _get_all_block_borders(self, interval_a, interval_b):
        """Return a list of block changes in both a and b

        :param interval_a: Interval
        :param interval_b: Interval
        :returns: list of all break_points
        :rtype: list

        """
        break_points = []
        for interval in [interval_a, interval_b]:
            offset = -interval.start_position.offset
            for region_path in interval.region_paths[:-1]:
                offset += self.blocks[region_path].length()
                break_points.append(offset)
            offset += interval.end_position.offset
            break_points.append(offset)

        unique_points = sorted(set(break_points))
        return unique_points

    @takes(Interval, Interval)
    def connect_intervals(self, interval_a, interval_b):
        pass

    def __str__(self):
        return "Graph: \n Blocks: %s\n Edges: %s" % (self.blocks, self.adj_list)

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def _get_reverse_edges(adj_list):
        reverse_edges = defaultdict(list)
        for block, edges in adj_list.items():
            for edge in edges:
                reverse_edges[edge].append(block)

        return reverse_edges

    def __eq__(self, other):

        if(len(self.blocks) != len(other.blocks)):
            return False

        if self.blocks != other.blocks:
            print("Different blocks")
            return False

        for adj in self.adj_list:
            #if adj in other.adj_list and self.adj_list[adj] != other.adj_list[adj]:
            if set(self.adj_list[adj]) != set(other.adj_list[adj]):
                #print("Different adj list for key %d" % (adj))
                return False
        for adj in other.adj_list:
            if set(self.adj_list[adj]) != set(other.adj_list[adj]):
            #if adj in self.adj_list and self.adj_list[adj] != other.adj_list[adj]:
                #print("Different adj list 2")
                return False


        #if self.adj_list != other.adj_list:
        #    return False



        return True

    @staticmethod
    def is_main_name(name):
        if "alt" not in name:
            return True
        if name.count("chr") > 1:
            return True
        return False

    @staticmethod
    def block_origin(name):
        """
        Checks if block is merged, alt or main based on the name
        :param name: block name
        :return: merged, alt or main
        """
        name = str(name)
        if name.count("chr") == 1 and not "alt" in name:
            return "main"
        elif name.count("chr") == 2:
            return "merged"
        elif "alt" in name:
            return "alt"
        return "main"

    @staticmethod
    def level_dict(blocks):
        level_mapping = {"alt": 0,
                         "main": 2,
                         "merged": 1}
        out = {}
        for b in blocks:
            out[b] = level_mapping[Graph.block_origin(b)]
        return out

    def find_critical_blocks(self, start_block):
        """Find all critical_blocks starting from
        start block

        :param start_block: block id of start_block
        :returns: List of critical block ids
        :rtype: list(str)

        """
        cur_block = start_block
        counter = 0
        critical_blocks = []
        visited = []

        while(self.adj_list[cur_block]):
            visited.append((cur_block, counter))
            # assert counter >= 0, visited
            if counter == 0:
                critical_blocks.append(cur_block)
            nexts = self.adj_list[cur_block]
            counter += (len(nexts) - 1)
            cur_block = nexts[0]
            counter -= (len(self.reverse_adj_list[cur_block])-1)

        if (counter == 0):
            critical_blocks.append(cur_block)

        return critical_blocks

    def find_previous_critical_block(self, block, critical_blocks=None):
        """Find previous critical block, going backwards from block.
        Traverse blocks and return once a block in critical_blocks
        is found.

        :param block: block to start from
        :param critical_blocks: list of critical blocks
        :returns: previous critical block
        :rtype: str

        """
        if critical_blocks is None:
            critical_blocks = self.critical_blocks

        if block in critical_blocks:
            return block
        cur_block = block
        while self.reverse_adj_list[cur_block]:
            prevs = self.reverse_adj_list[cur_block]
            if not prevs:
                raise Exception("No prevs: %s" % prevs)

            prev_block = prevs[0]

            if prev_block in critical_blocks:
                return prev_block

            cur_block = prev_block

        raise Exception("Did not find critical node from %s" % block)

    def find_next_critical_block(self, block, critical_blocks):
        if block in critical_blocks:
            return block
        cur_block = block
        while self.adj_list[cur_block]:
            prevs = self.adj_list[cur_block]
            prev_mains = [b for b in prevs]
            if not prev_mains:
                raise Exception("No prevs: %s" % prevs)

            prev_main_block = prev_mains[0]

            if prev_main_block in critical_blocks:
                return prev_main_block
            cur_block = prev_main_block

        raise Exception("No critical blocks found")

    def are_paralell(self, block_a, block_b, critical_blocks=None):
        if critical_blocks is None:
            critical_blocks = self.critical_blocks
        if block_a == block_b:
            return True
        critical_block_a = self.find_previous_critical_block(block_a, critical_blocks)
        critical_block_b = self.find_previous_critical_block(block_b, critical_blocks)
        return critical_block_a == critical_block_b

    def find_parallell_blocks(self, blocks, filter_func):
        """Find all parallell_blocks to block
        that satisfies filter_func

        :param block: region path id
        :param filter_func: predicate fucntion
        :returns: all parallell_blocks
        :rtype: list(str)

        """
        if isinstance(blocks, str):
            blocks = [blocks]
        prev_blocks = [b for b in self.reverse_adj_list[blocks[0]]
                       if filter_func(b)]
        assert len(prev_blocks) == 1, blocks
        start_block = prev_blocks[0]
        next_blocks = [b for b in self.adj_list[blocks[-1]]
                       if filter_func(b)]
        assert len(next_blocks) == 1, blocks
        end_block = next_blocks[0]
        cur_block = start_block
        parallell_blocks = []
        while True:
            next_blocks = [b for b in self.adj_list[cur_block]
                           if filter_func(b)]
            assert len(next_blocks) == 1
            cur_block = next_blocks[0]
            if cur_block == end_block:
                break
            parallell_blocks.append(cur_block)
        return parallell_blocks

    def find_all_critical_blocks(self):
        """Find critical blocks in the graph.
        I.e.  blocks that are traversed by all paths

        :returns: list of critical_blocks ids
        :rtype: list(str)

        """

        start_blocks = [block for block in self.blocks if
                        not self.reverse_adj_list[block]]

        critical_blocks = []
        for start_block in start_blocks:
            critical_blocks.extend(
                self.find_critical_blocks(start_block))

        return critical_blocks

    def n_edges_in(self, block):
        """
        Finds and returns the number of edges going in to a block
        :param block:
        :return: Returns the number of edges
        """
        n = 0
        for b in self.blocks:
            if block in self.adj_list[b]:
                n += 1
        return n

    def has_identical_structure(self, other):
        """
        Checks if this graph has identical
        structure (edges and blocks) to other graph.
        Size of region paths is ignores (and can be different).
        :param other: Graph to compare with
        :return: True if identical, otherwise False
        """

        if len(self.blocks) != len(other.blocks):
            print("Different number of blocks")
            return False

        # For every block, check that there exists
        # a block in other graph with the same number of
        # edges in and out
        other_blocks = list(other.blocks.keys()).copy()
        for b in self.blocks:
            #print("Checking %d" % b)
            match = False
            for ob in other_blocks:
                #print("   Checking %d" % ob)
                sim_out = len(self.adj_list[b]) == len(other.adj_list[ob])
                sim_in = self.n_edges_in(b) == other.n_edges_in(ob)
                if sim_out and sim_in:
                    # Remove from list to check, and check next (break)
                    other_blocks.remove(ob)
                    #print("      Match!")
                    match = True
                    break
            if not match:
                # No match for block b, return False
                #print("No match for block %d" % (b))
                return False

        return True
