from .interval import Interval, Position
from .multipathinterval import GeneralMultiPathInterval, \
    SingleMultiPathInterval, SimpleMultipathInterval
import pickle


class Translation(object):
    """
    Holds information necessary to translate coordinates and intervals
    from one graph to another

    >>> forward_dict = {"chr1": [Interval(0, 100, ["chr1-0", "chr1-1", "chr1-2"])],
    "alt": [Interval(0,100, ["chr1-1, "alt-0"])]}
    >>> reverse_dict = {"chr1-0": [Interval(0, 50, ["chr1"])],
    "chr1-1": [Interval(50, 110, ["chr1"]), Interval(0,50, ["alt"])],
    "chr1-2": [Interval(110, 210, ["chr1"])],
    "alt-0": [Interval(50, 100, ["alt"])]}
    >>> translation = Translation(forward_dict, reverse_dict, graph)

    """

    def __init__(self, translate_dict={}, reverse_dict={},
                 graph=None, block_lengths=None):
        """
        Init the translation object with two dicts. Each dict has
        region path IDs as keys and a list of intervals as values.

        :param translate_dict: Mapping of region paths in graph1 to intervals
        in graph2
        :param reverse_dict: Mapping of region paths in graph2 to intervals
        in graph 1
        """
        self.block_lengths = block_lengths
        self._a_to_b = translate_dict
        self._b_to_a = reverse_dict
        for k, intervals in translate_dict.items():
            for interval in intervals:
                assert interval.start_position != interval.end_position,\
                                                  "Empty interval in translate: (%s: %s)" % ("T", translate_dict)

        for k, intervals in reverse_dict.items():
            for interval in intervals:
                assert interval.start_position != interval.end_position,\
                                                  "Empty interval in translate: (%s: %s)" % ("R", reverse_dict)

        # Store graph references
        self.graph1 = None
        self.graph2 = None
        if len(self._a_to_b) > 0:
            self.graph2 = list(self._a_to_b.values())[0][0].graph
        if len(self._b_to_a) > 0:
            self.graph1 = list(self._b_to_a.values())[0][0].graph

        if graph is not None:
            self.graph1 = graph
        if self.graph1 and self.graph1.blocks:
            self.block_cls = list(self.graph1.blocks.values())[0].__class__
        else:
            from .graph import Block
            self.block_cls = Block

    @classmethod
    def make_name_translation(cls, trans_dict, graph):
        """
        Creates a copied version of trans_dict where region path have name IDs

        :param trans_dict:
        :param graph:
        :return: Returns a new translation object
        """
        for k in trans_dict.keys():
            assert k in graph.blocks, "%s not in %s" % (k, graph)

        rev_dict = {v: [Interval(Position(k, 0),
                                 Position(k, graph.blocks[k].length()),
                                 [k], graph)]
                    for k, v in trans_dict.items()}

        interval_dict = {k: [Interval(Position(v, 0),
                                      Position(v, graph.blocks[k].length()))]
                         for k, v in trans_dict.items()}

        return cls(interval_dict, rev_dict, graph)

    def to_file(self, file_name):
        """Writes translation object to pickle

        :param file_name: File name
        """
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def from_file(file_name):
        """Reads Translation object from pickle"""
        with open("%s" % file_name, "rb") as f:
            return pickle.loads(f.read())

    def translate(self, obj):
        """Check type of obj and translate accordingly
        Only for use where there is a one-to-one translation

        :param obj: Interval/Position
        :returns: translated object
        :rtype: Interval/Position

        """
        if isinstance(obj, Interval):
            ret = self.translate_interval(obj).get_single_path_intervals()
        elif isinstance(obj, Position):
            ret = self.translate_position(obj)
        else:
            raise ValueError("""Cannot translate object of type %s.
                Only Position and Interval are supportet""" % type(obj))
        assert len(ret) == 1
        return ret[0]

    def get_internal_edges(self, subgraph, edges, reverse_edges):
        """Find new edges gotten from splitting region paths
        under current translation

        :param subgraph: former subgraph
        :param edges: adjacency list to be updated
        :param reverse_edges: reverse adjacency list to be updated
        """
        for a in self._a_to_b:
            t_a = self._translations(a, inverse=False)
            for interval in t_a:
                for rp1, rp2 in zip(interval.region_paths[:-1],
                                    interval.region_paths[1:]):

                    assert rp1 != rp2
                    edges[rp1].append(rp2)
                    reverse_edges[rp2].append(rp1)

    def get_old_edges(self, subgraph):
        """Find all edges in subgraph that are still valid
        after this translation

        :param subgraph: former subgraph
        :returns: Adjacancy list of valid edges
        :rtype: defaultdict(list)

        """
        new_edges = subgraph.adj_list.copy()
        new_rev_edges = subgraph.reverse_adj_list.copy()
        for k, v in new_edges.items():
            new_edges[k] = v[:]
        for k, v in new_rev_edges.items():
            new_rev_edges[k] = v[:]

        for a in self._a_to_b:
            if a in new_edges:
                del new_edges[a]
            if a in new_rev_edges:
                del new_rev_edges[a]

        for a in self._a_to_b:
            for v in subgraph.reverse_adj_list[a]:
                if a not in new_edges[v]:
                    continue
                new_edges[v] = new_edges[v][:]
                new_edges[v].remove(a)
            for v in subgraph.adj_list[a]:
                if a not in new_rev_edges[v]:
                    continue
                new_rev_edges[v] = new_rev_edges[v][:]
                new_rev_edges[v].remove(a)

        return new_edges, new_rev_edges

    def get_external_edges(self, subgraph, edges, rev_edges):
        """Get new external edges, i.e edges between
        region paths that have been translated

        :param subgraph: subgraph to be translated
        :param edges: Adj list to update
        :rtype: None (updates edges)

        """
        for a in self._a_to_b:
            translations = self._translations(a, False)

            # Edges from a
            lasts = [i.region_paths[-1] for i in translations]
            for b in subgraph.adj_list[a]:
                firsts = [b]
                if b in self._a_to_b:
                    firsts = [i.region_paths[0]
                              for i in self._translations(b, False)]
                for l in lasts:
                    edges[l].extend([f for f in firsts if f != l])
                for f in firsts:
                    rev_edges[f].extend([l for l in lasts if l != f])

            # Edges to a
            my_firsts = [i.region_paths[0] for i in
                         translations]
            for b in subgraph.reverse_adj_list[a]:
                lasts = [b]
                if b in self._a_to_b:
                    lasts = [i.region_paths[-1]
                             for i in self._translations(b, False)]
                for l in lasts:
                    edges[l].extend([f for f in my_firsts if f != l])
                for f in my_firsts:
                    rev_edges[f].extend([l for l in lasts if f != l])

    def translate_subgraph(self, subgraph):
        """
        Translates a graph (forward). The graph has to be a subgraph of
        graph1 in the translation object.

        :param subgraph: Subgraph to translate.
        :return: Returns the translated subgraph
        """
        assert self.graph1 is not None,\
            "Cannot translate subgraph if graph1 is None"

        """
        Algorithm: For every interval crossing edge in subgraph, translate.
        Add every region path in translated interval and every edge
        Also, translate every region path and add edges and region paths found
        """
        edge_list_add = []

        edges, rev_edges = self.get_old_edges(subgraph)
        self.get_external_edges(subgraph, edges, rev_edges)
        self.get_internal_edges(subgraph, edges, rev_edges)
        for k, v in edges.items():
            edges[k] = list(set(v))
        for k, v in rev_edges.items():
            rev_edges[k] = list(set(v))

        new_blocks, edge_list_add = self._translate_subgraph_blocks(
            subgraph, [])
        return self.graph1.__class__(new_blocks, edges, rev_adj_list=rev_edges)

    def translate_interval(self, interval, inverse=False):
        """
        Translate an interval between the two coordinate systems.

        :param interval: Interval
        :param inverse: wheter to use the inverse translation
        :return: Returns an interval. If inverse is True,
            a list of intervals is returned.
        """
        """
        Algorithm:
            Translate start to new start
            Translate end to new end
            Find all new region paths
            Put this into a general multipath interval
        """
        trans_dict = self._b_to_a if inverse else self._a_to_b

        if not any(rp in trans_dict for rp in interval.region_paths):
            return SingleMultiPathInterval(interval,
                                           self._get_other_graph(inverse))

        if (interval.length() == 0):
            start_poses = self.translate_position(interval.
                                                  start_position,
                                                  inverse)
            return SimpleMultipathInterval(
                [Interval(sp, sp, [sp.region_path_id],
                          graph=self._get_other_graph(inverse),
                          direction = interval.direction)
                 for sp in start_poses])

        return self._translate_interval_nontrivial(interval, inverse)

    def translate_position(self, position, inverse=False):
        """
        Translates a position

        :param position: Position
        :param inverse: If True, translate back
        :return: Returns the translated Position
        :rtype: Position
        """
        trans_dict = self._b_to_a if inverse else self._a_to_b
        if position.region_path_id not in trans_dict:
            return [position]

        # Get interval for region path. Select first region path. Count offset.
        intervals = self._translations(position.region_path_id, inverse)
        positions = []
        for interval in intervals:
            if self.block_lengths is not None:
                rp_lens = [self.block_lengths[rp] for
                           rp in interval.region_paths]

            else:
                if interval.rp_lens_tmp is not None:
                    rp_lens = interval.rp_lens_tmp
                else:
                    rp_lens = [self._translations(rp, inverse=not inverse)[0].length()
                           for rp in interval.region_paths]
                    interval.rp_lens_tmp = list(rp_lens)

            found_pos = interval.get_position_from_offset(
                position.offset, rp_lens)
            positions.append(found_pos)

        return positions

    def interval_has_translation(self, interval, inverse=False):
        for rp in interval.region_paths:
            if not self.region_path_has_translation(rp, inverse):
                return False
        return True

    def region_path_has_translation(self, rp_id, inverse=False):
        """
        :return: Returns true if the region path has a translation
        """
        if not inverse and rp_id in self._a_to_b:
            return True
        elif inverse and rp_id in self._b_to_a:
            return True

        if not inverse and self.graph2 is None:
            raise Exception("Graph2 cannot be None when checking if forward translation of rp has translation")
        if inverse and self.graph1 is None:
            raise Exception("Graph1 cannot be None when checking if inverse translation of rp has translation")

        g = self.graph1 if inverse else self.graph2
        if rp_id in g.blocks:
            return True

        return False

    def __add__(self, other):
        """
        Combine (add) two translations.

        :param other: Another translation
        :return: Combined Translation
        :rtype: Translation
        """
        """
        Algorithm
            FOr every (region path => inter) in a, translate inter to c
            by using other.translate_interval(inter). Replace in a.
            For every (region path => inter) in c, translate recursviely back
            to all intervals using a.translate(inter, inverse). Replace in c
        """
        assert self.graph1 is not None,\
            "Graph1 cannot be 1 when adding translations"

        for intervals in self._b_to_a.values():
            for interval in intervals:
                assert interval.graph == self.graph1, \
                    "interval %s has graph %s \n!=\n %s \nTranslation: %s" % (
                        interval, interval.graph, self.graph1, self)

        new_trans = Translation(self._a_to_b, self._b_to_a, graph=self.graph1)
        region_paths = set(self._a_to_b.keys()).union(set(other._a_to_b.keys()))
        valid_region_paths = region_paths.intersection(set(self.graph1.blocks))

        # Find new forward translation dict
        new_translate_dict = {}
        for t in valid_region_paths:
            translated = other.translate_interval(
                self._translations(t)[0])
            new_translate_dict[t] = [
                i.copy() for i in translated.get_single_path_intervals()]

        new_trans._a_to_b = new_translate_dict
        i_region_paths = set(self._b_to_a.keys()).union(set(other._b_to_a.keys()))        
        valid_i_region_paths = [rp for rp in i_region_paths
                                if rp not in other._a_to_b]

        new_translate_dict = {}
        for t in valid_i_region_paths:
            translated = []
            for inter in other._translations(t, inverse=True):
                intervals = self.translate_interval(
                    inter, inverse=True).get_single_path_intervals()
                for translated_inter in intervals:
                    translated.append(translated_inter.copy())

            new_translate_dict[t] = translated

        new_b_to_a = new_translate_dict

        # Set b_to_a graph to graph1
        copied_graph = self.graph1
        for intvs in new_b_to_a.values():
            for intv in intvs:
                intv.graph = copied_graph

        new_trans.graph1 = copied_graph

        new_trans._b_to_a = new_b_to_a
        new_trans._assert_is_valid()

        # Important that intervals graphs match trans graph
        if len(list(new_trans._b_to_a.values())) > 0:
            assert(new_trans.graph1 == list(new_trans._b_to_a.values())[0][0].graph)

        # Find graph2 for new_trans
        graph2 = new_trans.translate_subgraph(new_trans.graph1)
        new_trans.graph2 = graph2

        # Update a_b intervals
        for intvs in new_trans._a_to_b.values():
            for intv in intvs:
                intv.graph = graph2

        return new_trans

    def translate_rp(self, rp, inverse=False):
        """Translate region path, by dict lookup

        :param rp: region path id
        :returns: translated interval
        :rtype: list(Interval)

        """
        if rp in self._a_to_b:
            return self._a_to_b[rp]

        rp_interval = Interval(Position(rp, 0),
                               Position(rp, self.graph1.blocks[rp].length()),
                               graph=self.graph1)
        return [rp_interval]

    def set_graph2(self, graph2):
        """Set graph2 ("other graph") and update all Interval objects
        in forward dict so that they have this graph.

        :param graph2: Graph
        """
        self.graph2 = graph2
        for k, v in self._a_to_b.items():
            for interval in v:
                interval.graph = graph2

    def __eq__(self, other):
        """Check if two translations are equal i.e
        both forward and reverse dicts are equal"""
        eq = self._a_to_b == other._a_to_b
        eq &= self._b_to_a == other._b_to_a
        return eq

    def copy(self):
        """Make new Translation object with copies
        of forward and reverse translation dicts

        :rtype: Translation

        """

        new_ab = {}
        new_ba = {}
        for a in self._a_to_b:
            new_ab[a] = [i.copy() for i in self._a_to_b[a]]

        for a in self._b_to_a:
            new_ba[a] = [i.copy() for i in self._b_to_a[a]]

        g = self.graph1
        if g is not None:
            g = self.graph1.copy()

        return Translation(new_ab, new_ba, graph=g)

    def __repr__(self):
        return "reverse translations: \n  " + str(self._b_to_a) +\
               "\nforward translation: \n   " + str(self._a_to_b)


    __str__ = __repr__

    def _translations(self, rp, inverse=False):
        dict = self._b_to_a if inverse else self._a_to_b
        if rp in dict:
            intervals = dict[rp]
            return intervals
        # Not in dict, means translate to itself (return interval covering
        # itself)
        if self.graph1 is not None and self.graph2 is not None:
            g = self.graph1 if inverse else self.graph2
        else:
            g = self.graph1 if self.graph2 is None else self.graph2
        if g is None:
            length = self.block_lengths[rp]
        else:
            length = g.blocks[rp].length()
        i = [Interval(0, length, [rp], g)]
        i[0].set_length_cache(length)
        return i

    def _get_other_graph(self, inverse):
        if inverse:
            return self.graph1
        else:
            return self.graph2

    def _assert_is_valid(self):
        rps = self._get_range()
        for rp in rps:
            assert rp in self._b_to_a, "%s not in %s" % (rp, self._b_to_a.keys())

    def _copy_blocks(self, subgraph):
        new_blocks = subgraph.blocks.copy()
        for b in self._a_to_b:
            if b in new_blocks:
                del new_blocks[b]
        return new_blocks

    def _translate_subgraph_blocks(self, subgraph, edge_list_add):

        new_blocks = self._copy_blocks(subgraph)

        # Add blocks sthat should be translated
        for block in self._a_to_b:

            if block not in subgraph.blocks:
                continue

            translated = self._translations(block, inverse=False)
            assert len(translated) <= 1, \
                "Only translations to max 1 interval supported. %d returned" \
                % (len(translated))
            translated = translated[0]
            for rp in translated.region_paths:
                if rp not in new_blocks:
                    new_blocks[rp] = self.block_cls(
                        self._translations(rp, inverse=True)[0].length())

            edge_list_add.extend(translated.get_adj_list())  # Add these later

        return new_blocks, edge_list_add

    def _find_rps_to_add_within_translated_interval(
            self, region_path, intervalt, interval, inverse):
        # Find out which of these rps should be added
        # Do not add region paths before the original interval
        intervalt_offset = intervalt.start_position.offset
        offset = 0
        new_region_paths = list(intervalt.region_paths[1:-1])
        # Check only first
        for rp in intervalt.region_paths:
            if inverse:
                assert intervalt.graph is not None, intervalt
                assert rp in intervalt.graph.blocks, \
                    "region path %s in interval %s does not exist in graph %s" % (
                        rp, intervalt, intervalt.graph)
                length = intervalt.graph.blocks[rp].length()
            else:
                length = self._translations(rp, True)[0].length()

            # If first region path
            is_first = interval.region_paths[0] == region_path
            is_tmp = offset + length - intervalt_offset <= interval.start_position.offset
            if is_first and is_tmp:
                offset += length
                continue

            # If last region path
            is_last = interval.region_paths[-1] == region_path
            if is_last and offset - intervalt_offset >= interval.end_position.offset:
                offset += length
                continue
            new_region_paths.append(rp)

            offset += length

        return new_region_paths

    def _translate_interval_nontrivial(self, interval, inverse=False):
        # New faster version: Get all rps, but remove the ones before start and after end rp
        new_starts = self.translate_position(interval.start_position, inverse)
        # Hack: Convert to inclusive end coordinate
        interval.end_position.offset -= 1
        new_ends = self.translate_position(interval.end_position, inverse)
        # Hack: Convert back to exclusive
        interval.end_position.offset += 1
        for new_end in new_ends:
            if new_end != interval.end_position:
                new_end.offset += 1

        is_simple = len(new_ends) == 1 and len(new_starts) == 1

        if len(new_ends) == len(new_starts):
            new_region_paths = [[] for i in range(0, len(new_starts))]
            for region_path in interval.region_paths:
                intervals = self._translations(region_path, inverse)
                for i in range(0, len(intervals)):
                    new_region_paths[i].extend(intervals[i].region_paths)

            # Remove all rps before start rp and after end rp
            new_region_paths_sliced = []
            for i in range(0, len(new_starts)):
                start_rp = new_starts[i].region_path_id
                end_rp = new_ends[i].region_path_id
                start_index = new_region_paths[i].index(start_rp)
                end_index = new_region_paths[i].index(end_rp)
                new_region_paths_sliced.extend(
                    new_region_paths[i][start_index:end_index+1])

            new_region_paths = new_region_paths_sliced
        else:
            new_region_paths = []

            for region_path in interval.region_paths:
                # Find all region paths that this region path follows
                intervals = self._translations(region_path, inverse)
                is_simple = is_simple and len(intervals) == 1
                for intervalt in intervals:
                    # Find out which of these rps should be added
                    # Do not add region paths before the original interval
                    new_region_paths.extend(
                        self._find_rps_to_add_within_translated_interval(
                            region_path, intervalt, interval, inverse))

        if is_simple:
            return SingleMultiPathInterval(
                Interval(new_starts[0], new_ends[0], new_region_paths,
                         graph=self._get_other_graph(inverse),
                         direction=interval.direction))

        return GeneralMultiPathInterval(
                new_starts,
                new_ends,
                new_region_paths,
                self._get_other_graph(inverse))

    def _get_range(self):
        rps = []
        for intervals in self._a_to_b.values():
            for interval in intervals:
                rps.extend(interval.region_paths)
        return rps
