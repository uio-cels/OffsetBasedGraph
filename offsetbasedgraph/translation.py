from .util import takes
from .interval import Interval, Position
from .multipathinterval import GeneralMultiPathInterval, SingleMultiPathInterval, SimpleMultipathInterval
import pickle


class Translation(object):

    def __init__(self, translate_dict={}, reverse_dict={}, graph=None):
        """
        Init the translation object with two dicts. Each dict has
        region path IDs as keys and a list of intervals as values.
        :param translate_dict: Mapping of region paths in graph1 to intervals
        in graph2
        :param reverse_dict: Mapping of region paths in graph2 to intervals
        in graph 1
        :return:
        """
        self._a_to_b = translate_dict
        self._b_to_a = reverse_dict
        for k, intervals in translate_dict.items():
            for interval in intervals:
                assert interval.start_position != interval.end_position, "Empty interval in translate: (%s: %s)" % ("T", translate_dict)

        for k, intervals in reverse_dict.items():
            for interval in intervals:
                assert interval.start_position != interval.end_position, "Empty interval in translate: (%s: %s)" % ("R", reverse_dict)


        # Store graph references
        self.graph1 = None
        self.graph2 = None
        if len(self._a_to_b) > 0:
            self.graph2 = list(self._a_to_b.values())[0][0].graph
        if len(self._b_to_a) > 0:
            self.graph1 = list(self._b_to_a.values())[0][0].graph

        # Why this? Removed temporarily
        """
        if self.graph1 is None:
            self.graph1 = self.graph2
        if self.graph2 is None:
            self.graph2 = self.graph1
        """
        if graph is not None:
            self.graph1 = graph

        #if self.graph2 is not None:
            #print("Created translation having graph2")
            #print(self._a_to_b)
        #else:
            #print("Graph2 is none")

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

        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def from_file(file_name):

        with open("%s" % file_name, "rb") as f:
            return pickle.loads(f.read())

    def _translations(self, rp, inverse=False):
        dict = self._b_to_a if inverse else self._a_to_b
        if rp in dict:
            intervals = dict[rp]
            #print("in dict, returning %s" % (intervals))
            return intervals
        # Not in dict, means translate to itself (return interval covering
        # itself)
        # g = self.graph1 if self.graph1 is not None else self.graph2
        if self.graph1 is not None and self.graph2 is not None:
            g = self.graph1 if inverse else self.graph2
        else:
            g = self.graph1 if self.graph2 is None else self.graph2

        assert(g is not None)
        if rp not in g.blocks:
            print(inverse)
            if g == self.graph1:
                print("Chose graph 1")
            elif g == self.graph2:
                print("Chose graph 2")

            raise Exception("Region path %d is not in graph" % rp)

        length = g.blocks[rp].length()
        return [Interval(0, length, [rp], g)]

    def _get_other_graph(self, inverse):
        if inverse:
            return self.graph1
        else:
            return self.graph2

    def translate_rp(self, rp, inverse=False):
        """Translate region path, by dict lookup

        :param rp: region path id
        :returns: translated interval
        :rtype: Interval

        """
        if rp in self._a_to_b:
            return self._a_to_b[rp]

        rp_interval = Interval(Position(rp, 0),
                               Position(rp, self.graph1.blocks[rp].length()),
                               graph=self.graph1)
        return [rp_interval]

    def translate(self, obj):
        """Check type of obj and translate accordingly
        Only for use where on one-to-one translation

        :param obj: Interval/Position
        :returns: translated object
        :rtype: Interval/Position

        """
        if isinstance(obj, Interval):
            ret = self.translate_interval(obj).get_single_path_intervals()
        elif isinstance(obj, Position):
            ret = self.translate_position(obj)
        else:
            raise Exception("Cannot translate %s" % obj)
        assert len(ret) == 1
        return ret[0]

    def _assert_is_valid(self):
        rps = self._get_range()
        for rp in rps:
            assert rp in self._b_to_a, "%s not in %s" % (rp, self._b_to_a)

    def translate_subgraph(self, subgraph):
        """
        Translates a graph (forward). The graph has to be a subgraph of
        graph1 in the translation object.
        :param subgraph: Subgraph to translate.
        :return: Returns the translated subgraph
        """
        assert self.graph1 is not None, "Cannot translate subgraph if graph1 is None"
        set(self.graph1.blocks.keys())


        # Check that subgraph really is a subgraph
        for block in subgraph.blocks:
            assert block in self.graph1.blocks, "%s not in %s" % (block, self.graph1.blocks)

        for adj in subgraph.adj_list:
            assert adj in self.graph1.adj_list

        new_adj = {}
        new_blocks = {}

        from .graph import Graph, Block

        """
        Algorithm: For every interval crossing edge in subgraph, translate.
        Add every region path in translated interval and every edge
        Also, translate every region path and add edges and region paths found
        """
        edge_list_add = []

        # Translate all edges
        for edge in subgraph.get_edges_as_list():
            block1, block2 = edge
            interval = Interval(0, subgraph.blocks[block2].length(),
                                [block1, block2], subgraph)
            #print("Translating interval: " + str(interval))
            translated = self.translate_interval(interval).get_single_path_intervals()
            assert len(translated) <= 1, \
                "Only translations to one interval supported. %d returned" \
                % len(translated)
            translated = translated[0]
            edge_list_add.extend(translated.get_adj_list())  # Add these later

        # Translate all blocks
        for block in subgraph.blocks:
            translated = self._translations(block, inverse=False)
            assert len(translated) <= 1, \
                "Only translations to max 1 interval supported. %d returned" \
                % (len(translated))
            translated = translated[0]
            for rp in translated.region_paths:
                new_blocks[rp] = Block(self._translations(rp, inverse=True)[0].length())
                # translated.graph.blocks[rp].length())
            edge_list_add.extend(translated.get_adj_list())  # Add these later

        # Add all edges we have found
        for edge in edge_list_add:
            if edge[0] in new_adj:
                if edge[1] not in new_adj[edge[0]]:
                    new_adj[edge[0]].append(edge[1])
            else:
                new_adj[edge[0]] = [edge[1]]

        return Graph(new_blocks, new_adj)

    # @takes(Interval)
    def translate_interval(self, interval, inverse=False):
        """
        Translate an interval between the two coordinate systems.
        :param interval:
        :param inverse:
        :return: Returns an interval. If inverse is True, a list of intervals
        is returned.
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
            return SingleMultiPathInterval(interval, self._get_other_graph(inverse))

        if (interval.length() == 0):
            start_poses = self.translate_position(interval.start_position, inverse)
            return SimpleMultipathInterval([Interval(sp, sp, [sp.region_path_id], graph=self._get_other_graph(inverse))
                                            for sp in start_poses])

        new_starts = self.translate_position(interval.start_position, inverse)
        # Hack: Convert to inclusive end coordinate
        interval.end_position.offset -= 1
        new_ends = self.translate_position(interval.end_position, inverse)
        # Hack: Convert back to exclusive
        interval.end_position.offset += 1
        for new_end in new_ends:
            if new_end != interval.end_position:
                new_end.offset += 1

        new_region_paths = []
        is_simple = len(new_ends) == 1 and len(new_starts) == 1
        for region_path in interval.region_paths:
            # Find all region paths that this region path follows
            intervals = self._translations(region_path, inverse)
            is_simple = is_simple and len(intervals) == 1
            for intervalt in intervals:
                #print("Adding " + str(intervalt.region_paths))

                # Find out which of these rps should be added
                # Do not add region paths before the original interval
                intervalt_offset = intervalt.start_position.offset
                offset = 0
                for rp in intervalt.region_paths:
                    if inverse:
                        assert intervalt.graph is not None, "interval %s has graph None" % (intervalt)
                        assert rp in intervalt.graph.blocks, \
                                "region path %s in interval %s does not exist in graph %s" % (rp, intervalt, intervalt.graph)
                        length = intervalt.graph.blocks[rp].length()
                    else:
                        length = self._translations(rp, True)[0].length()

                    # If first region path
                    if interval.region_paths[0] == region_path and offset + length - intervalt_offset <= interval.start_position.offset:
                        #print("Cont on first")
                        #print(offset + length - intervalt_offset)
                        #print(length)
                        offset += length
                        continue
                    # If last region path
                    elif interval.region_paths[-1] == region_path and offset >= interval.end_position.offset:
                        offset += length
                        continue
                    new_region_paths.append(rp)

                    offset += length

                #new_region_paths.extend(intervalt.region_paths)

        if is_simple:
            #print("Is simple. New region paths: ")
            #print(new_region_paths)
            new_interval = SingleMultiPathInterval(
                Interval(new_starts[0], new_ends[0], new_region_paths,
                         graph=self._get_other_graph(inverse)))

        else:
            new_interval = GeneralMultiPathInterval(
                new_starts, new_ends,
                new_region_paths, self._get_other_graph(inverse))

        return new_interval

    # @takes(Position)
    def translate_position(self, position, inverse=False):
        """
        Translates a position
        :param position:
        :param inverse: If True, translate back
        :return: Returns the translated Position
        """
        # Why assert graph here?
        # assert self.graph1 is not None, "Graph1 is none"
        # assert self.graph2 is not None, "Graph2 is none"

        trans_dict = self._b_to_a if inverse else self._a_to_b
        if position.region_path_id not in trans_dict:
            return [position]

        # Get interval for region path. Select first region path. Count offset.
        intervals = self._translations(position.region_path_id, inverse)
        positions = []
        #print("Translating position %s, %d, region path: %d" % (str(position), inverse, position.region_path_id))
        for interval in intervals:
            if True or (not inverse and self.graph2 is None) \
                    or (inverse and self.graph1 is None):
                #print("Finding rp lens for interval %s" % interval)
                rp_lens = [self._translations(rp, inverse=not inverse)[0].length()
                           for rp in interval.region_paths]
            else:
                rp_lens = [interval.graph.blocks[rp].length() for rp in interval.region_paths] #[self._get_other_graph(inverse).blocks[rp].length() for rp in interval.region_paths]

            found_pos = interval.get_position_from_offset(
                position.offset, rp_lens)
            positions.append(found_pos)

        return positions

    def _get_range(self):
        rps = []
        for intervals in self._a_to_b.values():
            for interval in intervals:
                rps.extend(interval.region_paths)
        return rps

    def __add__(self, other):
        """
        Combine (add) two translations.
        :param other: Another translation
        :return:
        """

        """
        Algorithm
            FOr every (region path => inter) in a, translate inter to c
            by using other.translate_interval(inter). Replace in a.
            For every (region path => inter) in c, translate recursviely back
            to all intervals using a.translate(inter, inverse). Replace in c
        """
        assert other.graph1 is not None
        assert self.graph1 is not None

        # Maybe not necessary, but makes things simpler to require this:
        assert self.graph2 is not None
        # self.graph2 should be the same as other.graph1
        assert self.graph2.blocks == other.graph1.blocks, \
                """The translation added does not have graph1
                identical to first translation's graph
                %s\n!=\n%s""" % (self.graph2, other.graph1)

        # Assert intervals have correct graphs
        for intervals in self._a_to_b.values():
            for interval in intervals:
                assert interval.graph.blocks == self.graph2.blocks, \
                    "first translation object has interval %s with graph different than graph2" % (
                        interval)

        for intervals in self._b_to_a.values():
            for interval in intervals:
                assert interval.graph == self.graph1, \
                    "interval %s has graph %s \n!=\n %s \nTranslation: %s" % (
                        interval, interval.graph, self.graph1, self)

        for intervals in other._b_to_a.values():
            for interval in intervals:
                assert interval.graph == other.graph1, \
                    "%s \n != \n %s" % (interval.graph, other.graph1)
                assert interval.graph.blocks == self.graph2.blocks

        new_trans = Translation(self._a_to_b, self._b_to_a, graph=self.graph1)
        region_paths = set(self._a_to_b.keys()).union(set(other._a_to_b.keys()))
        valid_region_paths = region_paths.intersection(set(self.graph1.blocks))

        # Find new forward translation dict
        new_translate_dict = {}
        for t in valid_region_paths:
            translated = other.translate_interval(
                self._translations(t)[0])
            new_translate_dict[t] = [i.copy() for i in translated.get_single_path_intervals()]
        new_trans._a_to_b = new_translate_dict
        i_region_paths = set(self._b_to_a.keys()).union(set(other._b_to_a.keys()))        
        valid_i_region_paths = [rp for rp in i_region_paths if rp not in other._a_to_b]

        new_translate_dict = {}
        for t in valid_i_region_paths:
            translated = []
            for inter in other._translations(t, inverse=True):
                for translated_inter in self.translate_interval(inter, inverse=True).get_single_path_intervals():
                    translated.append(translated_inter.copy())

                """
                translated.extend(
                    self.translate_interval(inter, inverse=True).get_single_path_intervals())
                """
            new_translate_dict[t] = translated
        new_b_to_a = new_translate_dict
        # new_trans._a_to_b = new_translate_dict
 #
 #
 #        changed = True
 #
 #        # valid_i_region_paths = i_region_paths.intersection
 #        while(changed):
 #            changed = False
 #            for t in i_region_paths:  # other._b_to_a:
 #                print("t ", t)
 #                # For every block=>interval, map backwards to intervals
 #                new_intervals = []
 #                for inter in other._translations(t, inverse=True):
 #                    new_intervals.extend(
 #                        #other.translate_interval(inter, True).get_single_path_intervals()
 #                        self.translate_interval(inter, True).get_single_path_intervals()
 #                    )
 #                old = new_b_to_a[t]
 #                if all(i in old for i in new_intervals) \
 #                   and all(i in new_intervals for i in old):
 #                    changed = False  # If change, translate again
 #                print("new intervals")
 #                print(new_intervals)
 #                new_b_to_a[t] = new_intervals
 #
        # Set b_to_a graph to graph1
        copied_graph = self.graph1.copy()
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

    def __eq__(self, other):
        """Check if two translations are equal"""
        eq = self._a_to_b == other._a_to_b
        eq &= self._b_to_a == other._b_to_a
        return eq

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "b to a: \n  " + str(self._b_to_a) + "\na to b: \n   " + str(self._a_to_b)
        return "\n".join(["%s: %s" % (_id, ",".join(intervals))
                          for _id, intervals in
                          self._a_to_b.items()])
