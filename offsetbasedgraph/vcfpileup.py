import numpy as np
from collections import namedtuple, deque
from offsetbasedgraph.vcfgraph import Path
Area = namedtuple("Area", ["complete_nodes", "partial_nodes", "offsets"])


class DistanceOracle:
    def __init__(self, graph):
        self._graph = graph

    def get_shortest_path(self):
        node = 0
        path = [node]
        while len(self.graph._adj_list[node]):
            node = max(self.graph._adj_list[node])
            path.append(node)
        return Path(path, self._graph)

    def get_area_within_distance(self, position, d):
        node_id = position.node_id
        d = self._graph._node_lens(position.node_id)
        queue = deque([])


class VCFPileup:
    def __init__(self, graph):
        self._graph = graph
        self._snp_counter = np.zeros_like(graph._snps._snps)
        self._in_counter = np.zeros_like(graph._node_lens)
        self._out_counter = np.zeros_like(graph._node_lens)
        self._starts = []
        self._ends = []

    def _get_pileup_offset(position):
        pass

    def add_interval(self, interval):
        self._in_counter[interval.node_ids[1:]] += 1
        self._out_counter[interval.node_ids[:-1]] += 1
        self._starts.append(self._get_pileup_offset(interval.start))
        self._ends.append(self._get_pileup_offset(interval.end))

    def add_extension(self, position, n):
        area = self._oracle.get_nodes_within_distance(position, n)
        self._in_counter[area.complete_nodes] += 1
        self._in_counter[area.partial_nodes] += 1
        self._out_counter[area.complete_nodes] += 1
        ends = [self._get_pileup_offset(pos)
                for pos in zip(area.partial_nodes, area.offsets)]
        self._ends.extend(ends)
