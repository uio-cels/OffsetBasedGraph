from collections import deque, defaultdict
import numpy as np


class NodeDistances(object):
    def __init__(self, node_ids):
        self.distances = {node_id: defaultdict(int) for node_id in node_ids}
        self.distances.update({-node_id: defaultdict(int)
                               for node_id in node_ids})


class DistanceIndex(object):
    def __init__(self, graph, max_distance):
        self.graph = graph
        self.max_distance = max_distance
        self.adj_list = self._create_undirected_graph()

    def _create_undirected_graph(self):
        adj_list = defaultdict(list)
        for node_id, adjs in self.graph.adj_list.items():
            adj_list[node_id].extend(adjs)
        for node_id, adjs in self.graph.reverse_adj_list.items():
            adj_list[node_id].extend(adjs)
        for node_id, adjs in adj_list.items():
            adj_list[node_id] = list(set(adjs))
        return adj_list

    def create(self):
        self.distances = {}
        for node_id in self.graph.blocks.keys():
            self._find_close_to_node(node_id)
            self._find_close_to_node(-node_id)
        self._convert_distance_dicts()

    def _find_close_to_node(self, start_node_id):
        node_size = self.graph.node_size(start_node_id)
        stack = deque([(node_id, node_size) for node_id
                       in self.adj_list[start_node_id]])
        distances = {start_node_id: 0, -start_node_id: 0}
        while stack:
            node_id, distance = stack.pop()
            if distance > self.max_distance:
                continue
            if node_id in distances:
                if distance >= distances[node_id]:
                    continue
            distances[node_id] = distance
            next_nodes = self.adj_list[node_id]
            new_dist = distance + self.graph.node_size(node_id)
            for next_node in next_nodes:
                stack.append((next_node, new_dist))
            next_nodes = self.adj_list[-node_id]
            for next_node in next_nodes:
                stack.append((next_node, distance))
        self.distances[start_node_id] = distances

    def _convert_distance_dicts(self):
        self.covered_neighbours = {}
        self.partial_neighbours = {}
        for node_id, distances in self.distances.items():
            if node_id > 0:
                covered = [
                    abs(other_id) for other_id, d in list(distances.items())+list(self.distances[-node_id].items())
                    if d+self.graph.node_size(other_id) <= self.max_distance]
                self.covered_neighbours[node_id] = list(sorted(set(covered)))
        for node_id, distances in self.distances.items():
            partial = [(other_id, self.max_distance-d)
                       for other_id, d in distances.items()
                       if abs(other_id) not in self.covered_neighbours[abs(node_id)]]

            self.partial_neighbours[node_id] = list(
                sorted(partial, key=lambda x: x[0]))
