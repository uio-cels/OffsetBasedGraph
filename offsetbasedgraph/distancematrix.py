from collections import deque, defaultdict
import numpy as np


class NodeDistances(object):
    def __init__(self, distances=None, max_distance=None):
        self.distances = distances if distances is not None else {}
        self.max_distance = max_distance
        if max_distance is not None:
            self.distances = {node_id: distance for node_id, distance in
                              self.items() if distance < max_distance}

    def set_distance(self, node_id, distance):
        if node_id in self.distances and self.distances[node_id] <= distance:
            return False
        self.distances[node_id] = distance
        return True

    def __str__(self):
        return str(self.distances)

    __repr__ = __str__

    def __add__(self, distance):
        return self.__class__({node_id: prev_distance+distance
                               for node_id, prev_distance in
                               self.distances.items()},
                              self.max_distance)

    def __eq__(self, other):
        return self.distances == other.distances

    def update(self, other):
        for node_id, distance in other.items():
            self.set_distance(node_id, distance)

    def items(self):
        return self.distances.items()


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
        i = 0
        for node_id in self.graph.blocks.keys():
            if i % 10000 == 0:
                print(i)
            i += 1
            self._find_close_to_node(node_id)
            self._find_close_to_node(-node_id)
        self._convert_distance_dicts()

    def _find_close_to_node(self, start_node_id):

        node_size = self.graph.node_size(start_node_id)
        stack = deque([(node_id, node_size) for node_id
                       in self.adj_list[start_node_id]])

        distances = NodeDistances({start_node_id: 0,
                                   -start_node_id: node_size},
                                  self.max_distance)
        while stack:
            node_id, distance = stack.pop()
            if distance > self.max_distance:
                continue
            if not distances.set_distance(node_id, distance):
                continue
            new_dist = distance + self.graph.node_size(node_id)
            if node_id in self.distances:
                new_distances = self.distances[node_id] + distance
                distances.update(new_distances)
            else:
                next_nodes = self.adj_list[node_id]
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
