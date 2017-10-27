from collections import deque, defaultdict
from itertools import chain
import numpy as np
import pickle


class NodeDistances(object):
    def __init__(self, max_distance=None):
        self.distances = {}
        self.max_distance = max_distance

    @classmethod
    def from_dict(cls, distances, max_distance=None):
        obj = cls(max_distance)
        obj.distances = distances
        return obj

    def set_distance(self, node_id, distance):
        if node_id in self.distances and self.distances[node_id] <= distance:
            return False
        self.distances[node_id] = distance
        return True

    def __str__(self):
        return str(self.distances)

    __repr__ = __str__

    def __add__(self, distance):
        return self.from_dict({node_id: prev_distance+distance
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

    def clean(self):
        self.distances = {node_id: distance for node_id, distance in
                          self.items() if distance < self.max_distance}

    def clean_to_numpy(self):
        node_ids = np.array(list(self.distances.keys()), dtype="int")
        distances = np.array([self.distances[node_id] for node_id in node_ids])
        small_distances = distances <= self.max_distance
        return NodeDistancesNP(node_ids[small_distances],
                               distances[small_distances],
                               self.max_distance)


class NodeDistancesNP(object):
    def __init__(self, node_ids, distances, max_distance=None):
        self.node_ids = node_ids
        self.distances = distances
        self.max_distance = max_distance

    def __str__(self):
        return str(self.node_ids) + "\n" + str(self.distances)

    def __len__(self):
        return self.node_ids.size

    __repr__ = __str__

    def __add__(self, distance):
        return NodeDistancesNP(self.node_ids, self.distances + distance,
                               self.max_distance)

    def __eq__(self, other):
        return self.distances == other.distances and self.node_ids == other.node_ids

    def items(self):
        return zip(self.node_ids, self.distances)

    def get_nodes_with_margin(self, node_sizes):
        margin_sizes = self.distances+node_sizes[np.abs(self.node_ids)]
        self.__within_margin = margin_sizes < self.max_distance
        return self.node_ids[self.__within_margin]

    def get_remains(self):
        not_within_margin = np.logical_not(self.__within_margin)  # NOT
        remains = self.max_distance-self.distances[not_within_margin]
        remain_ids = self.node_ids[not_within_margin]
        return self.__class__(remain_ids, remains, self.max_distance)


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
            self._find_close_to_node(node_id)
            self._find_close_to_node(-node_id)
            if i % 1000 == 0:
                print("#", i)
                print(len(self.distances[node_id]))
                print(len(self.distances[-node_id]))

            i += 1

        self._convert_distance_dicts()

    def _find_close_to_node(self, start_node_id):
        node_size = self.graph.node_size(start_node_id)
        stack = deque([(node_id, node_size) for node_id
                       in self.adj_list[start_node_id]])

        distances = NodeDistances(self.max_distance)
                                  
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
            if -node_id in self.distances:
                new_distances = self.distances[-node_id] + (distance - self.graph.node_size(node_id))
                distances.update(new_distances)
            else:
                next_nodes = self.adj_list[-node_id]
                for next_node in next_nodes:
                    stack.append((next_node, distance))
        self.distances[start_node_id] = distances.clean_to_numpy()

    def convert_pos_and_neg_node_distances(self, node_id):
        distances = self.distances[node_id]
        neg_distances = self.distances[-node_id]
        pos_covered = np.abs(distances.get_nodes_with_margin(self.graph._node_sizes))
        neg_covered = np.abs(neg_distances.get_nodes_with_margin(self.graph._node_sizes))
        all_covered = np.unique(np.concatenate((pos_covered, neg_covered)))
        self.covered_neighbours[node_id] = np.sort(all_covered)
        pos_partial = distances.get_remains()
        neg_partial = neg_distances.get_remains()
        self.partial_neighbours[node_id] = np.array([
            (other_id, dist) for other_id, dist in sorted(pos_partial.items(), key=lambda x: x[0]) if
            abs(other_id) not in self.covered_neighbours[abs(node_id)]])
        self.partial_neighbours[-node_id] = np.array([
            (other_id, dist) for other_id, dist in sorted(neg_partial.items(), key=lambda x: x[0]) if
            abs(other_id) not in self.covered_neighbours[abs(node_id)]])
        del self.distances[node_id]
        del self.distances[-node_id]

    def _convert_distance_dicts(self):
        self.covered_neighbours = {}
        self.partial_neighbours = {}

        pos_node_ids = [node_id for node_id in self.distances if node_id > 0]
        for node_id in pos_node_ids:
            self.convert_pos_and_neg_node_distances(node_id)

    def to_file(self, file_name):
        with open("%s" % file_name+".cov", "wb") as f:
            pickle.dump(self.covered_neighbours, f)

        with open("%s" % file_name+".par", "wb") as f:
            pickle.dump(self.partial_neighbours, f)

    @classmethod
    def from_file(cls, graph, max_distance, file_name):
        obj = cls(graph, max_distance)
        with open("%s.cov" % file_name, "rb") as f:
            cov = pickle.loads(f.read())
            obj.covered_neighbours = cov

        with open("%s.par" % file_name, "rb") as f:
            par = pickle.loads(f.read())
            obj.partial_neighbours = par
        return obj
