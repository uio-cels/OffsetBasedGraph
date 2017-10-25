from collections import defaultdict
import operator
from .interval import Interval

class GeneralArea(object):
    def __init__(self, complete_region_paths, touched_areas):
        self.complete_region_paths = complete_region_paths
        self.touched_areas = touched_areas


def update_areas(areas, new_areas):
    for node_id, intervals in new_areas.items():
        if node_id not in areas:
            areas[node_id] = intervals
            continue
        old_intervals = areas[node_id]
        tuples = [(idx, (-1)**i) for i, idx in
                  enumerate(old_intervals + intervals)]
        tuples.sort(key=lambda x: 3*x[0]-x[1])
        cum_code = 0
        filtered = []
        for idx, code in tuples:
            if cum_code == 0:
                filtered.append(idx)
            cum_code += code
            if cum_code == 0:
                filtered.append(idx)
        areas[node_id] = filtered


class BaseGraphTraverser(object):
    def __init__(self, graph, direction=+1):
        self.graph = graph
        self.node_sizes = graph._node_sizes
        self.adj_list = graph.adj_list
        if direction < 0:
            self.adj_list = graph.reverse_adj_list


    def stop_traversing(self, current_node):
        raise NotImplementedError()


class GraphTraverserUsingSequence(BaseGraphTraverser):
    def __init__(self, graph, search_sequence, graph_sequence_retriever, direction=+1):
        super(GraphTraverserUsingSequence, self).__init__(graph, direction)
        self.search_sequence = search_sequence.replace("\n", "")
        self.graph_sequence_retriever = graph_sequence_retriever
        self.path_found = None
        self.valid_nodes_in_path = {}
        self.traversed_nodes = {}
        self.visited = {}

        self.end_nodes = {}

    def _stop_recursion(self, node_id, offset):
        node_size = self.graph.node_size(node_id)
        if self.graph_sequence_retriever.get_sequence_on_directed_node(node_id) \
                != self.search_sequence[offset:offset+node_size]:
            return True
        return False

    def search_from_node(self, node_id):
        self.extend_from_block(node_id)
        #self.valid_nodes_in_path[node_id] = True

    def _get_paths_by_backtracking(self):
        print("Backtracking")
        print("N end nodes: %d" % len(self.end_nodes))

        reversed_sequence = self.search_sequence[::-1]

        sorted_keys = sorted(self.end_nodes.items(), key=operator.itemgetter(1))
        print("Sorted keys")
        #print(sorted_keys)
        end_nodes_in_longest_path = sorted_keys[-1][0]

        nodes = []
        current_node = end_nodes_in_longest_path
        print("Backtracking from %d" % current_node)
        offset = 0
        while True:
            #print("  %d" % current_node)
            nodes.append(current_node)
            node_size = self.graph.node_size(current_node)
            print(node_size)
            correct_seq = reversed_sequence[offset:offset+node_size]
            graph_seq = self.graph_sequence_retriever.get_sequence_on_directed_node(-current_node)
            print(correct_seq)
            print(graph_seq)
            print(len(correct_seq))
            print("Last char: -%s-" % correct_seq[-1])
            print("First char: -%s-" % correct_seq[0])
            print(len(graph_seq))
            assert graph_seq == correct_seq

            print("Current node: %d" % current_node)
            next = None
            prev_candidates = []
            for previous in self.graph.reverse_adj_list[-current_node]:
                #print("   prev: %d" % previous)
                if -previous in self.traversed_nodes:
                    print("Checking candidate %d" % -previous)
                    prev_node_size = self.graph.node_size(previous)
                    correct_seq = reversed_sequence[offset+node_size:offset+node_size+prev_node_size]
                    seq_on_prev = self.graph_sequence_retriever.get_sequence_on_directed_node(previous)
                    print(correct_seq)
                    print(seq_on_prev)
                    if correct_seq == seq_on_prev:
                        prev_candidates.append(-previous)

            offset += node_size
            print("Len prev candidates")
            print(len(prev_candidates))
            assert len(prev_candidates) <= 1

            if len(prev_candidates) == 0:
                break
            next = prev_candidates[0]

            current_node = next

        return nodes[::-1]

    def get_nodes_found(self):
        return self.path
        return self._get_paths_by_backtracking()

    def get_interval_found(self):
        nodes = self._get_paths_by_backtracking()
        end_pos = self.graph.node_size(nodes[-1])
        return Interval(0, end_pos, nodes, self.graph)

    def extend_from_block(self, start_node_id):
        path = []
        stack = [(start_node_id, 0 , 0, [])]
        while stack:
            node_id, offset, prev_node, path = stack.pop(0)
            print("Checking node %d, offset %d, path: %s" % (node_id,  offset, path))
            node_size = self.graph.node_size(node_id)
            path = path.copy()
            path.append(node_id)
            if self._stop_recursion(node_id, offset):
                #print("Stopping at %d" % node_id)
                self.end_nodes[prev_node] = offset
                if offset < len(self.search_sequence):
                    print("Stopping at invalid end, offset %d" % offset)
                    if prev_node in self.traversed_nodes:
                        del self.traversed_nodes[prev_node]
            else:
                #print("Offset + node size: %d, len sequence: %d" % (node_size + offset, len(self.search_sequence)))
                self.traversed_nodes[node_id] = True

                stack.extend([(next_id, offset + node_size, node_id, path) for
                              next_id in self.adj_list[node_id]])

                if len(self.adj_list[node_id]) == 0:
                    #print("Found end node at %d, offset %d" % (node_id, offset + node_size))
                    self.end_nodes[node_id] = offset + node_size

                if node_size + offset == len(self.search_sequence):
                    print("Found end at node %d, offset %d" % (node_id, offset + node_size))
                    self.path = path
                    return

class GraphTraverser(object):
    def __init__(self, graph, direction=+1):
        self.graph = graph
        self.node_sizes = graph._node_sizes
        self.adj_list = graph.adj_list
        if direction < 0:
            self.adj_list = graph.reverse_adj_list

    def extend_from_block(self, start_node_id, start_length, visited):
        stack = [(start_node_id, start_length)]
        while stack:
            node_id, length = stack.pop(0)
            if visited[node_id] >= length:
                continue
            visited[node_id] = length
            remaining = length-self.node_sizes[abs(node_id)]
            # ]graph.node_size(node_id)
            stack.extend([(next_id, remaining) for
                          next_id in self.adj_list[node_id]])

    def extend_from_block_rec(self, node_id, length, visited):
        if node_id in visited and visited[node_id] >= length:
            return
        visited[node_id] = length
        if length <= self.graph.node_size(node_id):
            return
        remaining = length-self.graph.node_size(node_id)
        for next_id in self.adj_list[node_id]:
            self.extend_from_block(next_id, remaining, visited)
