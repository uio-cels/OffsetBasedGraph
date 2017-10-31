from collections import defaultdict
import operator
from .interval import Interval
from .graph import Graph

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


class GraphTravserserBetweenNodes(BaseGraphTraverser):

    def __init__(self, graph):
        self.graph = graph

    def get_greedy_subgraph_between_nodes(self, start, end, include_start_and_end=True):
        blocks = {}
        edges = defaultdict(list)

        stack = [(start, None)]
        i = 0
        while stack:
            current_node, prev = stack.pop()
            #print("Current node: %d, prev: %s" % (current_node, prev))

            if include_start_and_end or (current_node != start and current_node != end):
                blocks[current_node] = self.graph.blocks[current_node]

            if prev is not None:
                if include_start_and_end or (prev != start and current_node != end):
                    edges[prev].append(current_node)

            if current_node == end:
                continue

            nexts = self.graph.adj_list[current_node]
            stack.extend([(next_id, current_node) for next_id in nexts])


        return Graph(blocks, edges)



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
        correct_seq = self.search_sequence[offset:offset+node_size]
        print("Linear seq: %s" % correct_seq)
        graph_seq = self.graph_sequence_retriever.get_sequence_on_directed_node(node_id)
        print("Graph seq : %s" % graph_seq)

        if graph_seq != correct_seq:
            return True
        return False

    def search_from_node(self, node_id):
        self.extend_from_block(node_id)

    def get_nodes_found(self):
        return self.path

    def get_interval_found(self):
        nodes = self.path
        end_pos = self.graph.node_size(nodes[-1])
        return Interval(0, end_pos, nodes, self.graph)

    def extend_from_block(self, start_node_id):
        path = []
        stack = [(start_node_id, 0, 0, 0, [])]
        i = 0
        while stack:
            #if i % 10000 == 0:
            #    print("Node %i" % i)
            #i += 1
            node_id, offset, prev_node, n_nodes, path_to_here = stack.pop()
            #if node_id in [6512, 6508, 6511] or node_id < 0:
            print("Checking node %d, offset %d. Added from %d" % (node_id,  offset, prev_node))


            new_path_to_here = path_to_here.copy()
            new_path_to_here.append(node_id)
            print("    Path: %s" % new_path_to_here)
            node_size = self.graph.node_size(node_id)
            n_delete = len(path)-n_nodes
            #print("N nodes to here: %d" % n_nodes)
            #print("Deleting from end of paths: %d" % n_delete)
            if n_delete > 0:
                del path[-n_delete:]
            #print("Path after cutting: %s" % path)
            path.append(node_id)

            if self._stop_recursion(node_id, offset):
                print("Stopping at %d" % node_id)
                self.end_nodes[prev_node] = offset
            else:
                nexts = self.adj_list[node_id]
                print("     Next: %s" % nexts)
                stack.extend([(next_id, offset + node_size, node_id, n_nodes + 1, new_path_to_here) for
                              next_id in nexts])

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
