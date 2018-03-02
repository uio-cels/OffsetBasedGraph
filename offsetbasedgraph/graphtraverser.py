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
        self.adj_list = graph.adj_list
        if direction < 0:
            self.adj_list = graph.reverse_adj_list


    def stop_traversing(self, current_node):
        raise NotImplementedError()


class GraphTravserserBetweenNodes(BaseGraphTraverser):

    def __init__(self, graph):
        self.graph = graph
        self._visited = {}

    def get_snarl_subgraph(self, start, end, include_start_and_end=False, print_debug=False):
        blocks = self._get_blocks_in_snarl_subgraph(start, end, print_debug=print_debug)

        if include_start_and_end:
            blocks.append(start)
            blocks.append(end)

        #print("Blocks: %s" % blocks)

        blocks = list(set(blocks))
        block_dict = {}
        for block in blocks:
            block_dict[block] = self.graph.blocks[block]

        edges = defaultdict(list)
        for block in blocks:
            for edge in self.graph.adj_list[block]:
                if edge in block_dict:
                    edges[block].append(edge)


        return Graph(block_dict, edges, do_not_create_node_indices=True)

    def _get_blocks_in_snarl_subgraph(self, start, end, include_start_and_end=True, print_debug=False):
        # Subgraph is connnected to graph by only two other nodes (start and end). Start can have other edges out,
        # than those into subgraph
        blocks = []
        visited = {}
        visited_false_paths = {}
        stack = [(start, 0, 0, [])]
        i = 0

        did_find_end_node = False

        while stack:
            node_id, n_nodes_to_here, prev_node, path = stack.pop()
            new_path = path.copy()
            if print_debug:
                print("%d: Current node: %d. End: %d" % (i, node_id, end))

            i += 1

            if node_id == end:
                blocks.extend(new_path)
                did_find_end_node = True
                for n in path:
                    visited[n] = True
                if print_debug:
                    print("  Reached end")
                continue

            if node_id in visited:
                blocks.extend(new_path)
                if print_debug:
                    print("Already visited  ")
                continue

            if node_id != start:
                new_path.append(node_id)

            #visited[node_id] = True

            nexts = self.graph.adj_list[node_id]
            #if print_debug:
            #    print("Nexts: %s " % nexts)
            if len(nexts) == 0:
                if print_debug:
                    print("   Reached end, do not add path")
                continue  # Do not add path

            stack.extend([(next_id, n_nodes_to_here + 1, node_id, new_path) for next_id in nexts])

        assert did_find_end_node, "Did not find end node between %d and %d.\n" % (start, end)

        return list(set(blocks))



    def get_greedy_subgraph_between_nodes(self, start, end, include_start_and_end=True, print_debug=False):
        # Faster: Assumes everything after start node is the subgraph
        blocks = {}
        edges = defaultdict(list)

        stack = [(start, None)]
        i = 0
        while stack:
            current_node, prev = stack.pop()
            if print_debug:
                print("Current node: %d, prev: %s" % (current_node, prev))

            if include_start_and_end or (current_node != start and current_node != end):
                blocks[current_node] = self.graph.blocks[current_node]

            if prev is not None:
                if include_start_and_end or (prev != start and current_node != end):
                    edges[prev].append(current_node)

            if current_node in self._visited:
                #print("  Already visited %d" % current_node)
                continue

            self._visited[current_node] = True
            if current_node == end:
                continue

            nexts = self.graph.adj_list[current_node]
            if print_debug:
                print("  Nexts: %s" % nexts)

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
        #print("Linear seq: %s" % correct_seq)
        graph_seq = self.graph_sequence_retriever.get_sequence_on_directed_node(node_id)
        #print("Graph seq : %s" % graph_seq)

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
            node_id, offset, prev_node, n_nodes, path_to_here = stack.pop()
            #print("Checking node %d, offset %d. Added from %d" % (node_id,  offset, prev_node))


            new_path_to_here = path_to_here.copy()
            new_path_to_here.append(node_id)
            print("    Path: %s" % new_path_to_here)
            node_size = self.graph.node_size(node_id)
            n_delete = len(path)-n_nodes
            if n_delete > 0:
                del path[-n_delete:]

            path.append(node_id)

            if self._stop_recursion(node_id, offset):
                self.end_nodes[prev_node] = offset
            else:
                nexts = self.adj_list[node_id]
                stack.extend([(next_id, offset + node_size, node_id, n_nodes + 1, new_path_to_here) for
                              next_id in nexts])

                if node_size + offset == len(self.search_sequence):
                    self.path = path
                    return


class GraphTraverser(object):
    def __init__(self, graph, direction=+1):
        self.graph = graph
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
            remaining = length-self.graph.node_size(node_id)
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


class LengthGraphTraverser(object):
    """
    Finds shortest distance to all nodes within given distance away
    """
    def __init__(self, graph, distance, direction=+1):
        self.graph = graph
        self.adj_list = graph.adj_list
        self.distance = distance
        if direction < 0:
            self.adj_list = graph.reverse_adj_list

    def extend_from_block(self, node_id):
        visited = {}
        stack = [(node_id, 0)]
        while stack:
            node_id, distance_to_node = stack.pop(0)
            #print("  Extending from %d (distance so far: %d)" % (node_id, distance_to_node))
            dist_to_next = distance_to_node + self.graph.node_size(node_id)
            if dist_to_next > self.distance:
                continue

            for next_node in self.adj_list[node_id]:
                if next_node in visited:
                    if visited[next_node] <= dist_to_next:
                        continue
                stack.append((next_node, dist_to_next))
                visited[next_node] = dist_to_next

        return visited