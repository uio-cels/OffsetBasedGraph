from collections import defaultdict
import operator

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
        self.search_sequence = search_sequence
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
        #self.extend_from_block_recursive(node_id, 0)
        self.extend_from_block(node_id)
        self.valid_nodes_in_path[node_id] = True

    def _get_paths_by_backtracking(self):
        print("Backtracking")
        print("N end nodes: %d" % len(self.end_nodes))

        sorted_keys = sorted(self.end_nodes.items(), key=operator.itemgetter(1))
        end_nodes_in_longest_path = sorted_keys[-1][0]


        nodes = []
        current_node = end_nodes_in_longest_path
        while True:
            print("  %d" % current_node)
            nodes.append(current_node)
            next = None
            for previous in self.graph.reverse_adj_list[-current_node]:
                print("   prev: %d" % previous)
                if -previous in self.traversed_nodes:
                    next = -previous
                    break

            if next is None:
                break

            current_node = next

        return nodes[::-1]

    def get_nodes_found(self):
        return self._get_paths_by_backtracking()
        print("Valid paths:")
        print(self.valid_nodes_in_path)

        nodes = []
        for node, val in self.traversed_nodes.items():
            if val:
                nodes.append(node)

        return nodes

    def extend_from_block(self, start_node_id):
        stack = [(start_node_id, 0 , 0)]
        while stack:
            node_id, offset, prev_node = stack.pop(0)
            print("Checking node %d, offset %d" % (node_id,  offset))
            node_size = self.graph.node_size(node_id)

            """
            if node_id in self.visited:
                print(" already visited")
                continue

            self.visited[node_id] = True
            """

            if self._stop_recursion(node_id, offset):
                print("Stopping at %d" % node_id)
                self.end_nodes[prev_node] = offset
            else:
                print("Offset + node size: %d, len sequence: %d" % (node_size + offset, len(self.search_sequence)))
                self.traversed_nodes[node_id] = True

                stack.extend([(next_id, offset + node_size, node_id) for
                              next_id in self.adj_list[node_id]])
                if node_size + offset == len(self.search_sequence):
                    pass

                if len(self.adj_list[node_id]) == 0:
                    print("Found end node at %d" % node_id)
                    self.end_nodes[node_id] = offset + node_size

    def extend_from_block_recursive(self, node_id, offset):
        print("Checking from %d, offset %d" % (node_id, offset))
        #print("Checking from %d, offset %d, nodes found: %s" % (node_id, offset, nodes_found))
        node_size = self.graph.node_size(node_id)

        if self.path_found:
            print("  Already found a path")
            return

        if self._stop_recursion(node_id, offset):
            print("Stopping at %d" % node_id)
            if node_size + offset == len(self.search_sequence):
                print("Found path")
                self.path_found = True
            return False

        #new_nodes_found = nodes_found.copy()
        #new_nodes_found.append(node_id)

        n_valid_paths = 0
        for next_id in self.adj_list[node_id]:
            if self.extend_from_block_recursive(next_id, offset + node_size):
                n_valid_paths += 1
                self.valid_nodes_in_path[next_id] = True

        if n_valid_paths > 0:
            return True
        else:
            if node_size + offset == len(self.search_sequence):
                self.path_found = True
                print("  Found path, sequence stopping")
            print("  No new valid edges")
            return True

        return False



class GraphTraverser(BaseGraphTraverser):
    def __init__(self, graph, visited, direction=+1):
        super(GraphTraverser, self).__init__(graph, direction)
        self.visited = visited

    def extend_from_block(self, start_node_id, start_length):
        stack = [(start_node_id, start_length)]
        while stack:
            node_id, length = stack.pop(0)
            if self.visited[node_id] >= length:
                continue
            self.visited[node_id] = length
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

#    def shift(self, position, offset):
#        node_id = position.region_path_id
#        start_offset = position.offset
#        cur_node_size = self.graph.node_size(node_id)
#        new_offset = start_offset+offset
#        if cur_node_size > new_offset and new_offset >= 0:
#            return [Position(node_id, new_offset)]
#
#        if new_offset >= 0:
#            new_offset -= cur_node_size
#            shifted_positions = []
#            adj_list = self.graph.adj_list if offset >= 0 else self.graph.reverse_adj_list
#        for next_node in adj_list[node_id]:
#            next_start = 0 if offset >= 0 else self.graph.node_size(next_node)
#            shifted_positions.extend(
#                self.shift(Position(next_node, next_start),
#                           new_offset))
#        return shifted_positions
#
#    def guided_shift(self, interval, offset):
#        interval_length = interval.length()
#        if interval_length > abs(offset):
#            if offset < 0:
#                offset = interval_length-offset
#            return [interval.get_position_from_offset(offset)]
#
#        if offset < 0:
#            return self.shift(interval.start_position,
#                              offset+interval_length)
#
#        return self.shift(interval.end_position,
#                          offset-interval_length)
#
#    def get_areas_from_point(self, point, length):
#        start_node = point.region_path_id
#        start_offset = point.offset
#        node_length = self.graph.node_size(start_node)
#        new_offset = start_offset + length
#        if node_length > new_offset and new_offset >= 0:
#            interval = [start_offset, new_offset]
#            if length < 0:
#                interval = interval[::-1]
#            return {start_node: interval}
#
#        if length < 0:
#            all_areas = {start_node: [0, point.offset]}
#        else:
#            all_areas = {start_node: [point.offset, node_length]}
#        adj_list = self.graph.adj_list if length >= 0 else self.graph.reverse_adj_list
#        if length >= 0:
#            new_offset -= node_length
#        for next_node in adj_list[start_node]:
#            next_offset = 0 if length >= 0 else self.graph.node_size(next_node)
#            areas = self.get_areas_from_point(
#                Position(next_node, next_offset),
#                new_offset)
#            update_areas(all_areas, areas)
#        return all_areas
