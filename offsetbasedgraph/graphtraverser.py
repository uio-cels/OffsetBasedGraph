from offsetbasedgraph import Graph, Position


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


class GraphTraverser(object):
    def __init__(self, graph, direction=+1):
        self.graph = graph
        self.adj_list = graph.adj_list
        if direction < 0:
            self.adj_list = graph.reverse_adj_list

    def extend_from_block(self, node_id, length, visited):
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
