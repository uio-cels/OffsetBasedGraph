from offsetbasedgraph import Graph, Position


class GeneralArea(object):
    def __init__(self, complete_region_paths, touched_areas):
        self.complete_region_paths = complete_region_paths
        self.touched_areas = touched_areas


def update_areas(areas, new_areas):
    print("#", areas, new_areas)
    for node_id, intervals in new_areas.items():
        if node_id not in areas:
            areas[node_id] = intervals
            continue
        old_intervals = areas[node_id]
        tuples = [(idx, (-1)**i) for i, idx in enumerate(old_intervals + intervals)]
        tuples.sort(key=lambda x: 3*x[0]-x[1])
        print(tuples)
        cum_code = 0
        filtered = []
        for idx, code in tuples:
            if cum_code == 0:
                filtered.append(idx)
            cum_code += code
            if cum_code == 0:
                filtered.append(idx)
        areas[node_id] = filtered
    print("=", areas)


class GraphTraverser(object):
    def __init__(self, graph):
        self.graph = graph

    def shift(self, position, offset):
        node_id = position.region_path_id
        start_offset = position.offset
        cur_node_size = self.graph.node_size(node_id)
        if cur_node_size > start_offset + offset:
            return [Position(node_id, start_offset + offset)]
        new_offset = offset-(cur_node_size-start_offset)
        shifted_positions = []
        for next_node in self.graph.adj_list[node_id]:
            shifted_positions.extend(
                self.shift(Position(next_node, 0),
                           new_offset))

        return shifted_positions

    def guided_shift(self, interval, offset):
        interval_length = interval.length()
        if interval_length > offset:
            return [interval.get_position_from_offset(offset)]

        return self.shift(interval.end_position,
                          offset-interval_length)

    def get_areas_from_point(self, point, length):
        start_node = point.region_path_id
        start_offset = point.offset
        node_length = self.graph.node_size(start_node)
        if node_length > start_offset + length:
            areas = {start_node: [start_offset, start_offset+length]}
            print("##", areas)
            return areas
        all_areas = {start_node: [point.offset, node_length]}
        for next_node in self.graph.adj_list[start_node]:
            areas = self.get_areas_from_point(
                Position(next_node, 0),
                length-(node_length-start_offset))
            update_areas(all_areas, areas)
        return all_areas
