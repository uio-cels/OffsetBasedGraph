from offsetbasedgraph import *


def get_simple_graph():
    main_interval = LinearInterval("hg38", "main", 0, 50)
    alt_interval = LinearInterval("hg38", "alt", 0, 30)
    region_paths = [
        RegionPath("main", {"main": main_interval}),
        RegionPath("alt", {"alt": alt_interval})]

    graph = OffsetBasedGraph("simplegraph")
    graph.blocks = {rp.id: rp for rp in region_paths}
    return graph

main_intervals = [
    LinearInterval("hg38", "main", 0, 10),
    LinearInterval("hg38", "main", 10, 20),
    LinearInterval("hg38", "main", 20, 30),
    LinearInterval("hg38", "main", 30, 40),
    LinearInterval("hg38", "main", 40, 50)]

alt_intervals = [
    LinearInterval("hg38", "alt", 0, 10),
    LinearInterval("hg38", "alt", 10, 20),
    LinearInterval("hg38", "alt", 20, 30)]

region_paths = [RegionPath("main" + str(i), {"main": li})
                for i, li in enumerate(main_intervals)]
region_paths[1].linear_references["alt"] = alt_intervals[0]
region_paths[3].linear_references["alt"] = alt_intervals[-1]


block_edges = {rp1.id: [rp2.id] for rp1, rp2 in
               zip(region_paths[:-1], region_paths[1:])}

region_paths.append(
    RegionPath("alt", {"alt": alt_intervals[1]}))
block_edges[region_paths[1].id].append("alt")
block_edges["alt"] = region_paths[3].id

blocks = {rp.id: rp for rp in region_paths}
graph = OffsetBasedGraph("test")
graph.blocks = blocks
graph.block_edges = block_edges
graph._create_block_index()
