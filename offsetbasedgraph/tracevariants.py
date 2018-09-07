def interval_to_variants_func(reference, graph):
    reference_nodes = reference.nodes_in_interval()
    get_next_node = next_node_func(graph, reference)

    def interval_to_variants(interval):
        if interval.region_paths[0] < 0:
            assert np.all(interval.region_paths < 0)
            interval = interval.get_reverse()
            nodes = interval.region_paths
        variant_nodes = [node-graph.min_node for node in nodes if node not in reference_nodes] # Offset nodes to array idx
        variant_edges = [edge for edge in zip(nodes[:-1], nodes[1:])
                         if all(node in reference for node in edge) and
                         node[1] != get_next_node(node[0])]
        return variant_nodes, variant_edges


def get_variants_from_intervals(reference, graph, intervals):
    interval_to_variants = interval_to_variants_func(reference, graph)
    return (interval_to_variants(interval) for interval in intervals)


def get_ids_from_variants(variant_maps, variants_list):
    return [nodes_edges_to_variant_ids(variants[0], variants[1], variant_maps) for variants in variants_list]
