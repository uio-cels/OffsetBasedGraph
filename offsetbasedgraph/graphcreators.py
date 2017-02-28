import os
from .translation import Translation
from .graph import Graph, Block
from .interval import Interval
from offsetbasedgraph.sequences import get_sequence_ucsc
import csv


def merge_flanks(intervals, final_trans, new_graph, name_translation):
    """
    Merges the start and end flanks represented in intervals on the given graph
    :param intervals: List of start and end flanks to merge [start, start, end, end]
    Intervals should be on first graph (i.e. name_translation.graph1)
    :param final_trans: Trans that will be updated
    :param new_graph: Current graph
    :param name_translation: Translation from human readable names to numeric IDs. Can be an empty translation
    :return: Returns the new graph and translation as a tuple
    """
    # Merge start flank of alt locus with main
    merge_intervals = intervals[0:2]
    merge_intervals = [name_translation.translate(i) for i in merge_intervals]
    merge_intervals = [final_trans.translate(i) for i in merge_intervals]

    for intv in merge_intervals:
        intv.graph = new_graph
    if merge_intervals[0].length() > 0:
        new_graph, trans = new_graph.merge(merge_intervals)
        final_trans += trans
    else:
        # Only connect by edge
        new_graph, trans = new_graph.connect_postitions(
            new_graph.prev_position(merge_intervals[0].start_position),
            merge_intervals[1].start_position
        )
        final_trans += trans
        final_trans.graph2 = new_graph.copy()
        final_trans.graph2._update_a_b_graph(final_trans._a_to_b, new_graph)

    # Merge end flank of alt locus with main

    merge_intervals = intervals[2:4]

    if merge_intervals[0].length() > 0:
        merge_intervals = [final_trans.translate(name_translation.translate(i))
                           for i in merge_intervals]
        for intv in merge_intervals:
            intv.graph = new_graph

        new_graph, trans = new_graph.merge(merge_intervals)
        final_trans += trans
    else:
        # Change position 1 back for alt loci
        ig = name_translation.graph1
        merge_intervals[1].start_position = \
            ig.prev_position(merge_intervals[1].start_position)

        merge_intervals = [final_trans.translate(name_translation.translate(i))
                           for i in merge_intervals]

        for intv in merge_intervals:
            intv.graph = new_graph

        # Only connect by edge
        new_graph, trans = new_graph.connect_postitions(
            merge_intervals[1].start_position,
            new_graph.next_position(merge_intervals[0].start_position)
        )

        assert len(new_graph.adj_list[merge_intervals[1].start_position.region_path_id]) > 0

        #print("=== end flank trans ===")
        #print(trans)

        final_trans += trans
        final_trans.graph2 = new_graph.copy()
        final_trans.graph2._update_a_b_graph(final_trans._a_to_b, new_graph)
        #print("No end flank")

    return new_graph, final_trans


def connect_without_flanks(graph, alt_loci_fn, name_translation,
                           filter_alt_loci=[]):
    """
    Connects the alternative loci in the given file to the grch38 graph,
    without flanks.
    :param alt_loci_fn: Filename of file containing alternative loci.
    One alt locus on each line.
    Four columns: alt_locus_id  chr chr_start   chr_stop
    :param filter_alt_loci: If not empty, only these alt loci will be connected
    :return: Returns the new graph
    """
    from .GRCH38 import AltLoci

    alt_loci = AltLoci.from_file(alt_loci_fn, filter_alt_loci)
    new_graph = graph
    final_trans = Translation(graph=graph)
    final_trans.graph2 = graph

    for alt_locus in alt_loci.alt_loci:
        if len(filter_alt_loci) > 0 and alt_locus.name not in filter_alt_loci:
            continue

        new_graph, final_trans = merge_flanks(
            [alt_locus.main_start_flank, alt_locus.start_flank,
             alt_locus.main_end_flank, alt_locus.end_flank],
            final_trans, new_graph, name_translation)
    return new_graph, final_trans


def convert_to_numeric_graph(graph):
    a_to_b = {k: i for i, k in enumerate(graph.blocks.keys())}
    trans = Translation.make_name_translation(a_to_b, graph)
    new_graph = trans.translate_subgraph(graph)
    trans.graph2 = new_graph

    # ?:
    for intervals in trans._a_to_b.values():
        for interval in intervals:
            interval.graph = new_graph

    return new_graph, trans


def create_initial_grch38_graph(chrom_sizes_fn):
    """
    Creates an initial grch38 graph with no connected blocks
    :param chrom_sizes_fn: Filename of chrom sizes file
    :return: Returns a Graph with no connected blocks
    """
    blocks = {}
    with open(chrom_sizes_fn, 'r') as csvfile:
        chroms = csv.reader(csvfile, delimiter='\t')
        for i, chrom in enumerate(chroms):
            blocks[chrom[0]] = Block(int(chrom[1]))

    return Graph(blocks, {})


def convert_to_text_graph(graph, name_translation, numeric_translation):
    for key in numeric_translation._b_to_a.keys():
        assert key in graph.blocks, "%s not in %s" % (key, graph)
    new_dict = {}

    # Set ids for rps in trans dict
    for i, key in enumerate(numeric_translation._b_to_a):
        rps = []

        # Get all region paths mapping to key
        for interval in numeric_translation._b_to_a[key]:
            rps.extend(interval.region_paths)
            new_id = str(i) + "".join(
                (name_translation._b_to_a[rp][0].region_paths[0] for rp in rps)
            )

        new_dict[key] = new_id

    # Set ids for rps not in trans dict
    for n_id in name_translation._b_to_a:
        if n_id not in numeric_translation._a_to_b:
            new_dict[n_id] = name_translation._b_to_a[n_id][0].region_paths[0]

    a_to_b = new_dict
    trans = Translation.make_name_translation(a_to_b, graph)
    new_graph = trans.translate_subgraph(graph)
    trans.graph2 = new_graph
    return new_graph, trans


def convert_cigar_graph_to_text(align_graph, name_translation,
                                numeric_translation):

    def create_name(blocks):
        names = [name_translation._b_to_a[block][0].region_paths[0]
                 for block in blocks]
        names = ["A" if "alt" in name else "M" for name in names]
        if len(names) == 2:  # HACK
            return "MA"

        return "".join(names)

    contributors = {block: [] for block in align_graph.blocks}
    for block, intervals in numeric_translation._b_to_a.items():
        for interval in intervals:
            contributors[block].extend(interval.region_paths)

    for block, conts in contributors.items():
        if not conts:
            contributors[block] = [block]

    a_to_b = {block: create_name(conts)+str(block)
              for block, conts in contributors.items()}
    graph1 = numeric_translation.graph2
    trans = Translation.make_name_translation(
        a_to_b, graph1)

    new_graph = trans.translate_subgraph(graph1)
    trans.set_graph2(new_graph)
    return new_graph, trans


def grch38_graph_to_numeric(original_grch38_graph):
    # Convert to numeric graph
    num_ab = {}
    num_ba = {}
    i = 0
    graph = original_grch38_graph.copy()
    for b in sorted(list(graph.blocks), reverse=True):
        i += 1
        num_ab[b] = [Interval(0, graph.blocks[b].length(), [i])]
        num_ba[i] = [Interval(0, graph.blocks[b].length(), [b])]
    trans = Translation(num_ab, num_ba, graph=original_grch38_graph)
    trans.graph2 = trans.translate_subgraph(original_grch38_graph)
    original_grch38_graph._update_a_b_graph(num_ab, trans.graph2)
    original_grch38_graph._update_a_b_graph(num_ba, original_grch38_graph)
    graph = trans.graph2

    return graph, trans


def merge_alt_using_cigar(original_numeric_grch38_graph,
                          trans, alt_id, ncbi_alignments_dir):
    """
    Uses a cigar string from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/
    and merges the alt locus.

    :param original_grch38_graph: Original numeric grch38 graph,
    built by calling the method create_initial_grch38_graph and then
    grch38_graph_to_numeric.
    :param trans: Translation object from original grch38 graph to
    numeric (returned by grch38_graph_to_numeric)
    :param alt_id: Alt id (e.g. chr2_KI270774v1_alt)
    :return: Returns the new graph and a translation object between
    the original grch38 graph and the new graph
    """
    # Find position of alt locus
    main_chr = ""
    main_start = 0
    main_end = 0
    alt_start = 0
    alt_end = 0

    try:
        with open(os.path.join(ncbi_alignments_dir, "%s.alignment" % alt_id)) as f:
            d = f.read().split(",")
            cigar = d[-1]
            # All numbers are 1-indxed and inclusive end (?)
            # Convert to 0-indexed and exclusive end
            main_start = int(d[0]) - 1
            main_end = int(d[1]) - 1 + 1  # Was inclusive
            alt_start = int(d[2]) - 1
            alt_end = int(d[3]) - 1 + 1
    except:
        print("Could not open alignment file alt_alignments/%s.alignment" %
              os.path.join(ncbi_alignments_dir, alt_id))
        return trans, trans.graph2
    main_chr = alt_id.split("_")[0]

    # Get sequences

    # uscs is 0 indexed and inclusive??
    alt_seq = get_sequence_ucsc(alt_id, alt_start + 1, alt_end)

    main_seq = get_sequence_ucsc(main_chr, main_start + 1, main_end)

    assert len(alt_seq) == (alt_end - alt_start)
    assert len(main_seq) == (main_end - main_start)

    from .cigar_align import clean_cigar, align_cigar
    cleaned_cigar = clean_cigar(cigar, alt_seq, main_seq)
    alt_id = trans.translate_rp(alt_id)[0].region_paths[0]
    main_id = trans.translate_rp(main_chr)[0].region_paths[0]

    trans = align_cigar(cleaned_cigar,
                        Interval(main_start, main_end, [main_id]),
                        Interval(alt_start, alt_end, [alt_id]),
                        original_numeric_grch38_graph)

    new_graph = trans.translate_subgraph(original_numeric_grch38_graph)
    trans.set_graph2(new_graph)

    return trans, new_graph


def _merge_alt_using_cigar(original_grch38_graph, trans, alt_id, cigar,
                           alt_seq, main_seq, main_chr, main_pos_start,
                           main_pos_end, alt_start, alt_end):
    """

    :param original_grch38_graph: Numeric grch38 graph (created by calling grch38_graph_to_numeric)
    :param trans: Trans returned by grch38_graph_to_numeric()
    :param alt_id:
    :param cigar:
    :param alt_seq:
    :param main_seq:
    :param main_chr:
    :param main_pos_start: Position of first base pair after alt loci starts
    :param main_pos_end: Position of last base pair before alt loci merges
    :param alt_length:
    :return:
    """

    # Put sequences into numpy arrays
    import numpy as np
    alt_seq = np.array(list(alt_seq))
    main_seq = np.array(list(main_seq))

    length = main_pos_end - main_pos_start

    """
    Algorithm:
        go through cigar. For every cigar element, modify graph.
        If cigar is m,
    """
    cigar = cigar.split(" ")
    main_offset = main_pos_start
    alt_offset = alt_start
    graph = original_grch38_graph.copy()

    for c in cigar:
        first = c[0]
        if first == "M" or first == "D" or first == "I":
            type = first
            n = int(c[1:])
        else:
            type = c[-1]
            n = int(c[:-1])

        if type == "I":
            alt_offset += n
            continue
        if type == "D":
            main_offset += n
            continue

        # If sequences are identical, then merge
        assert type == "M"
        seq_on_alt = alt_seq[alt_offset-alt_start:alt_offset-alt_start+n]
        seq_on_main = main_seq[main_offset-main_pos_start:main_offset+n-main_pos_start]
        prev = 0
        for i in range(0, len(seq_on_alt)):
            if seq_on_alt[i] != seq_on_main[i] or i == n - 1:
                if i <= prev:
                    continue

                if i == n-1:
                    i += 1

                # Something is different, merge from prev to i
                m_start = main_offset + prev
                m_end = main_offset + i
                a_start = alt_offset + prev
                a_end = alt_offset + i
                intv1 = trans.translate(Interval(m_start, m_end, [main_chr],
                                                 original_grch38_graph))
                intv2 = trans.translate(Interval(a_start, a_end, [alt_id],
                                                 original_grch38_graph))
                graph, mtrans = graph.merge([intv1, intv2])

                trans = trans + mtrans
                trans.graph2 = graph
                prev = i+1
        main_offset += n
        alt_offset += n

    # Should be correct if cigar is correct:

    if alt_offset - alt_start != len(alt_seq):
        print("alt offset %d != length of sequence %d" % (alt_offset - alt_start, len(alt_seq)))
    else:
        print("Alt offset correct")

    if main_offset - main_pos_start != len(main_seq):
        print("main offset %d != length of sequence %d" % (main_offset - main_pos_start, len(main_seq)))
    else:
        print("Main offset correct!")

    assert alt_offset - alt_start == len(alt_seq)
    assert main_offset - main_pos_start == len(main_seq)

    # Final translation should make final graph from original
    assert graph == trans.graph2 # (original_grch38_graph)

    return trans, graph
