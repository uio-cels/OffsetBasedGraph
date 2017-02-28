from .interval import Interval
from .graphcreators import *
from .genematcher import GeneMatchings
import csv
from offsetbasedgraph import CriticalPathsMultiPathInterval,\
    FuzzyMultipathInterval
from collections import defaultdict
from .gene import GeneList, Gene, MultiPathGene, FuzzyGene
import sys


def create_subgraph_around_alt_locus(graph, trans, alt_locus, padding=200000,
                                     alt_loci_fn='grch38_alt_loci.txt'):
    """
    Takes a translation from an original grch38 graph
    (created by create_initial_grch38_graph)
    and uses a translation from that graph to current graph
    in order to filter out the subgraph around a particular alt locus.
    :param graph: Graph to create subgraph from
    :param trans: Translation from original grch38 graph to current graph
    :param alt_locus: alt locus ID to create subgraph around
    :param padding: Padding. Number if basepairs to include around alt locus
    :return: Returns subgraph, a translation object and
    start position of the subgraph
    """

    # For now, only use trans object to find which intervals in graph
    # that can be used to create subgraph by calling
    # Graph.create_subgraph_from_intervals()
    # todo: This can be rewritten more efficiently later

    original_graph = trans.graph1
    assert graph == trans.graph2, "trans.graph2 should be identical to graph"
    alt_locus_interval = Interval(
        0, original_graph.blocks[alt_locus].length(), [alt_locus])

    from offsetbasedgraph.GRCH38 import AltLoci
    alt_loci_info = AltLoci.from_file(alt_loci_fn, [alt_locus])
    alt_locus_info = alt_loci_info.lookup[alt_locus]

    main_interval_start = alt_locus_info.start - 10
    main_interval_end = alt_locus_info.end + 10
    main_interval = Interval(main_interval_start,
                             main_interval_end,
                             [alt_locus_info.chrom])

    # Translate interval, get all intervals
    translated_to = []
    for orig_interval in [main_interval, alt_locus_interval]:
        trans_interval = trans.translate(orig_interval)
        translated_to.append(trans_interval)

    return graph.create_subgraph_from_intervals(
        translated_to, padding, alt_locus)


def parse_genes_file(genes_fn):
    """
    Parses a file containing genes (on the format of UCSC),
    and returns a list of dicts
    :param genes_fn: File name
    :return: Returns a list of genes, one dict for each gene
    """
    genes = []
    with open(genes_fn) as f:
        g = csv.DictReader(f, delimiter='\t')
        for gene in g:
            genes.append(gene)

    return genes


def get_gene_objects_as_intervals(fn, graph=None):
    """
    Returns a dict. Keys are gene names and values are intervals
    representing the gene
    :param fn: File name of file containing genes (on the format of UCSC)
    :return:
    """
    genes = parse_genes_file(fn)
    return [Gene.from_dict(gene) for gene in genes]


def analyze_genes_on_merged_graph(genes, translation):
    """
    Find similarities between gene annotations on alt-loci and
    main chromosomes when represented on a merged graph
    Translate the genes to the merged graph and perform analysis
    :param genes: list of genes on original graph
    :param translation: translation from original to merged graph
    """
    cashed = False
    translation.block_lengths = None
    if cashed:
        translated = GeneList.from_pickle("trans_genes").gene_list
    else:
        translated = [gene.translate(translation) for gene in genes]
        gene_list = GeneList(translated)
        gene_list.to_file("trans_genes")

    graph = translated[0].transcription_region.graph
    graph.critical_blocks = graph.find_all_critical_blocks()
    alt_genes = []
    main_genes = []
    for gene, t_gene in zip(genes, translated):
        if "alt" in gene.chrom:
            alt_genes.append(t_gene)
        else:
            main_genes.append(t_gene)
    alt_genes = GeneList(alt_genes)
    main_genes = GeneList(main_genes)
    matchings = GeneMatchings(alt_genes, main_genes)
    print(matchings)


def get_alt_loci_positions(alt_loci_file_name):
    """Create a dict of alt loci positions read from
    alt_loci_file_name

    :param alt_loci_file_name: file name
    :returns: dict{alt_loci_id: info (dict)}
    :rtype: dict{str: dict}

    """

    f = open(alt_loci_file_name, "r")
    alt_loci = {}
    for line in f.readlines():
        if line.startswith("#"):
            continue
        l = line.split()
        alt_locus_id = l[0]
        main_chr = l[1]
        start = int(l[2])
        end = int(l[3])
        length = int(l[4])
        alt_loci[alt_locus_id] = {
                "main_chr": main_chr,
                "start": start,
                "end": end,
                "length": length
             }

    f.close()
    return alt_loci


def create_gene_dicts(genes, alt_loci_fn, identical_names=False):
    """
    Takes a list of genes, and creates an alt loci dict and a gene name dict.
    The dicts index the genes on chromosome/alt_loci id
    :param genes: list of genes
    :return: alt loci dict, name dict, paralell_dict
    :rtype: defaultdict, defaultdict, defaultdict
    """
    gene_name_dict = defaultdict(list)
    for g in genes:
        gene_name_dict[g.name].append(g)

    # Find all genes on alt loci
    alt_loci_genes = defaultdict(list)
    chrom_genes = defaultdict(list)
    n = 0
    for g in genes:
        n += 1
        chrom = g.transcription_region.region_paths[0]
        if "alt" in chrom:
            alt_loci_genes[chrom].append(g)
        else:
            # Index non-alt loci genes using  offset of first exon
            chrom_genes[chrom].append(g)

    alt_infos = get_alt_loci_positions(alt_loci_fn)
    parallell_to_alt_loci_genes = defaultdict(list)  # Genes on  main path
    for alt_locus in alt_loci_genes:
        main_genes = chrom_genes[alt_locus.split("_")[0]]
        # Find main genes parallell to alt locus
        alt_info = alt_infos[alt_locus]
        for g in main_genes:
            start = g.transcription_region.start_position.offset
            end = g.transcription_region.end_position.offset

            if (start >= alt_info["start"] and start <= alt_info["end"]) or \
                    (end <= alt_info["end"] and end >= alt_info["start"]):
                parallell_to_alt_loci_genes[alt_locus].append(g)

    return alt_loci_genes, gene_name_dict, parallell_to_alt_loci_genes


def translate_to_fuzzy_interval(gene, trans):
    """Translate the gene to a fuzzy interval on graph2 of trans

    :param gene: Gene
    :param trans: Translation
    :rtype: FuzzyGene

    """
    transcription_region = trans.translate(
        gene.transcription_region)
    exons = [trans.translate(exon) for exon in gene.exons]
    f_t_region = FuzzyMultipathInterval.from_interval(
        transcription_region, 5)
    f_exons = [FuzzyMultipathInterval.from_interval(exon, 5)
               for exon in exons]
    f_gene = FuzzyGene(gene.name, f_t_region, f_exons, gene.coding_region,
                       gene.strand)
    return f_gene


def translate_single_gene_to_aligned_graph(gene, trans):
    """Translate the gene to a multipath gene on trans.graph2

    :param gene: Gene
    :param trans: Translation
    :rtype: MultiPathGene

    """
    start = gene.transcription_region.start_position
    end = gene.transcription_region.end_position
    end.offset = end.offset - 1
    new_start = trans.translate(start)
    new_end = trans.translate(end)
    end.offfset = end.offset + 1
    new_end.offset = new_end.offset + 1
    critical_intervals = [trans.translate(exon) for exon in gene.exons]

    mpinterval = CriticalPathsMultiPathInterval(
        new_start,
        new_end,
        critical_intervals
    )

    return MultiPathGene(gene.name, mpinterval)


def _analyse_fuzzy_genes_on_graph(genes_list, genes_against, graph):
    """Takes a list of mp genes and a graph
    Returns number of equal exons and equal genes
    """

    matches = [(gene, against_gene) for gene in genes_list
               for against_gene in genes_against if
               gene == against_gene]

    for a, b in matches:
        print("-----")
        print(a)
        print(b)
    return len(matches)


def analyze_multipath_genes_for_alt(alt_id, alt_loci_genes, main_genes,
                                    graph, name_trans, ncbi_alignments_dir):
    """Find genes with the same fuzzy multipath interval representation

    :param alt_id: Alt locus id
    :param alt_loci_genes: Genes on alt loci
    :param main_genes: Genes on main chromosome paralell to alt loci
    :param graph: Graph
    :param name_trans: Translation
    :param ncbi_alignments_dir: Directory to find alt-loci info
    :returns: Number of genes on alt-loci with match on main chromsome
    :rtype: int

    """
    print("Analysing genes on alt locus %s" % alt_id)
    genes_here = alt_loci_genes[alt_id]
    if not (genes_here and main_genes[alt_id]):
        return 0

    trans, complex_graph = merge_alt_using_cigar(graph, name_trans, alt_id,
                                                 ncbi_alignments_dir)
    full_trans = name_trans + trans

    # Find candidates on main path to check against:
    genes_against = [g.copy() for g in main_genes[alt_id]]
    genes_against_translated = []
    n = 0

    def report_progress(result, i, n):
        sys.stdout.write('\r  Translating main genes: ' +
                         str(round(100 * i / max(1, n))) +
                         ' % finished ' + ' ' * 20)
        sys.stdout.flush()
        return result

    n = len(genes_against)
    genes_against_translated = [report_progress(
        translate_to_fuzzy_interval(mg, full_trans), i, n)
                        for i, mg in enumerate(genes_against)]
    print()
    n = len(genes_here)
    genes_here_translated = [report_progress(
        translate_to_fuzzy_interval(mg, full_trans), i, n)
                             for i, mg in enumerate(genes_here)]
    print()
    return _analyse_fuzzy_genes_on_graph(
        genes_here_translated,
        genes_against_translated,
        complex_graph)


def fuzzy_gene_analysis(genes, text_graph, ncbi_alignments_dir,
                        alt_loci_filename):
    """Find number of genes on alt-loci that can be represented
    by the fuzzy multipath interval of a gene on a main chromsome
    when using a complex graph

    :param genes: List of genes
    :param text_graph: Graph
    :param ncbi_alignments_dir: Directory alignments
    :param alt_loci_filename: File with alt loci info
    """
    print("Readinng in genes")
    alt_loci_genes, gene_name_dict, main_genes = create_gene_dicts(
        genes, alt_loci_filename)
    graph, name_trans = grch38_graph_to_numeric(text_graph)
    equal_total = 0
    for b in text_graph.blocks:
        if "alt" not in b:
            continue
        equal_total += analyze_multipath_genes_for_alt(
            b, alt_loci_genes, main_genes, graph, name_trans,
            ncbi_alignments_dir)

        print()

    print("RESULTS:")
    print("%d genes on alternative loci have identical representation as at least one gene from the main chromosome." % equal_total)
    print("In total %d genes on alt loci" % len(alt_loci_genes))

