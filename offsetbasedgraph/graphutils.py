from .interval import Interval
from .graphcreators import *
from .translation import Translation
from .genematcher import GeneMatcher, GeneMatchings
import csv
from genutils import flanks
from offsetbasedgraph.sequences import get_sequence_ucsc
from offsetbasedgraph import CriticalPathsMultiPathInterval,\
    FuzzyMultipathInterval
from collections import defaultdict
from .gene import GeneList, Gene, MultiPathGene, FuzzyGene
import sys


def _connect_without_flanks(graph, alt_loci_fn, name_translation):
    """
    Connects the alternative loci in the given file to the grch38 graph,
    without flanks.
    :param alt_loci_fn: Filename of file containing alternative loci.
    One alt locus on each line.
    Four columns: alt_locus_id  chr chr_start   chr_stop
    :return: Returns the new graph
    """
    f = open(alt_loci_fn)
    new_graph = graph
    final_trans = Translation(graph=graph)
    final_trans.graph2 = graph
    for line in f.readlines():
        if line.startswith("#"):
            continue
        l = line.split()
        alt_locus_id = l[0]
        main_chr = l[1]
        start = int(l[2])
        end = int(l[3])
        length = int(l[4])

        intervals = flanks.get_flanks(alt_locus_id, length,
                                      main_chr, start-1, end)

        n_starts = len([b for b in new_graph.blocks if not
                        new_graph.reverse_adj_list[b]])
        new_graph, final_trans = merge_flanks(intervals, final_trans,
                                              new_graph, name_translation)

        new_n_starts = len([b for b in new_graph.blocks if not
                            new_graph.reverse_adj_list[b]])
        assert new_n_starts == n_starts-1, line

    f.close()
    return new_graph, final_trans


def create_subgraph_around_alt_locus(graph, trans, alt_locus, padding=200000,
                                     alt_loci_fn='grch38_alt_loci.txt'):
    """
    Takes a translation from an original grch38 graph (created by create_initial_grch38_graph)
    and uses a translation from that graph to current graph
    in order to filter out the subgraph around a particular alt locus.
    :param graph: Graph to create subgraph from
    :param trans: Translation from original grch38 graph to current graph
    :param alt_locus: alt locus ID to create subgraph around
    :param padding: Padding. Number if basepairs to include around alt locus
    :return: Returns subgraph, a translation object and start position of the subgraph
    """

    # For now, only use trans object to find which intervals in graph
    # that can be used to create subgraph by calling
    # Graph.create_subgraph_from_intervals()
    # todo: This can be rewritten more efficiently later

    original_graph = trans.graph1
    assert graph == trans.graph2 , "trans.graph2 should be identical to graph"
    alt_locus_interval = Interval(0, original_graph.blocks[alt_locus].length(), [alt_locus])

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
        #print("Orig interval")
        #print(orig_interval)
        trans_interval = trans.translate(orig_interval)
        #print("Including %s" % trans_interval)
        translated_to.append(trans_interval)

    return graph.create_subgraph_from_intervals(translated_to, padding, alt_locus)


def parse_genes_file(genes_fn):
    """
    Parses a file containing genes (on the format of UCSC), and returns a list of dicts
    :param genes_fn: File name
    :return: Returns a list of genes, one dict for each gene
    """
    genes = []
    with open(genes_fn) as f:
        g = csv.DictReader(f, delimiter='\t')
        for gene in g:
            genes.append(gene)

    return genes


def get_genes_as_intervals(fn, graph):
    """
    Returns a dict. Keys are gene names and values are intervals
    representing the gene
    :param fn: File name of file containing genes (on the format of UCSC)
    :return:
    """
    genes = parse_genes_file(fn)
    out = {}
    for gene in genes:
        chrom = gene["chrom"]
        start = int(gene["txStart"])
        end = int(gene["txEnd"])
        name = gene["name"]
        out[name] = Interval(start, end, [chrom], graph)

    return out


def get_gene_objects_as_intervals(fn, graph=None):
    """
    Returns a dict. Keys are gene names and values are intervals
    representing the gene
    :param fn: File name of file containing genes (on the format of UCSC)
    :return:
    """
    genes = parse_genes_file(fn)
    return [Gene.from_dict(gene) for gene in genes]


def find_unequal_sibling_genes(main_genes, alt_genes):
    """TODO:
    Categories:
        Completely equal
        Equal up to base pairs
        Fully parallel
        Starting ending equally, up to base pairs
        Starting/Ending equally
        Starting/Ending equally up to base pairs

    :param main_genes: 
    :param alt_genes: 
    :returns: 
    :rtype: 

    """
    graph = main_genes[0].transcription_region.graph
    main_dict = defaultdict(list)
    gene_categories = defaultdict(list)
    for gene in main_genes:
        main_dict[gene.name].append(gene)
    gene_matches = []
    for gene in alt_genes:
        if gene.name not in main_dict:
            continue
        gene_matches.append(GeneMatcher(gene, main_dict[gene.name]))

    return gene_matches


def classify_alt_gene(gene):
    def is_merged(name):
        return name.count("chr") > 1

    def was_alt(name):
        return "alt" in name

    def is_varying(name):
        return was_alt(name) and not is_merged(name)

    rps = gene.transcription_region.region_paths
    if not any(is_varying(rp) for rp in rps):
        return "FLANK"
    if is_merged(rps[0]) and is_merged(rps[-1]):
        return "BUBBLE"
    elif is_merged(rps[0]):
        return "START"
    elif is_merged(rps[-1]):
        return "END"
    if len(rps) == 1:
        return "ALT"
    else:
        raise Exception("Couldnt classify %s" % rps)


def classify_alt_genes(genes):
    """Group alt genes by their position on the alt loci
    i.e, on flank, varying or a combination.

    :param genes: list of  alt genes
    :returns: dict of gene-lists for each category
    :rtype: dict[str] = list(Gene)

    """

    rp_lists = [gene.transcription_region.region_paths for gene in genes]
    classificaitons = []
    for rps in rp_lists:
        class_lists = defaultdict(list)

    for gene, classi in zip(genes, classificaitons):
        class_lists[classi].append(gene)

    return class_lists


def analyze_genes_on_merged_graph(genes, translation):
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
    main_chr_dict = defaultdict(list)

    alt_dict = defaultdict(list)
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
    return

    print("N:", sum(len(v) for v in alt_dict.values()))
    gene_matches = []
    for chrom, genes in alt_dict.items():
        gene_matches.extend(
            find_unequal_sibling_genes(main_chr_dict[chrom], genes))

    print("M", len(gene_matches))
    for gene_category in ["ALT", "FLANK", "START", "END", "BUBBLE"]:
        category_mathes = [m for m in gene_matches
                           if m.category == gene_category]
        score_dict = defaultdict(int)
        for match in category_mathes:
            score_dict[match.score] += 1
        print(gene_category, score_dict)


def find_exon_duplicates(genes, translation):
    """Find and count duplicate genes on flanks of
    alt loci

    :param genes: genes on original graph
    :param translation: translation object
    """
    translation.block_lengths = None
    translated = [gene.translate(translation) for gene in genes]
    gene_list = GeneList(translated)
    gene_list.to_file("trans_genes")
    print([g.transcription_region for g in genes])
    print([g.transcription_region for g in translated])
    # translated = GeneList.from_file("trans_genes").gene_list
    graph = translated[0].transcription_region.graph
    graph.critical_blocks = graph.find_all_critical_blocks()

    main_chr_dict = defaultdict(list)

    # class_lists = classify_alt_genes(genes)
    alt_dict = defaultdict(list)
    gene_categories = defaultdict(list)

    for gene, t_gene in zip(genes, translated):
        if "alt" in gene.chrom:
            alt_dict[gene.chrom.split("_")[0]].append(t_gene)
        else:
            main_chr_dict[gene.chrom].append(t_gene)

    print("N:", sum(len(v) for v in alt_dict.values()))
    s = 0
    gene_matches = []
    for chrom, genes in alt_dict.items():
        gene_matches.extend(
            find_unequal_sibling_genes(main_chr_dict[chrom], genes))
        continue
        for k, v in new_things.items():
            gene_categories[k].extend(v)
        continue
        for gene in genes:
            for main_gene in main_chr_dict[chrom]:
                if gene == main_gene:
                    print(main_gene.name, gene.name)
                    s += 1


def blast_test():
    try:
        from StringIO import StringIO  # Python 2
    except ImportError:
        from io import StringIO

    q = "data/tmp/sequence_chr2_KI270769v1_alt_1_120616.fasta"
    r = "data/tmp/sequence_chr4_KI270785v1_alt_1_119912.fasta"

    from Bio import pairwise2
    qs = open(q).read()
    rs = open(r).read()
    alignments = pairwise2.align.globalxx(qs, rs)
    print(alignments)


def get_alt_loci_positions(alt_loci_file_name):
    try:
        f = open(alt_loci_file_name, "r")
    except:
        print("Error: %s does not exist." % alt_loci_file_name)
        raise

    alt_loci = {}
    for line in f.readlines():
        print(line)
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

    return alt_loci


def create_gene_dicts(genes, alt_loci_fn,
                      identical_names=False):
    """
    Takes a list of genes, and creates an alt loci dict and a gene name dict
    :param genes:
    :return: alt loci dict, name dict
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
    # Takes a list of mp genes and a graph
    # Returns number of equal exons and equal genes

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
    from offsetbasedgraph.graphutils import create_gene_dicts
    print("Readinng in genes")

    alt_loci_genes, gene_name_dict, main_genes = create_gene_dicts(
        genes, alt_loci_filename)

    # alt loci genes are only genes on alt loci (nothing on main)
    # exon_dict contains only genes on main, index by offset of first exon

    # For every alt loci create complex graph, translate genes and analyse them
    graph, name_trans = grch38_graph_to_numeric(text_graph)

    equal_total = 0
    # for b in text_graph.blocks:
    for b in text_graph.blocks:
        if "alt" not in b:
            continue
        equal_total += analyze_multipath_genes_for_alt(
            b, alt_loci_genes, main_genes, graph, name_trans,
            ncbi_alignments_dir)
        print(equal_total)
        print()

    print(equal_total)
