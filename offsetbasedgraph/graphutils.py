from .interval import Interval, Position
from .graph import Graph, Block
from .translation import Translation
import csv
import pickle
from genutils import flanks
from gendatafetcher.sequences import get_sequence_ucsc

from collections import defaultdict


class GeneList(object):
    def __init__(self, gene_list):
        self.gene_list = gene_list

    def to_file(self, file_name):
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def from_file(file_name):

        with open("%s" % file_name, "rb") as f:
            return pickle.loads(f.read())


class Gene(object):

    def __init__(self, name, transcription_region, exons):
        self.name = name
        self.transcription_region = transcription_region
        self.exons = exons
        self.chrom = transcription_region.region_paths[0]

    @classmethod
    def from_dict(cls, attr_dict):
        """Create Gene object from gene dict obtained
        from csv dict reader

        :param attr_dict: gene dict
        :returns: Gene object
        :rtype: Gene

        """
        chrom = attr_dict["chrom"]
        transcription_region = Interval(
            Position(chrom, int(attr_dict["txStart"])),
            Position(chrom, int(attr_dict["txEnd"])), [chrom]
            )

        exon_starts = [int(i) for i in attr_dict["exonStarts"].split(",")[:-1]]
        exon_ends = [int(i) for i in attr_dict["exonEnds"].split(",")[:-1]]
        exons = [Interval(start, end, [chrom]) for start, end in
                 zip(exon_starts, exon_ends)]
        return cls(attr_dict["name"], transcription_region, exons)

    def translate(self, T):
        """Translate transcription_region and exons and return
        new gene

        :param T: Translation object
        :returns: Translated gene
        :rtype: Gene

        """

        t_transcription_region = T.translate(self.transcription_region)
        t_exons = [T.translate(exon) for exon in self.exons]
        return Gene(self. name, t_transcription_region, t_exons)

    def __str__(self):
        exon_string = "\n\t".join(str(exon) for exon in self.exons)
        return "Gene(%s: %s \n %s)" % (self.name,
                                       self.transcription_region, exon_string)

    def __eq__(self, other):
        """Check if genes are equal up to name

        :param other: other gene
        :returns: Wheter genes are equal
        :rtype: bool

        """
        if not self.transcription_region == other.transcription_region:
            return False

        if len(self.exons) != len(other.exons):
            return False

        for exon in self.exons:
            if exon not in other.exons:
                return False

        return True


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
            new_id = str(i) + "".join((name_translation._b_to_a[rp][0].region_paths[0] for rp in rps))

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
    #if final_trans is not None:
    #print("=== Flanking intervals ===")
    #print(intervals)
    #print('\n'.join([str(i) for i in intervals]))
    # Merge start flank of alt locus with main
    merge_intervals = intervals[0:2]
    merge_intervals = [name_translation.translate(i) for i in merge_intervals]
    merge_intervals = [final_trans.translate(i) for i in merge_intervals]
    #print("=== Flanking intervals after translate===")
    #print('\n'.join([str(i) for i in merge_intervals]))
    #print("=== Graph ===")
    #print(merge_intervals[0].graph)
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

        #print("=== end flank trans ===")
        #print(trans)

        final_trans += trans
        final_trans.graph2 = new_graph.copy()
        final_trans.graph2._update_a_b_graph(final_trans._a_to_b, new_graph)
        #print("No end flank")

    return new_graph, final_trans


def connect_without_flanks(graph, alt_loci_fn, name_translation):
    """
    Connects the alternative loci in the given file to the grch38 graph,
    without flanks.
    :param alt_loci_fn: Filename of file containing alternative loci.
    One alt locus on each line. Four columns: alt_locus_id  chr chr_start   chr_stop
    :return: Returns the new graph
    """
    f = open(alt_loci_fn)
    n_flanks = 0
    new_graph = graph
    orig_graph = graph.copy()
    final_trans = Translation(graph=graph)
    final_trans.graph2 = graph
    for line in f.readlines():
        if line.startswith("#"):
            continue
        l = line.split()
        alt_locus_id = l[0]
        #print("Connecting %s" % (alt_locus_id))
        main_chr = l[1]
        start = int(l[2])
        end = int(l[3])
        length = int(l[4])

        intervals = flanks.get_flanks(alt_locus_id, length,
                                      main_chr, start-1, end)
        new_graph, final_trans = merge_flanks(intervals, final_trans,
                                              new_graph, name_translation)

    f.close()
    #print("NUMBER OF FLANKS: %d" % n_flanks)
    return new_graph, final_trans


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


def get_gene_objects_as_intervals(fn, graph):
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
    * Write diff for intervals
    * Write diff for genes
    * Numeric/

    :param main_genes: 
    :param alt_genes: 
    :returns: 
    :rtype: 

    """
    graph = main_genes[0].transcription_region.graph
    graph.critical_blocks = graph.find_all_critical_blocks()
    main_dict = defaultdict(list)
    for gene in main_genes:
        main_dict[gene.name].append(gene)
    for gene in alt_genes:
        if gene.name in main_dict and gene not in main_dict[gene.name]:
            codes = []
            for m_gene in main_dict[gene.name]:
                m_code = m_gene.transcription_region.diff(
                    gene.transcription_region)
                if all(c == "E" for c in m_code):
                    codes.append("("+"".join(m_code)+")")
            if codes:
                print(gene.transcription_region.region_paths)
                print(gene.name, "".join(codes))


def find_exon_duplicates(genes, translation):
    """Find and count duplicate genes on flanks of
    alt loci

    :param genes: genes on original graph
    :param translation: translation object
    """
    translated = GeneList.from_file("trans_genes").gene_list
    # [gene.translate(translation) for gene in genes]
    # gene_list = GeneList(translated)
    # gene_list.to_file("trans_genes")
    main_chr_dict = defaultdict(list)
    alt_dict = defaultdict(list)

    for gene, t_gene in zip(genes, translated):
        if "alt" in gene.chrom:
            alt_dict[gene.chrom.split("_")[0]].append(t_gene)
        else:
            main_chr_dict[gene.chrom].append(t_gene)

    s = 0
    for chrom, genes in alt_dict.items():
        print(chrom)
        find_unequal_sibling_genes(main_chr_dict[chrom], genes)
        continue
        for gene in genes:
            for main_gene in main_chr_dict[chrom]:
                if gene == main_gene:
                    print(main_gene.name, gene.name)
                    s += 1

    print(s, "/", sum(len(v) for v in alt_dict.values()))


def blast_test():



    from Bio.Blast import NCBIWWW
    from Bio.Blast.Applications import NcbiblastpCommandline
    from Bio.Blast import NCBIXML
    try:
        from StringIO import StringIO  # Python 2
    except ImportError:
        from io import StringIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    q = "data/tmp/sequence_chr2_KI270769v1_alt_1_120616.fasta"
    r = "data/tmp/sequence_chr4_KI270785v1_alt_1_119912.fasta"

    from Bio import pairwise2
    qs = open(q).read()
    rs = open(r).read()
    alignments = pairwise2.align.globalxx(qs, rs)
    print(alignments)
    """
    # Run BLAST and parse the output as XML
    output = NcbiblastpCommandline(query=q, cmnd='blastn', subject=r, outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))

    # Print some information on the result
    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
    """

def merge_alt_using_cigar(original_grch38_graph, alt_id):
    """
    Uses a cigar string from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/
    and merges the alt locus.

    :param original_grch38_graph: Original grch38 graph, built by calling the method create_initial_grch38_graph
    :param alt_id: Alt id (e.g. chr2_KI270774v1_alt)
    :return: Returns the new graph and a translation object between the original grch38 graph and the new graph
    """

    try:
        with open("alt_alignments/%s.alignment" % alt_id) as f:
            cigar = f.read()
    except:
            print("Could not open alignment file alt_alignments/%s.alignment" % alt_id)


     # Find position of alt locus
    main_chr = ""
    start = 0
    end = 0
    length = 0

    alt_loci_fn = "grch38_alt_loci.txt"
    f = open(alt_loci_fn)
    for line in f.readlines():
        if line.startswith("#"):
            continue
        l = line.split()
        alt_locus_id = l[0]
        if alt_locus_id == alt_id:
            main_chr = l[1]
            start = int(l[2])
            end = int(l[3])
            length = int(l[4])
            break

    # Get sequences
    alt_seq = get_sequence_ucsc(alt_id, 1, length)
    main_seq = get_sequence_ucsc(main_chr, start, end)

    assert len(alt_seq) == length

    print(alt_seq[-10:])
    print(main_seq[-10:])

    print(main_chr, start, end, length)

    return _merge_alt_using_cigar(original_grch38_graph, alt_id, cigar, alt_seq, main_seq, main_chr, start, end, length)

def _merge_alt_using_cigar(original_grch38_graph, alt_id, cigar, alt_seq, main_seq, main_chr, main_pos_start, main_pos_end, alt_length):
    """

    :param original_grch38_graph:
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

    """
    Algorithm:
        go through cigar. For every cigar element, modify graph.
        If cigar is m,
    """
    print(cigar)
    cigar = cigar.split(" ")
    alt_offset = 0
    main_offset = main_pos_start

    ab = {}
    ba = {}
    print("=== Merge alt ===")
    graph = original_grch38_graph.copy()

    # Convert to numeric graph
    num_ab = {}
    num_ba = {}
    i = 0
    for b in graph.blocks:
        i += 1
        num_ab[b] = [Interval(0, graph.blocks[b].length(), [i])]
        num_ba[i] = [Interval(0, graph.blocks[b].length(), [b])]
    trans = Translation(num_ab, num_ba, graph=original_grch38_graph)
    trans.graph2 = trans.translate_subgraph(original_grch38_graph)
    original_grch38_graph._update_a_b_graph(num_ab, trans.graph2)
    original_grch38_graph._update_a_b_graph(num_ba, original_grch38_graph)
    graph = trans.graph2


    print("== Numeric graph ===")
    print(graph)
    """
    graph, t = graph.connect_postitions(
                        trans.translate(Position(main_chr, main_pos_start - 1)),
                        trans.translate(Position(alt_id, 0)))
    trans += t
    graph, t = graph.connect_postitions(
                        trans.translate(Position(alt_id, alt_length - 1)),
                        trans.translate(Position(main_chr, main_pos_end + 1)))
    trans += t
    """

    print("=== Graph before cigar ===")
    print(graph)

    print(cigar)
    for c in cigar:

        first = c[0]
        if first == "M" or first == "D" or first == "I":
            type = first
            n = int(c[1:])
        else:
            type = c[-1]
            n = int(c[:-1])

        print("Type/len: %s/%d" % (type, n))

        if type == "M":
            #print("Match")
            # If sequences are identical, then merge
            if alt_seq[alt_offset:alt_offset+n] == main_seq[main_offset-main_pos_start:main_offset+n-main_pos_start]:
                print("Merging match")
                intv1 = trans.translate(Interval(main_offset, main_offset+n, [main_chr], graph))
                intv2 = trans.translate(Interval(alt_offset, alt_offset+n, [alt_id], graph))
                #print("Merging %s with %s" % (intv1, intv2))
                graph, mtrans = graph.merge([intv1, intv2])

                trans = trans + mtrans
                trans.graph2 = graph

                #print(trans)
                #print(graph)

            else:
                print("Not merging, different")
                #continue
                # Just connect ?
                #graph.connect_postitions()

            main_offset += n
            alt_offset += n
            #print("Offsets alt/main: %d/%d" % (alt_offset, main_offset))

        elif type == "I":
            alt_offset += n
        elif type == "D":
            main_offset += n

    # Should be correct if cigar is correct:
    assert alt_offset == len(alt_seq)
    assert main_offset - main_pos_start == len(main_seq)

    # Final translation should make final graph from original
    assert graph == trans.translate_subgraph(original_grch38_graph)


    return trans, graph

