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
        assert isinstance(gene_list, list)
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

    def length(self):
        return self.transcription_region.length()

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

    def to_file_line(self):
        cur_rp = self.transcription_region.region_paths[0]
        exon_strings = []
        for exon in self.exons:
            rps = [rp for rp in exon.region_paths if rp != cur_rp]
            offsets = [exon.start_position.offset, exon.end_position.offset]
            exon_string = ",".join(str(e) for e in rps+offsets)
            exon_strings.append(exon_string)
            cur_rp = exon.region_paths[-1]
        transcript_string = ",".join(self.transcription_region.region_paths+[
            str(self.transcription_region.start_position.offset),
            str(self.transcription_region.end_position.offset)])
        return "\t".join([self.name,
                          transcript_string,
                          " ".join(exon_strings)])

    @classmethod
    def parse_transcript_region(cls, transcript_string):
        parts = transcript_string.split(",")
        start_offset, end_offset = [int(n) for n in parts[-2:]]
        region_paths = parts[:-2]
        return Interval(start_offset, end_offset, region_paths)

    @classmethod
    def parse_exons(cls, exon_string, start_rp):
        exons = exon_string.split(" ")
        for exon in exons:
            start_offset, end_offset = [int(n) for n in parts[-2:]]

    @classmethod
    def from_file_line(cls, line):
        name, transcript_string, exons_string = line.split("\t")
        transcription_region = cls.parse_transcript_region(transcript_string)
        exons = cls.parse_exons(exons_string)

    def start_and_end_diffs(self, other):
        my_interval = self.transcription_region
        other_interval = other.transcription_region
        start_diff = 0
        end_diff = 0
        if my_interval.region_paths[0] == other_interval.region_paths[0]:
            start_diff = abs(my_interval.start_position.offset - other_interval.start_position.offset)
        else:
            start_diff = -1

        if my_interval.region_paths[-1] == other_interval.region_paths[-1]:
            end_diff = abs(my_interval.end_position.offset - other_interval.end_position.offset)
        else:
            end_diff = -1
        return (start_diff, end_diff)

    def exon_diffs(self, other):
        if len(self.exons) != len(other.exons):
            return None
        diff_sum = 0
        for exon1, exon2 in zip(self.exons, other.exons):
            if exon1.region_paths == exon2.region_paths:
                s_diff = exon1.start_position.offset-exon2.start_position.offset
                e_diff = exon1.end_position.offset-exon2.end_position.offset
                diff_sum += abs(s_diff)+abs(e_diff)
            else:
                diff_sum += abs(exon1.length()-exon2.length())

        return diff_sum/len(self.exons)

    def contains(self, other, tolerance=0):
        if not self.transcription_region.contains(
                other.transcription_region, tolerance):
            return False
        cur_i = 0
        for exon in other.exons:
            while not self.exons[cur_i].contains(exon, tolerance):
                cur_i += 1
                if cur_i == len(self.exons):
                    return False
            cur_i += 1
        return True

    def is_cut_version(self, other, tolerance=0):
        """Check if other is made from cutting self
        at the edge of the alt loci

        :param other: alt-gene
        :param tolerance: number of base-pair tolerance
        :returns: whether other is cut
        :rtype: bool

        """
        # Check transcription_region

        f_tx_region = other.transcription_region.filter_to_main()
        if f_tx_region is None:
            return False
        f_exons = [e.filter_to_main() for e in other.exons]
        f_exons = [e for e in f_exons if e is not None]
        f_gene = Gene("f" + other.name, f_tx_region, f_exons)
        return self.contains(f_gene, tolerance)


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
        for interval in intervals:
            print(interval)
            
        n_starts = len([b for b in new_graph.blocks if not
                        new_graph.reverse_adj_list[b]])
        new_graph, final_trans = merge_flanks(intervals, final_trans,
                                              new_graph, name_translation)

        new_n_starts = len([b for b in new_graph.blocks if not
                            new_graph.reverse_adj_list[b]])
        assert new_n_starts == n_starts-1, line
    f.close()
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
    for gene in alt_genes:
        scores = []
        if gene.name not in main_dict:
            continue
        if gene in main_dict[gene.name]:
            gene_categories[11].append((gene, gene))
            continue
        for m_gene in main_dict[gene.name]:
            diffs = gene.start_and_end_diffs(m_gene)
            m_code = m_gene.transcription_region.diff(
                gene.transcription_region)
            exon_diffs = gene.exon_diffs(m_gene)
            if m_gene.contains(gene, 5):
                scores.append(9.1)
                continue
            elif m_gene.is_cut_version(gene):
                # print("------------------")
                # print(gene.to_file_line())
                # print(m_gene.to_file_line())
                scores.append(8.1)
                continue
            if exon_diffs is None:
                scores.append(-3)
                continue
            if all(c == "E" for c in m_code):
                if diffs[0] == 0 and diffs[1] == 0:
                    scores.append(10)
                    assert exon_diffs <= 40, "%s\n%s\n%s" % (exon_diffs, gene.to_file_line(), m_gene.to_file_line())
                    break
                elif max(diffs) < 5:
                    scores.append(9)
                    assert exon_diffs <= 2, "%s\n%s\n%s" % (exon_diffs, gene.to_file_line(), m_gene.to_file_line())
                    continue
                else:
                    scores.append(4)
                continue
            if m_code[0] == "E" and m_code[-1] == "E":
                if diffs[0] == 0 and diffs[1] == 0:
                    scores.append(8)
                    assert exon_diffs <= 2
                elif max(diffs) < 5:
                    scores.append(7)
                    assert exon_diffs <= 2
                else:
                    scores.append(3)
                continue
            if m_code[0] == "E":
                if diffs[0] == 0:
                    scores.append(6.1)
                    if not exon_diffs <= 2:
                        pass # print("%s\n%s\n%s" % (exon_diffs, gene.to_file_line(), m_gene.to_file_line()))
                elif diffs[0] < 5:
                    scores.append(5.1)
                    assert exon_diffs <= 2
                else:
                    scores.append(2.1)
                continue

            if m_code[-1] == "E":
                if diffs[1] == 0:
                    scores.append(6)
                    if not exon_diffs <= 2:
                        pass # print("%s\n%s\n%s" % (exon_diffs, gene.to_file_line(), m_gene.to_file_line()))
                elif diffs[1] < 5:
                    scores.append(5)
                    assert exon_diffs <= 2
                else:
                    scores.append(2)
                continue
            if all(c == "P" for c in m_code):
                if exon_diffs < 5:
                    scores.append(4.5)
                else:
                    scores.append(1)
                continue
            if "N" in m_code:
                scores.append(0)
                continue
            scores.append(-1)
        gene_scores = zip(scores, main_dict[gene.name])
        gene_scores = list(sorted(gene_scores, key=lambda x: x[0]))
        score = gene_scores[-1]
        if score[0] == 0:
            prev_blocks = graph.reverse_adj_list[gene.transcription_region.region_paths[0]]
            prev_critical = graph.find_previous_critical_block(gene.transcription_region.region_paths[0])
            if prev_blocks and False:
                print("______________________________")
                print(gene.length(), gene.to_file_line())
                print(prev_blocks, prev_critical)
                for m_gene in main_dict[gene.name]:
                    print(m_gene.length(), m_gene.to_file_line())
                    prev_blocks = graph.reverse_adj_list[m_gene.transcription_region.region_paths[0]]
                    prev_critical = graph.find_previous_critical_block(m_gene.transcription_region.region_paths[0])
                    print(prev_blocks, prev_critical)
                print(m_code)
        gene_categories[gene_scores[-1][0]].append((gene, gene_scores[-1][1]))
    return gene_categories


def find_exon_duplicates(genes, translation):
    """Find and count duplicate genes on flanks of
    alt loci

    :param genes: genes on original graph
    :param translation: translation object
    """
    # translated = GeneList.from_file("trans_genes").gene_list
    translated = [gene.translate(translation) for gene in genes]
    gene_list = GeneList(translated)
    gene_list.to_file("trans_genes")
    main_chr_dict = defaultdict(list)
    alt_dict = defaultdict(list)
    gene_categories = defaultdict(list)
    graph = translated[0].transcription_region.graph
    graph.critical_blocks = graph.find_all_critical_blocks()
    for gene, t_gene in zip(genes, translated):
        if "alt" in gene.chrom:
            alt_dict[gene.chrom.split("_")[0]].append(t_gene)
        else:
            main_chr_dict[gene.chrom].append(t_gene)
    print("N:", sum(len(v) for v in alt_dict.values()))
    s = 0
    for chrom, genes in alt_dict.items():
        new_things = find_unequal_sibling_genes(main_chr_dict[chrom], genes)
        for k, v in new_things.items():
            gene_categories[k].extend(v)
        continue
        for gene in genes:
            for main_gene in main_chr_dict[chrom]:
                if gene == main_gene:
                    print(main_gene.name, gene.name)
                    s += 1

    for category, v in gene_categories.items():
        print(category, len(v))

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


def grch38_graph_to_numeric(original_grch38_graph):
    # Convert to numeric graph
    num_ab = {}
    num_ba = {}
    i = 0
    graph = original_grch38_graph.copy()
    for b in graph.blocks:
        i += 1
        num_ab[b] = [Interval(0, graph.blocks[b].length(), [i])]
        num_ba[i] = [Interval(0, graph.blocks[b].length(), [b])]
    trans = Translation(num_ab, num_ba, graph=original_grch38_graph)
    trans.graph2 = trans.translate_subgraph(original_grch38_graph)
    original_grch38_graph._update_a_b_graph(num_ab, trans.graph2)
    original_grch38_graph._update_a_b_graph(num_ba, original_grch38_graph)
    graph = trans.graph2

    return graph, trans


def merge_alt_using_cigar(original_numeric_grch38_graph, trans, alt_id):
    """
    Uses a cigar string from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/
    and merges the alt locus.

    :param original_grch38_graph: Original numeric grch38 graph, built by calling the method create_initial_grch38_graph and then grch38_graph_to_numeric.
    :param trans: Translation object from original grch38 graph to numeric (returned by grch38_graph_to_numeric)
    :param alt_id: Alt id (e.g. chr2_KI270774v1_alt)
    :return: Returns the new graph and a translation object between the original grch38 graph and the new graph
    """

    # Find position of alt locus
    main_chr = ""
    main_start = 0
    main_end = 0
    alt_start = 0
    alt_end = 0

    try:
        with open("alt_alignments/%s.alignment" % alt_id) as f:
            d = f.read().split(",")
            cigar = d[-1]
            # All numbers are 1-indxed and inclusive end (?)
            # Convert to 0-indexed and exclusive end
            main_start = int(d[0]) - 1
            main_end = int(d[1]) - 1 + 1  # Was inclusive
            alt_start = int(d[2]) - 1
            alt_end = int(d[3]) - 1 + 1
    except:
        print("Could not open alignment file alt_alignments/%s.alignment" % alt_id)
        return trans, trans.graph2

    main_length = main_end - main_start
    alt_length = alt_end - alt_start


    """
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
    """
    main_chr = alt_id.split("_")[0]


    # Get sequences
    print("Fetching alt sequence")
    alt_seq = get_sequence_ucsc(alt_id, alt_start + 1, alt_end)  # uscs is 0 indexed and inclusive??
    print("Fetching main sequence")
    main_seq = get_sequence_ucsc(main_chr, main_start + 1, main_end)

    #print(alt_seq[0:10])
    #print(main_seq[0:10])
    ##print(alt_seq[-10:])
    #print(main_seq[-10:])
    #print(len(alt_seq))
    #print(alt_end - alt_start)
    #print(main_end - main_start)

    assert len(alt_seq) == (alt_end - alt_start)
    assert len(main_seq) == (main_end - main_start)

    #print(alt_seq[-10:])
    #print(main_seq[-10:])

    #print(main_chr, start, end, length)

    return _merge_alt_using_cigar(original_numeric_grch38_graph, trans, alt_id, cigar, alt_seq, main_seq, main_chr, main_start, main_end, alt_start, alt_end)


def _merge_alt_using_cigar(original_grch38_graph, trans, alt_id, cigar, alt_seq, main_seq, main_chr, main_pos_start, main_pos_end, alt_start, alt_end):
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

    length = main_pos_end - main_pos_start
    print("Merging %s" % alt_id)

    """
    Algorithm:
        go through cigar. For every cigar element, modify graph.
        If cigar is m,
    """
    #print(cigar)
    cigar = cigar.split(" ")
    #print("Alt start: %d" % alt_start)
    main_offset = main_pos_start
    alt_offset = alt_start
    print("alt offset: %d" % alt_offset)

    ab = {}
    ba = {}
    #print("=== Merge alt ===")
    graph = original_grch38_graph.copy()




    #print("== Numeric graph ===")
    #print(graph)
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


    #print("=== Graph before cigar ===")
    #print(graph)

    #print(cigar)
    for c in cigar:

        first = c[0]
        if first == "M" or first == "D" or first == "I":
            type = first
            n = int(c[1:])
        else:
            type = c[-1]
            n = int(c[:-1])

        #print("Type/len: %s/%d" % (type, n))

        if type == "M":
            print("Match %d" % (n))
            # If sequences are identical, then merge
            seq_on_alt = alt_seq[alt_offset-alt_start:alt_offset-alt_start+n]
            seq_on_main = main_seq[main_offset-main_pos_start:main_offset+n-main_pos_start]

            print(seq_on_alt)
            print(seq_on_main)
            """
            print(seq_on_alt[0:10])
            print(seq_on_alt[-10:])
            print(seq_on_main[0:10])
            print(seq_on_main[-10:])
            """
            prev = 0
            match = True
            for i in range(0, len(seq_on_alt)):


                if seq_on_alt[i] != seq_on_main[i] or i == n - 1:
                    if i <= prev:
                        continue
                    print("Merging from %d to %d" % (prev, i))
                    # Something is different, merge from prev to i
                    m_start = main_offset + prev
                    m_end = main_offset + prev + i
                    a_start = alt_offset + prev
                    a_end = alt_offset + prev + i
                    print("Merging alt %d,%d with main %d,%d" % (a_start, a_end, m_start, m_end))

                    intv1 = trans.translate(Interval(m_start, m_end, [main_chr], original_grch38_graph))
                    intv2 = trans.translate(Interval(a_start, a_end, [alt_id], original_grch38_graph))

                    print("Merging %s with %s" % (intv1, intv2))
                    graph, mtrans = graph.merge([intv1, intv2])

                    trans = trans + mtrans
                    trans.graph2 = graph

                    prev = i+1


            """
            if seq_on_alt == seq_on_main:
                print("Merging match")
                #print(trans)
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
            """
            main_offset += n
            alt_offset += n
            #print("Offsets alt/main: %d/%d" % (alt_offset, main_offset))

        elif type == "I":
            print("Insertion")
            alt_offset += n
        elif type == "D":
            print("Deletion")
            main_offset += n

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

    #print("=== Final trans ===")
    #print(trans)

    # Final translation should make final graph from original
    assert graph == trans.graph2 # (original_grch38_graph)

    return trans, graph
