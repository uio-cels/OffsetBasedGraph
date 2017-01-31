from .interval import Interval, Position
from .graph import Graph, Block
from .translation import Translation
import csv
import pickle
from genutils import flanks


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
        print("Connecting %s" % (alt_locus_id))
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
    graph.critical_blocks = graph.find_all_critical_blocks()
    main_dict = defaultdict(list)
    gene_categories = defaultdict(list)
    for gene in main_genes:
        main_dict[gene.name].append(gene)
    for gene in alt_genes:
        scores = []
        if gene.name in main_dict and gene not in main_dict[gene.name]:
            for m_gene in main_dict[gene.name]:
                diffs = gene.start_and_end_diffs(m_gene)
                m_code = m_gene.transcription_region.diff(
                    gene.transcription_region)
                if all(c == "E" for c in m_code):
                    if diffs[0] == 0 and diffs[1] == 0:
                        print(gene.to_file_line())
                        print(m_gene.to_file_line())
                        scores.append(10)
                        break
                    elif max(diffs) < 5:
                        scores.append(9)
                        continue
                    else:
                        scores.append(4)
                    continue
                if m_code[0] == "E" and m_code[-1] == "E":
                    if diffs[0] == 0 and diffs[1] == 0:
                        scores.append(8)
                    elif max(diffs) < 5:
                        scores.append(7)
                    else:
                        scores.append(3)
                    continue
                if m_code[0] == "E":
                    if diffs[0] == 0:
                        scores.append(6.1)
                    elif diffs[0] < 5:
                        scores.append(5.1)
                    else:
                        scores.append(2.1)
                    continue

                if m_code[-1] == "E":
                    if diffs[1] == 0:
                        scores.append(6)
                    elif diffs[1] < 5:
                        scores.append(5)
                    else:
                        scores.append(2)
                    continue
                if all(c == "P" for c in m_code):
                    if abs(gene.length()-m_gene.length()) < 5:
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
            gene_categories[gene_scores[-1][0]].append((gene, gene_scores[-1][1]))
    return gene_categories


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
    gene_categories = defaultdict(list)
    for gene, t_gene in zip(genes, translated):
        if "alt" in gene.chrom:
            alt_dict[gene.chrom.split("_")[0]].append(t_gene)
        else:
            main_chr_dict[gene.chrom].append(t_gene)

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

    # print(s, "/", sum(len(v) for v in alt_dict.values()))
