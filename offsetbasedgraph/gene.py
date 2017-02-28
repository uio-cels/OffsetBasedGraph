import csv
import pickle
from collections import defaultdict
from .interval import Interval, Position


class GeneList(object):
    def __init__(self, gene_list):
        assert isinstance(gene_list, list)
        self.gene_list = gene_list
        self.lookup = defaultdict(list)
        for gene in self.gene_list:
            self.lookup[gene.name].append(gene)

    def to_file(self, file_name):
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    def translate(self, T):
        t_gene_list = [gene.translate(T) for gene in self.gene_list]
        return GeneList(t_gene_list)

    @staticmethod
    def from_pickle(file_name):
        with open("%s" % file_name, "rb") as f:
            o = pickle.loads(f.read())
            if isinstance(o, list):
                o = GeneList(o)
            return o

    @classmethod
    def from_file(cls, file_name):
        with open(file_name) as f:
            reader = csv.DictReader(f, delimiter="\t")
            gene_list = [Gene.from_dict(d) for d in
                         reader]
        return cls(gene_list)

    def __eq__(self, other):
        if len(self.gene_list) != len(other.gene_list):
            return False

        for g in self.gene_list:
            if g not in other.gene_list:
                return False

        for g in other.gene_list:
            if g not in self.gene_list:
                return False

        return True


class MultiPathGene(object):
    def __init__(self, name, multipath_interval):
        self.name = name
        self.interval = multipath_interval


class GeneBase(object):
    def __init__(self, name, transcription_region, exons,
                 coding_region, strand):
        self.name = name
        self.transcription_region = transcription_region
        self.coding_region = coding_region
        self.exons = exons
        self.strand = strand
        self.chrom = transcription_region.region_paths[0]
        self.graph = self.transcription_region.graph
        self.transcript_length = sum(exon.length() for exon in exons)

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

        return all(e1 == e2 for e1, e2 in zip(self.exons, other.exons))


class FuzzyGene(GeneBase):
    def __str__(self):
        return "%s\t%s" % (self.name, str(self.transcription_region))

    __repr__ = __str__


class Gene(GeneBase):

    def copy(self):
        coding_region = None
        if self.coding_region is not None:
            coding_region = self.coding_region.copy()
        return Gene(self.name, self.transcription_region.copy(),
                    [ex.copy() for ex in self.exons],
                    coding_region,
                    self.strand,
                    )

    ###### Should be elswere
    def multiple_alt_loci(self):
        """
        Returns true if gene stretches over multiple alt loci.
        NB: Works only when gene is represented on a graph where
        alt loci blocks have names with alt loci in them
        """
        from .graph import Graph
        unique_loci = []
        for b in self.transcription_region.region_paths:
            if Graph.block_origin(b) == "alt":
                unique_loci.append(b)

        if len(set(unique_loci)) > 1:
            return True
        else:
            return False

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

        coding_region = Interval(
            Position(chrom, int(attr_dict["cdsStart"])),
            Position(chrom, int(attr_dict["cdsEnd"])), [chrom]
            )

        exon_starts = [int(i) for i in attr_dict["exonStarts"].split(",")[:-1]]
        exon_ends = [int(i) for i in attr_dict["exonEnds"].split(",")[:-1]]
        exons = [Interval(start, end, [chrom]) for start, end in
                 zip(exon_starts, exon_ends)]
        strand = attr_dict["strand"]

        return cls(attr_dict["name"], transcription_region, exons,
                   coding_region, strand)

    def translate(self, T):
        """Translate transcription_region and exons and return
        new gene

        :param T: Translation object
        :returns: Translated gene
        :rtype: Gene

        """

        t_transcription_region = T.translate(self.transcription_region)
        t_coding_region = self.coding_region  # T.translate(self.coding_region)
        t_exons = [T.translate(exon) for exon in self.exons]
        return Gene(self. name, t_transcription_region, t_exons,
                    t_coding_region, self.strand)

    def get_transcript_length(self):
        """Return combined length of exons

        :returns: transcript length
        :rtype: int

        """
        return sum(exon.length() for exon in self.exons)

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

        return all(e1 == e2 for e1, e2 in zip(self.exons, other.exons))
#
#        for exon in self.exons:
#            if exon not in other.exons:
#                return False
#
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
                s_diff = exon1.start_position.offset - exon2.start_position.offset
                e_diff = exon1.end_position.offset-exon2.end_position.offset
                diff_sum += abs(s_diff)+abs(e_diff)
            else:
                diff_sum += abs(exon1.length()-exon2.length())

        return diff_sum/len(self.exons)

    def approxEquals(self, other, tolerance=0):
        if not len(self.exons) == len(other.exons):
            return False
        my_regions = [self.transcription_region] + self.exons
        other_regions = [other.transcription_region] + other.exons
        return all(my_reg.approx_equals(other_reg) for
                   my_reg, other_reg in zip(my_regions, other_regions))

    def contains(self, other, tolerance=0):
        """Checks if self contains other. Allows
        mismatch in start and end coordinated up to 
        tolerance

        :param other: Other gene
        :param tolerance: bps allowed discrepancy
        :returns: Wheter self contains other
        :rtype: bool

        """
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

    def get_contained_pairs(self, other, tolerance=0):
        pairs = []
        for exon in other.exons:
            containers = [s_exon for s_exon in self.exons
                          if s_exon.contains(exon, tolerance)]
            if not containers:
                continue
            assert len(containers) == 1
            pairs.append((containers[0], exon))
        return pairs

    def partition_region_paths(self, region_paths=[]):
        """Partition the gene into one part for each
        region paths. Splits the transcription_region and
        the exons, but NB! not the coding region

        :returns: list of gene-parts
        :rtype: list(Gene)

        """
        rps = self.transcription_region.region_paths
        if len(rps) == 0:
            return [self]

        transcription_regions = self.transcription_regions.partition_region_paths()
        exon_partitions = chain([exon.partition_region_paths()
                                 for exon in self.exons])

        rps = [tr.region_paths[0] for tr in transcription_regions]
        exons_per_rp = {rp: [e for e in exon_partitions
                             if e.region_paths[0] == rp]
                        for rp in rps}
        gene_partitions = []
        for rp, tr in zip(rps, transcription_regions):
            gene_partitions.append(
                Gene(self.name, tr, exons_per_rp[rp],
                     self.coding_region, self.strand))

        return gene_partitions

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
