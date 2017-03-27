import csv
import pickle
from collections import defaultdict
from .interval import Interval, Position


class GeneList(object):
    """ Container for collection of genes"""

    def __init__(self, gene_list):
        assert isinstance(gene_list, list)
        self.gene_list = gene_list
        self.lookup = defaultdict(list)
        for gene in self.gene_list:
            self.lookup[gene.name].append(gene)

    def to_file(self, file_name):
        """Write gene list to pickle"""
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def from_pickle(file_name):
        """Read gene list from pickle"""
        with open("%s" % file_name, "rb") as f:
            o = pickle.loads(f.read())
            if isinstance(o, list):
                o = GeneList(o)
            return o

    @classmethod
    def from_file(cls, file_name):
        """Read gene list from csv file
        
        :param file_name: csv file name
        :returns: list of genes
        :rtype: GeneList
        """
        with open(file_name) as f:
            reader = csv.DictReader(f, delimiter="\t")
            gene_list = [Gene.from_dict(d) for d in
                         reader]
        return cls(gene_list)

    def to_tsv_file(self, file_name):
        header = "# Name\t TranscriptionRegion \t CodingRegion \t Exons \CodingStatus"
        lines = [GeneIO(gene).to_file_line() for
                 gene in self.gene_list]
        with open(file_name, "w") as f:
            f.write(header + "\n" + "\n".join(lines))

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

    def translate(self, translation):
        translated_genes = [gene.translate(translation)
                            for gene in self.gene_list]
        return GeneList(translated_genes)

    def __str__(self):
        lines = [GeneIO(gene).to_file_line() for gene in self.gene_list]
        return "GeneList:\n" + "\n".join(lines)


class MultiPathGene(object):
    def __init__(self, name, multipath_interval):
        self.name = name
        self.interval = multipath_interval


class GeneBase(object):
    def __init__(self, name, transcription_region, exons,
                 coding_region, strand, cds_status):
        self.name = name
        self.transcription_region = transcription_region
        self.coding_region = coding_region
        self.exons = exons
        self.strand = strand
        self.chrom = transcription_region.region_paths[0]
        self.graph = self.transcription_region.graph
        self.transcript_length = sum(exon.length() for exon in exons)
        self.cds_status = cds_status

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
    """ Gene representation that includes in adddition to the specified
    intervals all intervals that deviate with less than a threshold from
    those.
    """

    def __str__(self):
        return "%s\t%s" % (self.name, str(self.transcription_region))

    __repr__ = __str__


class Gene(GeneBase):

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
        cds_status = (attr_dict["cdsStartStat"], attr_dict["cdsEndStat"])

        return cls(attr_dict["name"], transcription_region, exons,
                   coding_region, strand, cds_status)

    def translate(self, T):
        """Translate transcription_region and exons and return
        new gene

        :param T: Translation object
        :returns: Translated gene
        :rtype: Gene

        """

        t_transcription_region = T.translate(self.transcription_region)
        t_coding_region = T.translate(self.coding_region)
        t_exons = [T.translate(exon) for exon in self.exons]
        return Gene(self. name, t_transcription_region, t_exons,
                    t_coding_region, self.strand, self.cds_status)

    def get_transcript_length(self):
        """Return combined length of exons

        :returns: transcript length
        :rtype: int

        """
        return sum(exon.length() for exon in self.exons)

    def length(self):
        """Return length of transcription region
        :rtype: int
        """
        return self.transcription_region.length()


    def approxEquals(self, other, tolerance=0):
        """Check if other is approximately equal to self.
        Difference between start and end positions of 
        transcription region as well as all exons must be
        less than tolerance

        :param other: Gene
        :param tolerance: allowed differnce (int)
        :rtype: bool

        """
        if not len(self.exons) == len(other.exons):
            return False
        my_regions = [self.transcription_region] + self.exons
        other_regions = [other.transcription_region] + other.exons
        return all(my_reg.approx_equals(other_reg) for
                   my_reg, other_reg in zip(my_regions, other_regions))

    def get_contained_pairs(self, other, tolerance=0):
        """Get list of all pairs (my, his) of exons in other
        that are contained in an exon in self

        :param other: Gene
        :param tolerance: tolerance for contains (int)
        :returns: list of pairs of (container, containee) exons
        :rtype: list( (Interval, Interval)...)

        """
        pairs = []
        for exon in other.exons:
            containers = [s_exon for s_exon in self.exons
                          if s_exon.contains(exon, tolerance)]
            if not containers:
                continue
            assert len(containers) == 1
            pairs.append((containers[0], exon))
        return pairs

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
        return True

    def copy(self):
        coding_region = None
        if self.coding_region is not None:
            coding_region = self.coding_region.copy()
        return Gene(self.name, self.transcription_region.copy(),
                    [ex.copy() for ex in self.exons],
                    coding_region,
                    self.strand,
                    )


class GeneIO(object):
    """Unfinished gene writer/reader"""
    def __init__(self, gene):
        self.gene = gene

    def to_file_line(self):
        cur_rp = self.gene.transcription_region.region_paths[0]
        exon_strings = []
        for exon in self.gene.exons:
            rps = [rp for rp in exon.region_paths if rp != cur_rp]
            offsets = [exon.start_position.offset, exon.end_position.offset]
            exon_string = ",".join(str(e) for e in rps+offsets)
            exon_strings.append(exon_string)
            cur_rp = exon.region_paths[-1]
        transcript_string = ",".join(
            self.gene.transcription_region.region_paths+[
                str(self.gene.transcription_region.start_position.offset),
                str(self.gene.transcription_region.end_position.offset)])

        coding_region_string = ",".join(
            self.gene.coding_region.region_paths+[
                str(self.gene.coding_region.start_position.offset),
                str(self.gene.coding_region.end_position.offset)])

        cds_status_string = " ".join(self.gene.cds_status)
        return "\t".join([self.gene.name,
                          transcript_string,
                          coding_region_string,
                          " ".join(exon_strings),
                          cds_status_string])

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
        transcription_region = cls.parse_transcript_region(
            transcript_string)
        exons = cls.parse_exons(exons_string)
