class LinearInterval(object):
    """An interval on a linear reference"""

    def __init__(self, genome_id, chromosome, start, end, strand="+"):
        self.genome_id = genome_id
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.gene_name = ""  # Used in some cases

    def length(self):
        return self.end-self.start

    def contains(self, other):
        if not (self.genome_id == other.genome_id):
            return False
        if not self.chromosome == other.chromosome:
            return False
        if self.start > other.start or self.end < other.end:
            return False
        return True

    def contains_position(self, chrom, offset):
        if not chrom == self.chromosome:
            return False
        return offset >= self.start and offset < self.end

    def distance_to_position(self, chrom, offset):
        assert chrom == self.chromosome
        if self.contains_position(chrom, offset):
            return 0
        return min(abs(self.start-offset), abs(self.end-offset))

    def intersects(self, other):
        if not (self.genome_id == other.genome_id):
            return False
        if not self.chromosome == other.chromosome:
            return False
        if (self.start < other.end and self.end > other.start) or \
                (other.start < self.end and other.end > self.start):
            return True

        return False

    def split(self, offsets):
        assert all(offset < self.length() for offset in offsets)
        starts = [self.start]+[self.start+offset for offset in offsets]
        ends = [self.start+offset for offset in offsets] + [self.end]

        return [LinearInterval(self.genome_id, self.chromosome,
                               start, end, self.strand)
                for start, end in zip(starts, ends)]

    def __add__(self, other):
        """TODO: check that they are continuous"""
        assert self.chromosome == other.chromosome
        assert self.genome_id == other.genome_id
        assert self.strand == other.strand
        start = min(self.start, other.start)
        end = max(self.end, other.end)
        name = self.gene_name + other.gene_name
        interval = LinearInterval(self.genome_id, self.chromosome,
                                  start, end, self.strand)
        interval.gene_name = name
        return interval

    def __eq__(self, other):
        if not self.chromosome == other.chromosome:
            return False
        if not self.start == other.start:
            return False
        if not self.end == other.end:
            return False

        return True

    def __repr__(self):
        return "Lin seg in species %s, %s [%d, %d] %s" % \
               (self.genome_id, self.chromosome,
                self.start, self.end, self.strand)

    def label(self):
        return "%s" % (self.chromosome)

    def __str__(self):
        return "%s from %d % d on chr %s" % (
            self.label(), self.start, self.end, self.chromosome)
