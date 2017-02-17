from offsetbasedgraph import Interval
from genutils import flanks


class AltLocus(object):

    def __init__(self, name, length, chrom, start, end, find_flanks=True):
        self.name = name
        self.length = length
        self.chrom = chrom
        self.start = start
        self.end = end
        self.alt_interval = Interval(0, length, [name])
        self.main_interval = Interval(start, end, [chrom])
        if find_flanks:
            self.find_flanks()

    @classmethod
    def from_file_line(cls, line):
        """Read a alt locus entry from a file line.
        Format:
        name\t chromsome \t start\t end\t length\t

        :param line: file line (str)
        :returns: alt locus object
        :rtype: AltLocus

        """
        l = line.split()
        name = l[0]
        chrom = l[1]
        start = int(l[2])
        end = int(l[3])
        length = int(l[4])
        return cls(name, length, chrom, start, end)

    def find_flanks(self):
        flank_intervals = flanks.get_flanks(
            self.name, self.length,
            self.chrom,
            self.start-1, self.end)
        self.start_flank = flank_intervals[1]
        self.end_flank = flank_intervals[3]
        self.varying_region = Interval(self.start_flank.end_position,
                                       self.end_flank.start_position)

        self.main_start_flank = flank_intervals[0]
        self.main_end_flank = flank_intervals[2]
        self.main_varying_region = Interval(
            self.main_start_flank.end_position,
            self.main_end_flank.start_position)

    def __str__(self):
        return "AltLoci(%s:%s)" % (self.name, self.length)

    def __repr__(self):
        return self.__str__()


class AltLoci(object):
    def __init__(self, alt_loci):
        self.alt_loci = alt_loci
        self.lookup = {alt_locus.name: alt_locus for alt_locus in alt_loci}

    @classmethod
    def from_file(cls, filename, filter_alt_loci=[]):
        with open(filename) as f:
            alt_loci = []
            for line in f.readlines():
                if line.startswith("#"):
                    continue

                if len(filter_alt_loci) > 0 and line.split()[0] not in filter_alt_loci:
                    continue

                alt_loci.append(AltLocus.from_file_line(line))

        return cls(alt_loci)


class Chromosome(object):

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.interval = Interval(0, length, [name])
