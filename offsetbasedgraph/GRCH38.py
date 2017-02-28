from offsetbasedgraph import Interval, Position
from offsetbasedgraph.sequences import get_sequence_ucsc


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
        flank_intervals = get_flanks(
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


def get_seq_from_file(filename):
    """ Read sequence from fasta formated file
    TODO: Should be done by fasta reader"""

    return "".join(
        [line.strip() for line in
         open(filename, "r").readlines()[1:]])


def get_start_flank_length(seq1, seq2):
    start_flank_length = 0
    for i in range(min(len(seq1), len(seq2))):
        start_flank_length = i
        if seq1[i] != seq2[i]:
            break

    return start_flank_length


def get_stop_flank_length(seq1, seq2):
    stop_flank_length = 0
    for i in range(min(len(seq1), len(seq2))):
        stop_flank_length = i
        if seq1[-i-1] != seq2[-i-1]:
            break

    return stop_flank_length


def get_identical_flanks(seq1, seq2):
    start_flank_length = get_start_flank_length(seq1, seq2)
    stop_flank_length = get_stop_flank_length(seq1, seq2)
    return start_flank_length, stop_flank_length


def get_split_list(start, start_flank_length, stop_flank_length, end):
    mid_start = start + start_flank_length
    mid_end = end - stop_flank_length
    return [start, mid_start, mid_end, end]


def get_intervals_from_split_list(split_list, region_name):
    intervals = []
    for start, stop in zip(split_list[:-1], split_list[1:]):
        intervals.append(
            Interval(Position(region_name, start),
                     Position(region_name, stop)))

    return intervals


def get_flanks(region_name, length, main_chromosome, chrom_start, chrom_end):
    """
    Get linear intervals for
    start_flank, varying_region, end_flank
    for alt_loci and the corresponing main region
    """
    # region_name = alt_info["name"]
    alt_seq = get_sequence_ucsc(
        region_name, 1, length)


    consensus_seq = get_sequence_ucsc(
        main_chromosome, chrom_start+1, chrom_end)

    """
    print(alt_seq[0:10])
    print(alt_seq[-10:])
    print(consensus_seq[0:10])
    print(consensus_seq[-10:])
    """

    start_flank_length, stop_flank_length = get_identical_flanks(
        alt_seq, consensus_seq)

    alt_coords = get_split_list(0, start_flank_length, stop_flank_length,
                                length)
    main_coords = get_split_list(chrom_start, start_flank_length,
                                 stop_flank_length, chrom_end)

    alt_intervals = get_intervals_from_split_list(alt_coords, region_name)
    main_intervals = get_intervals_from_split_list(main_coords,
                                                   main_chromosome)
    """
    if alt_intervals[2].length() == 0:
        alt_intervals[2].start_position.offset -= 1
        alt_intervals[2].end_position.offset -= 1
        main_intervals[2].start_position.offset -= 1
        main_intervals[2].end_position.offset -= 1
    """

    return [main_intervals[0], alt_intervals[0],
            main_intervals[2], alt_intervals[2]]