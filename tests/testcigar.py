import unittest
from collections import defaultdict
from offsetbasedgraph.cigar_align import get_match_cigar, clean_cigar, align_cigar
from offsetbasedgraph import Block, Graph, Interval, Translation


class TestCigar(unittest.TestCase):
    def test_match_string(self):
        seq1 = "ACGTCGTATC"
        seq2 = "ACGTTTCATG"
        fasit = [("M", 4), ("V", 3), ("M", 2), ("V", 1)]
        answer = get_match_cigar(seq1, seq2)
        self.assertEqual(answer, fasit)

    def test_clean_cigar(self):
        cigar = ("D1", "M3", "I2", "M4")
        seq_main = "".join(["A", "CGT", "CGTA"])
        seq_alt = "".join(["CAT", "GG", "CGGG"])
        fasit = [("D", 1), ("M", 1), ("V", 1), ("M", 1),
                 ("I", 2), ("M", 2), ("V", 2)]
        answer = clean_cigar(cigar, seq_alt, seq_main, 0, 0)
        self.assertEqual(answer, fasit)

    def test_align_cigar(self):
        cigar = [("D", 1), ("M", 1), ("V", 1), ("M", 1),
                 ("I", 2), ("M", 2), ("V", 2)]
        blocks = {100: Block(100),
                  101: Block(50)}

        graph = Graph(blocks, defaultdict(list))
        main_interval = Interval(40, 48, [100])
        alt_interval = Interval(20, 29, [101])
        MAIN = 100
        ALT = 101
        a_to_b = {100: [Interval(0, 52, [0,2,3,5,6,  8,10,11])],
                  101: [Interval(0, 21, [1,  3,4,6,7,8, 9,12])]}
        b_to_a = {0: [Interval(0, 40, [MAIN])],
                  1: [Interval(0, 20, [ALT])],
                  2: [Interval(40, 41, [MAIN])],
                  3: [Interval(41, 42, [MAIN]),
                      Interval(20, 21, [ALT])],
                  4: [Interval(21, 22, [ALT])],
                  5: [Interval(42, 43, [MAIN])],
                  6: [Interval(43, 44, [MAIN]),
                      Interval(22, 23, [ALT])],
                  7: [Interval(23, 25, [ALT])],
                  8: [Interval(44, 46, [MAIN]),
                      Interval(25, 27, [ALT])],
                  9: [Interval(27, 29, [ALT])],
                  10: [Interval(46, 48, [MAIN])],
                  11: [Interval(48, 100, [MAIN])],
                  12: [Interval(29, 50, [ALT])]}
        fasit = Translation(a_to_b, b_to_a)
        answer = align_cigar(cigar, main_interval, alt_interval, graph)
        for k in fasit._b_to_a:
            print(fasit._b_to_a[k])
            print(answer._b_to_a[k])

        self.assertEqual(answer, fasit)

if __name__ == "__main__":
    unittest.main()
