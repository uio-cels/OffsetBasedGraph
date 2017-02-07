import unittest
from collections import defaultdict
from offsetbasedgraph.cigar_align import get_match_cigar, clean_cigar, align_cigar
from offsetbasedgraph import Block, Graph, Interval, Translation
from offsetbasedgraph.graphutils import convert_cigar_graph_to_text


class TestCigar(unittest.TestCase):
    def test_match_string(self):
        seq1 = "ACGTCGTATC"
        seq2 = "ACGTTTCATG"
        fasit = [("M", 4), ("V", 3), ("M", 2), ("V", 1)]
        answer = get_match_cigar(seq1, seq2)
        self.assertEqual(answer, fasit)

    def test_clean_cigar(self):
        cigar = " ".join(("D1", "M3", "I2", "M4"))
        seq_main = "".join(["A", "CGT", "CGTA"])
        seq_alt = "".join(["CAT", "GG", "CGGG"])
        fasit = [("D", 1), ("M", 1), ("V", 1), ("M", 1),
                 ("I", 2), ("M", 2), ("V", 2)]
        answer = clean_cigar(cigar, seq_alt, seq_main)
        self.assertEqual(answer, fasit)

    def test_align_cigar(self):
        cigar = [("D", 1), ("M", 1), ("V", 1), ("M", 1),
                 ("I", 2), ("M", 2), ("V", 2)]

        MAIN = -2
        ALT = -1
        blocks = {MAIN: Block(100),
                  ALT: Block(50)}

        graph = Graph(blocks, defaultdict(list))

        main_interval = Interval(40, 48, [MAIN])
        alt_interval = Interval(20, 29, [ALT])
        a_to_b = {MAIN: [Interval(0, 52, [0, 2, 3, 5, 6, 8, 10, 11])],
                  ALT: [Interval(0, 21, [1, 3, 4, 6, 7, 8,  9, 12])]}
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

        self.assertEqual(answer, fasit)

    def _check_name_graph(self, graph):
        codes = ["M", "A", "M",
                 "MA", "A", "M",
                 "MA", "A",
                 "MA", "A", "M", "M", "A"]

        lens = [40, 20, 1,
                1, 1, 1,
                1, 2,
                2, 2, 2, 52, 21]

        blocks = {codes[i] + str(i): Block(lens[i]) for i in range(13)}
        adj_list = {0: [2], 1: [3], 2: [3],
                    3: [4, 5], 4: [6], 5: [6],
                    6: [7, 8], 7: [8],
                    8: [9, 10], 9: [12], 10: [11]}

        name_adj_list = {codes[k] + str(k):
                         [codes[i]+str(i) for i in v]
                         for k, v in adj_list.items()}
        fasit_graph = Graph(blocks, name_adj_list)

        self.assertEqual(graph, fasit_graph)

    def test_name_translation_cigar(self):
        cigar = [("D", 1), ("M", 1), ("V", 1), ("M", 1),
                 ("I", 2), ("M", 2), ("V", 2)]

        MAIN = -2
        ALT = -1

        blocks = {"chr": Block(100),
                  "alt": Block(50)}

        name_dict = {"chr": MAIN,
                     "alt": ALT}

        graph = Graph(blocks, defaultdict(list))
        name_trans = Translation.make_name_translation(name_dict, graph)
        graph2 = name_trans.translate_subgraph(graph)
        name_trans.set_graph2(graph2)

        main_interval = Interval(40, 48, [MAIN])
        alt_interval = Interval(20, 29, [ALT])
        numeric_trans = align_cigar(cigar, main_interval, alt_interval, graph2)

        graph3 = numeric_trans.translate_subgraph(graph2)
        numeric_trans.set_graph2(graph3)

        final_graph, name_trans2 = convert_cigar_graph_to_text(
            graph3,
            name_trans,
            numeric_trans)

        codes = ["M", "A", "M", "MA", "A", "M", "MA",
                 "A", "MA", "A", "M", "M", "A"]

        a_to_b = {i: codes[i]+str(i) for i in range(13)}
        fasit = Translation.make_name_translation(a_to_b, graph3)

        self.assertEqual(name_trans2, fasit)

        full_translation = name_trans+numeric_trans+name_trans2

        main_rps = [0, 2, 3, 5, 6, 8, 10, 11]
        main_name_rps = [codes[i] + str(i) for i in main_rps]
        alt_rps = [1, 3, 4, 6, 7, 8, 9, 12]
        alt_name_rps = [codes[i] + str(i) for i in alt_rps]

        MAIN = "chr"
        ALT = "alt"

        a_to_b = {MAIN: [Interval(0, 52, main_name_rps)],
                  ALT: [Interval(0, 21, alt_name_rps)]}

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

        b_to_a = {codes[k] + str(k): v for k, v in b_to_a.items()}
        fasit = Translation(a_to_b, b_to_a)
        answer = full_translation
        self.assertEqual(answer, fasit)
        final_graph = full_translation.translate_subgraph(graph)
        self._check_name_graph(final_graph)

if __name__ == "__main__":
    unittest.main()
