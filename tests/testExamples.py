import unittest
from offsetbasedgraph import Interval, Position, Graph, Translation, Block

class testExamples(unittest.TestCase):

    def _test_overlapping_alt_loci(self):
        chrom_file = "examples/chrom.sizes.test"
        alt_loci = "examples/alt_loci_test"

        import examples.gene_experiment as g
        graph = g.create_initial_grch38_graph(chrom_file)

        numeric_graph, name_translation = g.convert_to_numeric_graph(graph)

        self.assertEqual(len(graph.blocks), 3)
        self.assertEqual(len(graph.adj_list), 0)

        new_numeric_graph, numeric_translation = \
                g.connect_without_flanks(numeric_graph, alt_loci, name_translation)


        correct_graph_structure = Graph(
            {
                1: Block(1),
                2: Block(1),
                3: Block(1),
                4: Block(1),
                5: Block(1),
                6: Block(1),
                7: Block(1),
                8: Block(1),
                9: Block(1),
            },
            {
                1: [2, 8],
                2: [3, 9],
                3: [4],
                4: [5],
                5: [6],
                6: [7],
                9: [5],
                8: [6]
            }
        )



        self.assertTrue(correct_graph_structure.has_identical_structure(new_numeric_graph))

    def _test_merge_flanks_2x(self):

        # Merge flanks of two alt loci
        g = Graph({1: Block(8), 2: Block(3), 3: Block(5)}, {})

        # MERGE FIRST ALT LOCUS
        intervals = [
                        Interval(2, 3, [1], g),
                        Interval(0, 1, [2], g),
                        Interval(4, 5, [1], g),
                        Interval(2, 3, [2], g)
                    ]

        import examples.gene_experiment as e
        final_trans = Translation({}, {}, graph=g)
        final_trans.graph2 = g
        name_trans = Translation({}, {}, graph=g)
        new_graph, trans = e.merge_flanks(intervals, final_trans, g, name_trans)

        correct_structure = Graph(
            {
                1: Block(1),
                2: Block(1),
                3: Block(1),
                4: Block(1),
                5: Block(1),
                6: Block(1),
                7: Block(1),
            },
            {
                1: [2],
                2: [3, 6],
                3: [4],
                6: [4],
                4: [5]
            }
        )
        self.assertTrue(new_graph.has_identical_structure(correct_structure))

        # MERGE SECOND ALT LOCUS
        intervals = [
            Interval(1, 3, [1], g),
            Interval(0, 2, [3], g),
            Interval(5, 6, [1], g),
            Interval(4, 5, [3], g)
        ]
        new_graph, trans = e.merge_flanks(intervals, trans, new_graph, name_trans)

        correct_structure = Graph(
            {
                1: Block(1),
                2: Block(1),
                3: Block(1),
                4: Block(1),
                5: Block(1),
                6: Block(1),
                7: Block(1),
                8: Block(1),
                9: Block(1),
            },
            {
                1: [2],
                2: [3],
                3: [8, 4, 9],
                4: [5],
                5: [6],
                6: [7],
                8: [5],
                9: [6]
            }
        )
        self.assertTrue(new_graph.has_identical_structure(correct_structure))

    def _test_merge_flanks_2x_case2(self):

        # Merge flanks of two alt loci
        g = Graph({1: Block(8), 2: Block(3), 3: Block(5)}, {})

        # MERGE FIRST ALT LOCUS
        intervals = [
                        Interval(2, 3, [1], g),
                        Interval(0, 1, [2], g),
                        Interval(4, 5, [1], g),
                        Interval(2, 3, [2], g)
                    ]

        import examples.gene_experiment as e
        final_trans = Translation({}, {}, graph=g)
        final_trans.graph2 = g
        name_trans = Translation({}, {}, graph=g)
        new_graph, trans = e.merge_flanks(intervals, final_trans, g, name_trans)


        # MERGE SECOND ALT LOCUS
        intervals = [
            Interval(2, 4, [1], g),
            Interval(0, 2, [3], g),
            Interval(5, 6, [1], g),
            Interval(4, 5, [3], g)
        ]
        new_graph, trans = e.merge_flanks(intervals, trans, new_graph, name_trans)

        correct_structure = Graph(
            {
                1: Block(1),
                2: Block(1),
                3: Block(1),
                4: Block(1),
                5: Block(1),
                6: Block(1),
                7: Block(1),
                8: Block(1)
            },
            {
                1: [2],
                2: [3, 7],
                3: [8, 4],
                4: [5],
                5: [6],
                7: [4],
                8: [5]
            }
        )

        print(new_graph)
        self.assertTrue(new_graph.has_identical_structure(correct_structure))

    def test_merge_alt_using_cigar(self):

        # Case 1
        graph = Graph({
                "chr1": Block(30),
                "chr1_test_alt": Block(10)
            },
            {});

        from offsetbasedgraph.graphutils import  grch38_graph_to_numeric
        graph, trans = grch38_graph_to_numeric(graph)

        start = 10
        end = 19

        alt_seq = "CCCTGGGAAA"
        main_seq = "CCCGGGTAAA"
        cigar = "M3 I1 M3 1D M3"
        from offsetbasedgraph.graphutils import _merge_alt_using_cigar
        trans, new_graph = _merge_alt_using_cigar(graph, trans,
                               "chr1_test_alt",
                               cigar,
                               alt_seq,
                               main_seq,
                               "chr1",
                               start,
                               end,
                               0,
                               10)

        correct_structure = Graph(
            {
                1: Block(1),
                2: Block(1),
                3: Block(1),
                4: Block(1),
                5: Block(1),
                6: Block(1),
                7: Block(1),
             },
            {
                1: [2],
                2: [3, 6],
                3: [4, 7],
                4: [5],
                6: [3],
                7: [4]
            }
        )
        self.assertTrue(new_graph.has_identical_structure(correct_structure))







if __name__ == "__main__":
    unittest.main()

