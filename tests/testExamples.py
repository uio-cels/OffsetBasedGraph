import unittest
from offsetbasedgraph import Interval, Position, Graph, Translation, Block

class testExamples(unittest.TestCase):

    def test_overlapping_alt_loci(self):
        chrom_file = "examples/chrom.sizes.test"
        alt_loci = "examples/alt_loci_test"

        import examples.gene_experiment as g
        graph = g.create_initial_grch38_graph(chrom_file)
        print("Original graph")
        print(graph)
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
                9: [4],
                8: [6]
            }
        )

        self.assertTrue(correct_graph_structure.has_identical_structure(new_numeric_graph))

        """
        # Find first region path
        ngraph = new_numeric_graph
        first = None
        for rp in ngraph.blocks:
            for other_rp in ngraph.blocks:
                if other_rp != rp and rp in ngraph.adj_list[other_rp]:
                    continue
                first = rp
                break

        self.assertEqual(len(ngraph.blocks), 9)

        print(ngraph)

        assert first is not None
        print("First" , first)

        self.assertEqual(len(ngraph.adj_list[first]), 2)
        next_main = ngraph.adj_list[first][0]
        next_alt = ngraph.adj_list[first][1]

        if len(ngraph.adj_list[next_main]) == 1:
            # swap
            next_main = ngraph.adj_list[first][1]
            next_alt = ngraph.adj_list[first][0]

        self.assertEqual(len(ngraph.adj_list[next_main]), 2)
        self.assertEqual(len(ngraph.adj_list[next_alt]), 1)

        # Check path from next main
        self.assertEqual()


        next = ngraph.adj_list[first][0]
        self.assertEqual(len(ngraph.adj_list[next]), 2)
        """

if __name__ == "__main__":
    unittest.main()

