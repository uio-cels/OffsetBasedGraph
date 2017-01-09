import unittest
from offsetbasedgraph import Interval, Position, Graph, Translation, Block

class testExamples(unittest.TestCase):

    def test_overlapping_alt_loci(self):
        chrom_file = "examples/chrom.sizes.test"
        alt_loci = "examples/alt_loci_test"

        import examples.gene_experiment as g
        graph = g.create_initial_grch38_graph(chrom_file)
        numeric_graph, name_translation = g.convert_to_numeric_graph(graph)

        self.assertEqual(len(graph.blocks), 3)
        self.assertEqual(len(graph.adj_list), 0)

        new_numeric_graph, numeric_translation = \
                g.connect_without_flanks(numeric_graph, alt_loci, name_translation)

        # Find first region path
        ngraph = new_numeric_graph
        first = None
        for rp in ngraph.blocks:
            for other_rp in ngraph.blocks:
                if other_rp != rp and rp in ngraph.adj_list[other_rp]:
                    continue
                first = rp
                break

        print(ngraph)

        assert first is not None
        print("First" , first)

        self.assertEqual(len(ngraph.adj_list[first]), 1)

        next = ngraph.adj_list[first][0]
        self.assertEqual(len(ngraph.adj_list[next]), 2)


if __name__ == "__main__":
    unittest.main()

