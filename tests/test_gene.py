import unittest
from offsetbasedgraph import Graph, Block, Interval, Position

from offsetbasedgraph.gene import Gene, GeneList, FuzzyGene, MultiPathGene

class TestGene(unittest.TestCase):

    def test_gene_list_from_file(self):
        test_file = "tests/genes_chr1_KI270760v1_alt.txt"
        gene_list = GeneList.from_file(test_file)

        self.assertEqual(len(gene_list.gene_list), 3)

        for g in gene_list.gene_list:
            self.assertEqual(g.chrom, "chr1_KI270760v1_alt")

    def test_genelist_pickle(self):
        test_file = "tests/genes_chr1_KI270760v1_alt.txt"
        gene_list = GeneList.from_file(test_file)

        gene_list.to_file("tmp_gene_list")
        self.assertEqual(gene_list, GeneList.from_pickle("tmp_gene_list"))

    def test_gene_list_equal(self):
        test_file1 = "tests/genes_chr1_KI270760v1_alt.txt"
        test_file2 = "tests/genes_chr1_GL383518v1_alt.txt"
        test_file3 = "tests/genes_chr1_KI270760v1_alt_with_error.txt"
        gene_list1 = GeneList.from_file(test_file1)
        gene_list2 = GeneList.from_file(test_file2)
        gene_list3 = GeneList.from_file(test_file3)

        self.assertEqual(gene_list1, gene_list1)
        self.assertEqual(gene_list2, gene_list2)
        self.assertTrue(gene_list2 != gene_list1)
        self.assertTrue(gene_list1 != gene_list2)
        self.assertTrue(gene_list3 != gene_list1)

    def test_gene_equal(self):
        test_file1 = "tests/genes_chr1_KI270760v1_alt.txt"
        test_file2 = "tests/genes_chr1_KI270760v1_alt_with_error.txt"

        gene_list1 = GeneList.from_file(test_file1).gene_list
        gene_list2 = GeneList.from_file(test_file2).gene_list

        self.assertEqual(gene_list1[0], gene_list2[0])
        self.assertTrue(gene_list1[2] != gene_list2[2])
        self.assertTrue(gene_list1[1] != gene_list2[1])



if __name__ == "__main__":
    unittest.main()