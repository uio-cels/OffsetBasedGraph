import unittest
from dummygraph import *

simple_graph = get_simple_graph()


class TestGraph(unittest.TestCase):

    def test_alignment(self):
        graph = get_simple_graph()
        alignment1 = (
            LinearInterval("hg38", "main", 10, 20),
            LinearInterval("hg38", "alt", 0, 10))
        alignment2 = (
            LinearInterval("hg38", "main", 30, 40),
            LinearInterval("hg38", "alt", 20, 30))

        graph.include_alignments(
            [alignment1, alignment2])

        # for k, v in graph.blocks.items():
        #    print(k, v)

        first_block = graph.blocks["main"]
        self.assertEqual(first_block.linear_references,
                         {"main": LinearInterval("hg38", "main", 0, 10)})
        second_block = graph.get_next_blocks(first_block)[0]

        self.assertEqual(second_block.linear_references,
                         {"main": LinearInterval("hg38", "main", 10, 20),
                          "alt": LinearInterval("hg38", "alt", 0, 10)})
        third_blocks = graph.get_next_blocks(second_block)
        third_block_m = [b for b in third_blocks if "main"
                         in b.linear_references][0]

        self.assertEqual(third_block_m.linear_references,
                         {"main": LinearInterval("hg38", "main", 20, 30)})

        third_block_a = [b for b in third_blocks if "alt"
                         in b.linear_references][0]

        self.assertEqual(third_block_a.linear_references,
                         {"alt": LinearInterval("hg38", "alt", 10, 20)})

        fourth_block = graph.get_next_blocks(third_block_a)[0]
        self.assertEqual(fourth_block.linear_references,
                         {"main": LinearInterval("hg38", "main", 30, 40),
                          "alt": LinearInterval("hg38", "alt", 20, 30)})

if __name__ == "__main__":
    unittest.main()
