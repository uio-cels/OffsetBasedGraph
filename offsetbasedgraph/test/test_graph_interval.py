import unittest
from offsetbasedgraph import LinearInterval, GraphInterval
from dummygraph import graph, region_paths


class TestGraphInterval(unittest.TestCase):
    def test_main_from_linear(self):
        interval = LinearInterval("hg38", "main", 15, 35)
        graph_interval = GraphInterval.from_linear_interval(graph, interval)
        self.assertEqual(
            graph_interval.block_list,
            region_paths[1:4]
            )

    def test_alt_from_linear(self):
        interval = LinearInterval("hg38", "alt", 5, 25)
        graph_interval = GraphInterval.from_linear_interval(graph, interval)
        self.assertEqual(
            graph_interval.block_list,
            [region_paths[1], region_paths[-1], region_paths[3]]
            )


if __name__ == "__main__":
    unittest.main()
