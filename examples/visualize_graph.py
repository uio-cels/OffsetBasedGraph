

from offsetbasedgraph import Graph, Block
from offsetbasedgraph import VisualizeHtml



levels = {
    "chr1-1": 0,
    "chr1-merged": 1,
    "chr1-2": 0,
    "chr1-3": 0,
    "chr1-4": 0,
    "chr1_alt": 2,
    "chr1_alt3": 3,
    "chr1_alt2": 2
}

graph = Graph(
    {
        "chr1-1": Block(10),
        "chr1-2": Block(10),
        "chr1-3": Block(10),
        "chr1-4": Block(10),
        "chr1-merged": Block(2),
        "chr1_alt": Block(10),
        "chr1_alt2": Block(10),
        "chr1_alt3": Block(14),
    },
    {
        "chr1-1": ["chr1-2", "chr1_alt", "chr1_alt3"],
        "chr1_alt": ["chr1-merged"],
        "chr1-2": ["chr1-merged"],
        "chr1-merged": ["chr1_alt2", "chr1-3"],
        "chr1_alt2": ["chr1-4"],
        "chr1-3": ["chr1-4"],
        "chr1_alt3": ["chr1-3"]
    }
)

graph.start_block = "chr1-1"

from offsetbasedgraph.graphutils import Gene
from offsetbasedgraph import Interval

gene = Gene("test",
            Interval(3, 1, ["chr1-1", "chr1_alt", "chr1-merged"], graph),
            [Interval(4, 5, ["chr1-1"], graph)]
            )

print(gene)




v = VisualizeHtml(graph, 1, 10, 0, levels, "", 300, [gene])
print(v.get_wrapped_html())

# python3 visualize_graph.py > test.html | google-chrome test.html