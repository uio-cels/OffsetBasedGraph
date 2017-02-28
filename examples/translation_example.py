from offsetbasedgraph import Block, Graph, Interval, Translation
graph = Graph(
    {
        "chr1": Block(20),
        "chr1_alt1": Block(10),
    },
    {}
)

gene = Interval(0, 8, ["chr1"], graph)
print(gene)

merge_translation = Translation(
    {
        "chr1_alt1": [Interval(0, 5, ["merged_part", "chr1_alt1-unmerged"])],
        "chr1": [Interval(0, 15, ["merged_part", "chr1-unmerged"])]
    },
    {
        "merged_part": [Interval(0, 5, ["chr1"]), Interval(0, 5, ["chr1_alt1"], graph)],
        "chr1-unmerged": [Interval(5, 20, ["chr1"], graph)],
        "chr1_alt1-unmerged": [Interval(5, 10, ["chr1_alt1"], graph)]
    },
    graph=graph
)

translated_interval = merge_translation.translate(gene)
print(translated_interval)