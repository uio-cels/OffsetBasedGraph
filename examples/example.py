from gendatafetcher.ucscdb import DbWrapper
from offsetbasedgraph import OffsetBasedGraph, LinearInterval

"""
Example:
We get alternative loci from the UCSC databse, and create a graph for GRCh38
"""

# Initiate an empty graph with the name "grch38"
graph = OffsetBasedGraph("grch38")

# Create a connection to the ucsc database and get chrom sizes
db = DbWrapper()
chrom_sizes = db.get_chrom_lengths()

# Get all alt locis from the db
alt_loci_infos = db.get_alt_loci_infos(False)

# Add each non-alt chromosome to the graph
for chromosome, length in chrom_sizes.items():
    if "_" not in chromosome:
        graph.add_chromosome("hg38", chromosome, length)

# Now, add every alternative loci
for info in alt_loci_infos:
    main_interval = LinearInterval(
        "hg38", info["chrom"], info["chromStart"], info["chromEnd"])

    alt_interval = LinearInterval(
        "hg38", info["name"], 0, info["length"])

    graph.merge_lin_refs(main_interval, alt_interval)


# Print the blocks and edges
print(graph.blocks)
print(graph.block_edges)

