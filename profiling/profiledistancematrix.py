import cProfile
import pstats
import offsetbasedgraph as obg
from offsetbasedgraph.distancematrix import DistanceIndex

data_folder = "../../graph_peak_caller/tests/"


def run(graph, N):
    d_index = DistanceIndex(graph, N)
    d_index.create()
    d_index.to_file("distance_index")

if __name__ == "__main__":
    obgraph = obg.GraphWithReversals.from_file(data_folder + "obgraph")
    n = 500
    cProfile.run("run(obgraph, %s)" % n, "distance_matrix%s.profile" % n)
    p = pstats.Stats("distance_matrix%s.profile" % n)
    p.sort_stats("tottime").print_stats()
