

import networkx

# Test the statistics function
import pcd.util
import pcd.graphs
def test_graph_stats(g):
    pcd.util.graph_stats(g, level=2, _test=True)

# Complete graph
test_graph_stats(networkx.complete_graph(10))

# Random graphs:
test_graph_stats(networkx.gnp_random_graph(10, .9))
test_graph_stats(networkx.gnp_random_graph(50, .5))
test_graph_stats(networkx.gnp_random_graph(100, .1))

# Directed
test_graph_stats(networkx.bipartite_random_graph(100, 100, .25, directed=True))

# Real graph:
test_graph_stats(pcd.graphs.karate_club())

