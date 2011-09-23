# Richard Darst, August 2011

import pcd.graphs

graph = pcd.graphs.nussinov_256node(weight=-1)
#graph = pcd.graphs.relabel_nodes(graph)
G = pcd.Graph.fromNetworkX(graph, defaultweight=1)
G.minimize(.1)
#G.viz()
#exit(2)

print G.imatrix

#G.minimize(gamma=10)
#G.viz()
#exit()

from pcd import LogInterval
MR = pcd.MultiResolution()
MR.do(Gs=[G]*12, gammas=LogInterval(low=.06, high=12), trials=10)
MR.write("tmp-256node.txt")
#MR.viz()
