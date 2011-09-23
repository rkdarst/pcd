# Richard Darst, August 2011

from numpy import random
random.seed(0)

import pcd.graphs

graph = pcd.graphs.dolphins(weightFriend=-1)
#graph = pcd.graphs.relabel_nodes(graph)
G = pcd.Graph.fromNetworkX(graph, defaultweight=1)
G.minimize(.1)
#G.viz()
#exit(2)

print G.imatrix

#G.minimize(gamma=10)
#G.viz()
#exit()

MR = pcd.MultiResolution()
MR.do(Gs=[G]*12, logGammaArgs=dict(low=.007, high=11), trials=10)
MR.write("tmp-dolphins.txt")
#MR.viz()
