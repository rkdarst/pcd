# Richard Darst, August 2011

import pcd.models
import pcd.graphs

graph = pcd.graphs.dolphins(weightFriend=-1)
#graph = pcd.graphs.relabel_nodes(graph)
G = pcd.models.Graph.fromNetworkX(graph, defaultweight=1)
G.minimize(.1)
#G.viz()
#exit(2)

print G.interactions

#G.minimize(gamma=10)
#G.viz()
#exit()

MR = pcd.models.MultiResolution(low=.007, high=11)
MR.do(Gs=[G]*12, trials=10)
MR.write("tmp-dolphins.txt")
#MR.viz()
