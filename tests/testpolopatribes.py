# Richard Darst, August 2011

import numpy

import pcd.models
import pcd.graphs

# weights of -1, 2, 6 seem to make it match the paper... for q.
graph = pcd.graphs.polopa_tribes(weightAllied=-1, weightHostile=2)
#graph = pcd.graphs.relabel_nodes(graph)
G = pcd.models.Graph.fromNetworkX(graph, defaultweight=1)

print G.interactions

assert numpy.all(G.interactions - G.interactions.T == 0)

G.minimize(gamma=.6)
#G.viz()
print G._test()
print G.interactions
#from fitz import interactnow
#exit(2)

MR = pcd.models.MultiResolution(low=.01, high=10)
MR.do(Gs=[G]*12, trials=10)
MR.write("tmp-polopatribes.txt")
#MR.viz()
