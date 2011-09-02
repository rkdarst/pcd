# Richard Darst, August 2011

import numpy
import random
import networkx.algorithms

import pcd.models
import pcd.graphs
import pcd.util
pcd.util.logNumber = 100

# weights of -1, 2, 6 seem to make it match the paper... for q.
g = pcd.graphs.polopa_tribes(weightAllied=-1, weightHostile=2)
g_p2 = pcd.graphs.polopa_tribes_peter(fname="highland.txt")
g_p1 = pcd.graphs.polopa_tribes_peter()

assert networkx.algorithms.is_isomorphic(g, g_p2)
assert networkx.algorithms.is_isomorphic(g, g_p1, weighted=True)


G = pcd.models.Graph.fromNetworkX(g, defaultweight=1, diagonalweight=0)
G2 = pcd.models.Graph.fromNetworkX(g_p1, defaultweight=1, diagonalweight=0)
print G.interactions
assert numpy.all(G.interactions - G.interactions.T == 0)

while True:
    G.minimize(gamma=.75)
    #G2.minimize(gamma=.75)
    #G.viz()
    #print G._test()
    #print G.interactions
    #print G2.interactions
#from fitz import interactnow
    break
#exit(2)

#for a,b in [random.sample(range(G.N), 2) for _ in range(1000)]:
#    pcd.util.matrix_swap_basis(G.interactions, a, b)
MR = pcd.models.MultiResolution(low=.01, high=10)
MR.do(Gs=[G]*12, trials=10)
MR.write("tmp-polopatribes.txt")
#MR.viz()

#for a,b in [random.sample(range(G2.N), 2) for _ in range(1000)]:
#    pcd.util.matrix_swap_basis(G2.interactions, a, b)
#MR = pcd.models.MultiResolution(low=.01, high=10)
#MR.do(Gs=[G2]*12, trials=50)
#MR.write("tmp-polopatribes2.txt")
