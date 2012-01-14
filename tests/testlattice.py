# Richard Darst, August 2011

import numpy

import networkx
from networkx.generators.classic import grid_graph

import util
import pcd


g = grid_graph([20,20], periodic=True)
for a,b,d in g.edges(data=True):
    d['weight'] = -1

G = pcd.Graph.fromNetworkX(g, defaultweight=1)
print (numpy.sum(G.imatrix) - numpy.sum(G.imatrix.diagonal()))/(400**2-400.)
#exit()

assert numpy.all(G.imatrix - G.imatrix.T == 0)

MR = pcd.MultiResolution()
MR.run([G]*5, gammas=dict(low=.1, high=100, density=10))
MR.write('tmp-lattice.txt')
#MR.viz()

