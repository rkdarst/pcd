# Richard Darst, September 2011

from numpy import random
#random.seed(0)

import pcd
import pcd.graphs

gamma = 1
graph = pcd.graphs.dolphins(weightFriend=-1)
#graph = pcd.graphs.relabel_nodes(graph)

import networkx.drawing.layout as layout
layout = layout.spring_layout(graph)

G = pcd.Graph.fromNetworkX(graph, defaultweight=1, layout=layout)
G.minimize_trials(gamma=gamma, trials=10)

print G.q
print G.energy(gamma=gamma)
print sorted(G.cmtyHistogram().iteritems())

G.overlapMinimize(gamma=gamma)
print G.q
print G.energy(gamma=gamma)
print sorted(G.cmtyHistogram().iteritems())


G = pcd.graphs.fractalsquare(L=16)
G.minimize(gamma=.1)
G.overlapMinimize(gamma=.1)

G.savefig('blah.png', G._layout, hulls=True, base_radius=.1)

#G.viz()
#exit(2)

#G.minimize(gamma=10)
#G.viz()
#exit()

#MR = pcd.MultiResolution(low=.007, high=11)
#MR.do(Gs=[G]*12, trials=10)
#MR.write("tmp-dolphins.txt")
#MR.viz()
