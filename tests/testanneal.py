# Richard Darst, October 2011

from numpy import random
random.seed(0)

import pcd
import pcd.graphs

#graph = pcd.graphs.dolphins(weightFriend=-1)
##graph = pcd.graphs.relabel_nodes(graph)
#G = pcd.Graph.fromNetworkX(graph, defaultweight=1)


#g = pcd.graphs.karate_club()
#G = pcd.Graph.fromNetworkX(g, defaultweight=1)


G = pcd.graphs.bss2d_n240()


print "N nodes", G.N
gamma = .001
G.verbosity = 1.5
G.cmtyCreate()
G.trials(gamma, trials=10, minimize=G.anneal, attempts=100) # , Escale=100000
print G.energy(gamma)
G.trials(gamma, trials=10, minimize=G.minimize)
print G.energy(gamma)
#G.minimize(.1)
#G.anneal(.1)
