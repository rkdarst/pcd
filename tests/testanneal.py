# Richard Darst, October 2011

import os

import pcd
import pcd.graphs

dirname = 'tests-output/anneal/'
if not os.path.exists(dirname):
    os.makedirs(dirname)

#graph = pcd.graphs.dolphins(weightFriend=-1)
##graph = pcd.graphs.relabel_nodes(graph)
#G = pcd.Graph.fromNetworkX(graph, defaultweight=1)


#g = pcd.graphs.karate_club()
#G = pcd.Graph.fromNetworkX(g, defaultweight=1)


G = pcd.graphs.bss2d_n240_T050()


print "N nodes", G.N
gamma = .001
G.verbosity = 1.5
G.cmtyCreate()
G.trials(gamma, trials=10, minimizer=G.anneal, attempts=10, Escale=10000)
print G.energy(gamma)
G.trials(gamma, trials=10, minimizer=G.minimize)
print G.energy(gamma)
#G.minimize(.1)
#G.anneal(.1)


MR = pcd.MultiResolution(
    minimizer='trials',
    minimizerargs=dict(minimizer='anneal', Escale=10000, attempts=10,
                       trials=5,
                   betafactor=1.1),
    savefigargs=dict(fname=dirname+'gamma%(gamma)09.4f.png'))

MR.do([G]*2, logGammaArgs=dict(low=.01, high=10, density=2))
