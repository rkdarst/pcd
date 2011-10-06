# Richard Darst, September 2011

from numpy import random
import os
#random.seed(0)

import pcd
import pcd.graphs

outputdir = os.path.join("tests-output/overlap")
if not os.access(outputdir, os.F_OK):
    os.makedirs(outputdir)


def overlapDegree(G):
    return sum(len(G.cmtyContents(c)) for c in G.cmtys()) - G.N


#
# Part 1 of tests: non-periodic.
#
gamma = 1
graph = pcd.graphs.dolphins(weightFriend=-1)
#graph = pcd.graphs.relabel_nodes(graph)

import networkx.drawing.layout as layout
layout = layout.spring_layout(graph)

G = pcd.Graph.fromNetworkX(graph, defaultweight=1, layout=layout)
G.minimize_trials(gamma=gamma, trials=10)

print G.q
print G.energy(gamma=gamma)
print sorted(G.n_counts().iteritems())

G.overlapMinimize(gamma=gamma)
print G.q
print G.energy(gamma=gamma)
print sorted(G.n_counts().iteritems())


G = pcd.graphs.fractalsquare(L=16)
G.minimize(gamma=.1)
G.overlapMinimize(gamma=.1)

G.savefig(os.path.join(outputdir, 'blah.png'),
          G._layout, hulls=True, base_radius=.1)

#G.viz()
#exit(2)

#G.minimize(gamma=10)
#G.viz()
#exit()

#MR = pcd.MultiResolution(low=.007, high=11)
#MR.do(Gs=[G]*12, trials=10)
#MR.write("tmp-dolphins.txt")
#MR.viz()





#
# This tests periodic system for overlaps (mainly to test convex hulls
# and making them not overlap)
#

import support.gromacs
G = support.gromacs.get_test_G()
G.setOverlap(True)

G.minimize(gamma=.1)
#G.overlapMinimize(gamma=.1)

G.savefig(os.path.join(outputdir, 'amorphous.png'),
          G._layout, hulls=True, base_radius=.1, periodic=G._periodic)
print G.q
print G.n_counts()
print "ov degree:", overlapDegree(G)


G.minimize_trials(gamma=.1, trials=10)
G.savefig(os.path.join(outputdir, 'amorphous_10trials.png'),
          G._layout, hulls=True, base_radius=.1, periodic=G._periodic)
print G.q
print G.n_counts()
print "ov degree:", overlapDegree(G)
