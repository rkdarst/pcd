# Richard Darst, November 2011

import os
import random

import networkx
import numpy

import networkx

import pcd

if not os.access('tests-output/overlapmin/', os.F_OK):
    os.mkdir('tests-output/overlapmin/')

size = 25
number = 10

nodes = set(range(size*number))

cmtys = [ ]
for i in range(number):
    cmtys.append(set(range(size*i, size*i+size)))
for cnodes in cmtys:
    extraNodes = 0
    while len(cnodes) < size+10:
        n = random.choice(tuple(nodes))
        cnodes.add(n)
for c in cmtys:
    print len(c)
    print " ".join(str(x) for x in sorted(c))

g = networkx.Graph()

g.add_nodes_from(nodes)


for c in cmtys:
    nodes = list(c)
    for i, n1 in enumerate(nodes):
        for n2 in nodes[i+1:]:
            # iterating over pairs of nodes n1, n2
            if random.uniform(0,1) < .9:
                g.add_edge(n1, n2, weight=-1)

nodes = g.nodes()
for i,n1 in enumerate(nodes):
    for n2 in nodes[i+1:]:
        if random.uniform(0,1) < 0.15:
            g.add_edge(n1, n2, weight=-1)

G = pcd.Graph.fromNetworkX(g, defaultweight=1)

trials = 10 ; density = 10
if 'fast' in globals():
    trials = 2 ; density = 3

G.trials(gamma=1, trials=trials)
G.trials(gamma=1, trials=trials, minimizer='ovGreedy')


#nx.drawing.nx_pydot.pydot_layout(g)
#import networkx.drawing.layout as layout
import networkx.drawing.nx_pydot #layout as layout
#G.coords = layout.spring_layout(g)
g2 = g.copy()
for a,b,d in g2.edges_iter(data=True):
    if G.cmty[a] == G.cmty[b]:
        d['weight'] = 5
    else:
        d['weight'] = 1
    #d['weight'] = -d['weight']
G.coords = networkx.drawing.nx_pydot.pydot_layout(g2)
print G.coords

#G.savefig('tests-output/overlapmin/gamma1.png', hulls=True)
#G.viz()
#exit()

trials=10
density=10
if 'fast' in globals():
    trials = 2 ; density = 3


MR = pcd.MultiResolution(
    #savefigargs=dict(fname='tests-output/overlapmin/mr.png'),
    minimizerargs=dict(minimizer='ovGreedy', trials=trials),
    plotargs=dict(fname='tests-output/overlapmin/mr.png'),
    overlap=5)
MR.run([G]*5, gammas=dict(low=.001,high=100, density=density))

