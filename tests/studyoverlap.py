# Richard Darst, December 2011

import ctypes
import itertools
import numpy
import os
import random
import sys
import time
import types

import networkx

import pcd
import pcd.cmodels as cmodels
import pcd.multiresolution



def partition_nodes(nodes, Nsets, nodesPerSet):
    nodes = list(nodes)
    print len(nodes)
    random.shuffle(nodes)
    N = len(nodes)//Nsets
    N = min(N, nodesPerSet)
    sets = [ nodes[N*i:N*(i+1)] for i in range(Nsets) ]
    print [len(s) for s in sets]
    for s in sets:
        while len(s) < nodesPerSet:
            n = random.choice(nodes)
            if n in s:
                continue
            s.append(n)
    print [len(s) for s in sets]
    return sets

N = 1000
U = list(range(N))

level1 = None
size1 = 25
N1 = 40
P1 = .95
level1 = [ random.sample(U, size1) for _ in range(N1) ]
U = set().union(*level1)

level2 = None
size2 = 100
N2 = 10
P2 = .65
#level2 = [ random.sample(U, size2) for _ in range(N2) ]
level2 = partition_nodes(U, N2, size2)

PUniverse = 0

g = networkx.Graph()

def add_edges(g, points, prob, weight=-1):
    for i, x in enumerate(points):
        for j, y in enumerate(points[i+1:]):
            if random.uniform(0, 1) < prob:
                g.add_edge(x, y, weight=weight)

if level1:
    for points in level1:
        add_edges(g, points, P1)
if level2:
    for points in level2:
        add_edges(g, points, P2)
add_edges(g, list(g.nodes()), PUniverse)
print len(g)
#raw_input()

G = pcd.Graph.fromNetworkX(g)
G.make_sparse(default='auto')

if level1: level1 = [ [ G._nodeIndex[n] for n in ns ] for ns in level1 ]
if level2: level2 = [ [ G._nodeIndex[n] for n in ns ] for ns in level2 ]

if level1: print G.cmtyConnectivity(nodes=level1)
if level2: print G.cmtyConnectivity(nodes=level2)
#sys.exit(0)



# Make G0 be the "known" graph.
G0 = G.copy()
for c, contents in enumerate(level1):
#for c, contents in enumerate(level2):
    for n in contents:
        #n = G0._nodeIndex[i]
        G0.cmty[n] = c
assert numpy.all(G0.cmty != pcd.models.NO_CMTY)
G0.cmtyListInit()
#for c, contents in enumerate(level1):
for c, contents in enumerate(level2):
    for n in contents:
        #n = G0._nodeIndex[i]
        if not G0.cmtyContains(c, n):
            G0.cmtyListAddOverlap(c, n)



G2 = pcd.Graph(N=10)
print
G2.cmty[:] = 0
G2.cmty[0:2] = 1
G2.cmtyListInit()
G2.oneToOne = 0
G2.cmtyListAddOverlap(0, 0)
G2.cmtyListAddOverlap(0, 1)
G2.combineSubsets()
#sys.exit(0)



G.greedy(gamma=.1)
#print PO(G0, G)
G.ovgreedy(gamma=.1)
#print PO(G0, G)
#from fitz import interactnow
#sys.exit(0)

#sys.exit(0)
import pcd.F1
MR = pcd.MultiResolution(minimizerargs=dict(minimizer="ovfree",
                                            trials=5),
                         output='tests-output/studyoverlap/tmp7c.txt',
                         #G0=G0,
                         overlap=3,
                         #pairstyle='all',
                         )
MR.G0 = G0
MR.run([G]*3,
      #logGammaArgs=dict(low=.01, high=100, density=10),
      #logGammaArgs=dict(low=.01, high=1000, density=5),
      gammas=dict(low=.1, high=10, density=15),
      threads=1)
MR.write('tests-output/studyoverlap/tmp7c.txt')
MR.plot('tests-output/studyoverlap/tmp7c.png')

