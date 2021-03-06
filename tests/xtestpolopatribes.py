# Richard Darst, August 2011

import numpy
import random
import networkx.algorithms
import os

import pcd
import pcd.graphs
import pcd.util

fast = globals().get('fast', False)
if not os.access('tests-output/tribes/', os.F_OK):
    os.mkdir('tests-output/tribes/')

# weights of -1, 2, 6 seem to make it match the paper... for q.
g = pcd.graphs.polopa_tribes(weightAllied=-1, weightHostile=2)
g_p2 = pcd.graphs.polopa_tribes_peter(fname="highland.txt")
g_p1 = pcd.graphs.polopa_tribes_peter()

assert networkx.algorithms.is_isomorphic(g, g_p2)
assert networkx.algorithms.is_isomorphic(g, g_p1, weighted=True)

g_gml = g.copy()
for n, dat in g_gml.nodes(data=True):
    dat['name'] = n
for n1, n2, dat in g_gml.edges(data=True):
    w = dat['weight']
    if w == -1:
        dat['weight'] = 100
    else:
        dat['weight'] = -100
        g_gml.remove_edge(n1, n2)
import networkx
networkx.write_gml(g_gml, 'tribes.gml')

G = pcd.Graph.fromNetworkX(g, defaultweight=1, diagonalweight=0)
G2 = pcd.Graph.fromNetworkX(g_p1, defaultweight=1, diagonalweight=0)
print G.imatrix
assert numpy.all(G.imatrix - G.imatrix.T == 0)

from numpy import array
layout = {'GAVEV': array([ 0.53410183,  0.06217559]), 'KOHIK': array([ 0.       ,  0.1474209]), 'UHETO': array([ 0.14830981,  0.34515878]), 'KOTUN': array([ 0.37326864,  0.        ]), 'NAGAM': array([ 0.26417734,  0.41649365]), 'ASARO': array([ 0.56901665,  0.49511319]), 'MASIL': array([ 0.55606811,  0.41225591]), 'NOTOH': array([ 0.07688195,  0.24116678]), 'UKUDZ': array([ 0.70980807,  0.41850366]), 'ALIKA': array([ 1.        ,  0.17964178]), 'GEHAM': array([ 0.43121325,  0.39188928]), 'GAHUK': array([ 0.56559475,  0.2937324 ]), 'GAMA': array([ 0.3887513 ,  0.17512021]), 'OVE': array([ 0.80835312,  0.17731299]), 'NAGAD': array([ 0.28961018,  0.10020374]), 'SEUVE': array([ 0.27922518,  0.62217477])}
G.coords = dict((n, layout[G._nodeLabel[n]]) for n in range(G.N))

G.trials(gamma=1, trials=10)
G.savefig('tests-output/tribes/gamma1.png', base_radius=.02)
G.savefig('tests-output/tribes/gamma1.svg', base_radius=.02)

while True:
    break
    gamma = 3
    from matplotlib import pyplot
    pyplot.ion()
    pyplot.clf()
    G.cmtyCreate()
    G.minimize(gamma=gamma)
    print gamma, G.q, G.energy(gamma=gamma)
    #if G._layout is None:
    #    import networkx.drawing.layout as layout
    #    G._layout = layout.spring_layout(g)
    #    print G._layout
    if fast and (1 or G.q != 5 and G.energy(gamma=gamma) <= -21):
        g = G.viz()
        raw_input(',')
#exit(2)

nreplicas, ntrials, ldensity = 12, 100, 100
nreplicas, ntrials, ldensity = 12, 5, 30
if fast:
    nreplicas, ntrials, ldensity = 3, 10, 10

#for a,b in [random.sample(range(G.N), 2) for _ in range(1000)]:
#    pcd.util.matrix_swap_basis(G.imatrix, a, b)
from pcd import LogInterval
import pcd.F1

if fast:
    savefigargs = None
else:
    savefigargs = dict(fname="tests-output/tribes/gamma%(gamma)05.3f.png")
MR = pcd.MultiResolution(overlap=5,
        savefigargs=savefigargs
                         )
G.make_sparse('auto')
G.verbosity = -1

MRR = pcd.MRRunner(MR)
#raw_input('before MRR >')
MRR.do(Gs=[G]*nreplicas,
       #gammas=LogInterval(low=.01, high=10, density=ldensity).values(),
       gammas=dict(low=.01, high=10, density=ldensity),
       )
#raw_input('after MRR >')
MR.write("tests-output/tribes/tribes.txt")
f = MR.plot("tests-output/tribes/tribes.png",
            ax1items=('VI', 'In', #'s_F1'
                      ))
#MR.viz()
f.axes[1].set_ylim(ymin=0, ymax=6)
f.savefig("tests-output/tribes/tribes.png")

#for a,b in [random.sample(range(G2.N), 2) for _ in range(1000)]:
#    pcd.util.matrix_swap_basis(G2.imatrix, a, b)
#MR = pcd.MultiResolution(low=.01, high=10)
#MR.do(Gs=[G2]*12, trials=50)
#MR.write("tests-outputs/tribes/tmp-polopatribes2.txt")
