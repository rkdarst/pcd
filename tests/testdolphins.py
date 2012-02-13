# Richard Darst, August 2011

import networkx

import os
from numpy import random
random.seed(0)

import pcd.graphs

if not os.access('tests-output/dolphins', os.F_OK):
    os.mkdir('tests-output/dolphins')
dirname = 'tests-output/dolphins/'

coords = networkx.read_dot('/home/richard/dolphins-coords-messy.gv')

graph = pcd.graphs.dolphins(weightFriend=-1)
#graph = pcd.graphs.relabel_nodes(graph)
G = pcd.Graph.fromNetworkX(graph, defaultweight=1)
G.minimize(.1)
#G.viz()
#exit(2)
for node in graph.nodes():
    graph.node[node]['pos'] = coords.node[str(node)]['pos']

print G.imatrix

#G.minimize(gamma=10)
#G.viz()
for e in graph.edges():
    del graph.edge[e[0]][e[1]]['weight']

for gamma in ('1e-2', '4e-2', '1e-1', '2e-1', '7e-1'):
    G.cmtyCreate()
    G.trials(gamma=float(gamma), trials=100)
    G.remap()
    G.colors_to_networkx(graph)
    networkx.write_gml(graph.copy(), dirname+'dolphins_gamma%s.gml'%gamma)
#raise
#exit()
# gml2gv research/pcd/tests-output/dolphins/dolphins_gamma2e-1.gml | neato -Tpng -Gsplines=true -Goverlap=false -Epenwidth=2 -Npenwidth=2 -Nfontsize=15 | display

import pcd.F1
MR = pcd.MultiResolution(overlap=5)
MR.run(Gs=[G]*12, gammas=dict(low='auto', high=11, density=20))

MR.write("tests-output/dolphins/tmp-dolphins.txt")
f = MR.plot(fname='tests-output/dolphins/mr.png',
            ax1items=('VI', 'In', #'ov_N', #'entropy',
                      #'s_F1'
                      ),
            plotstyles={'H':dict(scale=.2)})
f.axes[0].set_ylim(ymin=0,)
f.savefig('tests-output/dolphins/mr.png')

#MR.viz()
