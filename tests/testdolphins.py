# Richard Darst, August 2011

import os
from numpy import random
random.seed(0)

import pcd.graphs

if not os.access('tests-output/dolphins', os.F_OK):
    os.mkdir('tests-output/dolphins')

graph = pcd.graphs.dolphins(weightFriend=-1)
#graph = pcd.graphs.relabel_nodes(graph)
G = pcd.Graph.fromNetworkX(graph, defaultweight=1)
G.minimize(.1)
#G.viz()
#exit(2)

print G.imatrix

#G.minimize(gamma=10)
#G.viz()
#exit()
import pcd.F1
MR = pcd.MultiResolution(overlap=5)
MR.run(Gs=[G]*12, gammas=dict(low='auto', high=11))

MR.write("tests-output/dolphins/tmp-dolphins.txt")
MR.plot(fname='tests-output/dolphins/mr.png',
        ax1items=('VI', 'In', 'ov_N', 'entropy', 's_F1'))
#MR.viz()
