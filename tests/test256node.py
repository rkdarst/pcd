# Richard Darst, August 2011

import os

import pcd.graphs

if not os.access("tests-output/256node", os.F_OK):
    os.mkdir("tests-output/256node")

G = pcd.graphs.nussinov_256node_G()
G.greedy(.1)
#G.viz()
#exit(2)

print G.imatrix

#G.minimize(gamma=10)
#G.viz()
#exit()

MR = pcd.MultiResolution(overlap=True)
runner = MR.runner()
runner.do(Gs=[G]*12, gammas=dict(low=.06, high=12))
MR.write("tests-output/256node/mr.txt")
#MR.viz()
