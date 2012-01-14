# Richard Darst, November 2011

from numpy import random
import os

import pcd
import pcd.graphs


outputdir = os.path.join("tests-output/overlapMR")
if not os.access(outputdir, os.F_OK):
    os.makedirs(outputdir)


G = pcd.graphs.bss2d_n240_T050()

#G.setOverlap(True)

#
# Test MR and overlaps
#
density=10
nreplicas=8
if 'fast' in globals():
    density=3
    nreplicas=3
MR = pcd.MultiResolution(
    overlap=5,
    plotargs=dict(
                  fname=os.path.join(outputdir,'MR.png'),
                  ax1items=('ov_N', 'entropy', 'In'),
                  plotstyles={'ov_N':dict(scale=10),
                              'In':dict(scale=3)}
                  ))
MR.run([G]*nreplicas, gammas=dict(low=.01, high=100, density=density))
