# Richard Darst, November 2011

from numpy import random
import os

import pcd
import pcd.graphs


outputdir = os.path.join("tests-output/overlapMR")
if not os.access(outputdir, os.F_OK):
    os.makedirs(outputdir)


G = pcd.graphs.bss2d_n240_T050()

G.setOverlap(True)

#
# Test MR and overlaps
#
MR = pcd.MultiResolution(
    plotargs=dict(fname=os.path.join(outputdir,'MR.png'),
                  ax1items=('NmiO', 'entropy', 'In'),
                  plotstyles={'NmiO':dict(scale=10),
                              'In':dict(scale=3)}
                  ))
MR.do([G]*8, logGammaArgs=dict(low=.01, high=100))
