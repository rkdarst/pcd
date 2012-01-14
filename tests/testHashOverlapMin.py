# Richard Darst, Nomember 2011

from numpy import random ; random.seed(50)
import sys

import pcd
import pcd.graphs

#G = pcd.graphs.bss2d_n240_T050()
G = pcd.graphs.bss2d_n5760_T040()
#G = pcd.graphs.bss2d_n23040_T040()

G2 = G.copy()
print "reducing"
G.make_sparse(cutoff=499.5,
              default=500)
#G._check_sparse(imatrixDefault=500)
#raw_input("waiting")
#G2.minimize(gamma=1)
#G.minimize(gamma=1)
#from fitz import interactnow

#G.cmtyCreate()
gamma=1
if len(sys.argv) > 1 and sys.argv[1] != 's':
    G2.trials(gamma=gamma, trials=10)
else:
    #G.trials(minimizer='overlapMinimize', gamma=gamma, trials=10)
    G.trials(minimizer='minimize', gamma=gamma, trials=10)


#MR = pcd.MultiResolution()
#MR.do([G]*5, logGammaArgs=dict(low=.001,high=100))
