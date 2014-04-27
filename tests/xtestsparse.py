# Richard Darst, Nomember 2011

from numpy import random ; random.seed(50)
import sys

import pcd
import pcd.graphs

G = pcd.graphs.bss2d_n240_T050()
#G = pcd.graphs.bss2d_n5760_T040()
#G = pcd.graphs.bss2d_n23040_T040()

G2 = G.copy()
print "reducing"
G.make_sparse(cutoff=499.1,
              default=500)
#G._check_sparse(imatrixDefault=500)
#raw_input("waiting")
#G2.minimize(gamma=1)
#G.minimize(gamma=1)
#from fitz import interactnow

#G.cmtyCreate()
gamma=.0001
if len(sys.argv) > 1 and sys.argv[1] != 's':
    G2.trials(gamma=gamma, trials=10)
else:
    G.trials(gamma=gamma, trials=10)
    #G.trials(gamma=gamma, minimizer='overlapMinimize', trials=10,
    #         initial='current')

MR = pcd.MultiResolution()
MR.run([G]*5, gammas=dict(low=.001,high=100),
      threads=2)

