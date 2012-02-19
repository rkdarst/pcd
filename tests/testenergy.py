# Richard Darst, February 2012

import pcd

import pcd.graphs
from pcd.util import floateq
from pcd import cmodels

def compareEnergy(G1, G2):
    for gamma in (.01, 1, 100):
        assert floateq(G1.energy(gamma),
                       G2.energy_sparse(gamma))
        for c in G1.cmtys():
            #assert floateq(G1.energy_cmty(gamma, c),
            #               G2.energy_cmty_sparse(gamma, c))
            for n in G1.cmtyContents(c):
                assert floateq(G1.energy_cmty_n(gamma, c, n),
                               #G2.energy_cmty_n_sparse(gamma, c, n))
                               cmodels.energy_cmty_n_sparse(G2, gamma, c, n))

        for c2 in G1.cmtys():
            assert floateq(G1.energy_cmty_cmty(gamma,c,c2),
                           cmodels.energy_cmty_cmty_sparse(G2,gamma,c,c2))

G = pcd.graphs.polopa_tribes_G()

G.trials(gamma=.1, trials=10)
G2 = G.copy()
#G2.make_sparse('auto')
G2.make_sparse(default=1, cutoff_op=lambda x,y:x!=y)

print G2.simatrixDefault, G2.srmatrixDefault
print G2.simatrixN
print G2.simatrixId
#G2.srmatrixDefault = 5
print G.srmatrixDefault
print G.hasSparse, G2.hasSparse

compareEnergy(G, G2)
print G.energy(.5), G2.energy_sparse(.5)

