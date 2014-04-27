# Richard Darst, July 2011

from math import exp, log

import pcd.graphs as graphs
import pcd.models as models
import pcd.util as util

approxequal = lambda x,y: abs(x-y) <= 5e-6*max(abs(x),abs(y))
gamma = 1.0

# Initialized to one community per list
G = graphs.random_graph(size=5)

# Do we have the right number of initial communites - one per particle?
assert G.q == 25
assert G.q == len(set(G.cmty))
G.minimize(gamma)
assert G.q == len(set(G.cmty))
G.cmtyCreate(randomize=False)

# Test energies
# energy = sum(energy_cmty)
assert G.energy(gamma) == sum(G.energy_cmty(gamma, c) for c in range(G.Ncmty))
G.minimize(gamma)
assert G.energy(gamma) == sum(G.energy_cmty(gamma, c) for c in range(G.Ncmty))
# test energy_cmty vs energy_cmty_n
for c in range(G.Ncmty):
    assert G.energy_cmty(gamma, c) == sum(G.energy_cmty_n(gamma, c, n)
                                          for n in G.cmtyContents(c))
# test energy_cmty_cmty
# FIXME: how can this be tested well?



# Test entropies
# Evenly distribution
G.cmtyCreate(randomize=False)
assert G.entropy_python == G.entropy_c
s = G.entropy
calculatedS = -(25 * (1/25.) * log(1/25., 2))
#print s, calculatedS, s-calculatedS
assert round(s-calculatedS, 10)==0, ("Entropy calculation: %s %s"%
                                     (s, calculatedS))
# Spiked distribution, all in one community, entropy should be zero.

G.cmty[:] = 0
G.cmtyListInit()
assert G.entropy == 0
assert G.entropy_python == G.entropy_c

# Reset to normal
G.cmtyCreate(randomize=False)

# Test mutual information with itself
mi = util.mutual_information(G, G)
print mi, s
assert approxequal(mi, s)
G.minimize(gamma=gamma)
assert approxequal(util.mutual_information(G, G), G.entropy), \
       (util.mutual_information(G, G), G.entropy)
assert approxequal(util.mutual_information_c(G, G),
                   util.mutual_information_python(G,G))


# Test a case where mutual information should be zero.
# FIXME this test case does not work
G.cmtyCreate(randomize=False)
G2 = G.copy()
print id(G), id(G2), id(G.cmty), id(G2.cmty)
for i in range(5):
    G.cmty[5*(i):5*(i+1)] = i
    G2.cmty[i::5] = i
    pass
print G.cmty
print G2.cmty
G.cmtyListInit()
G2.cmtyListInit()
print G.entropy, G2.entropy
print util.mutual_information_python(G, G2), \
      util.mutual_information_c(G, G2)
assert approxequal(util.mutual_information_c(G, G2),
                   util.mutual_information_python(G,G2))

# Mutual information test case 2 (is it nonzero?)
# FIXME this test case does not work
G2 = G.copy()
G .cmtyCreate(randomize=False)
G2.cmtyCreate(randomize=True)
G .cmtyListInit()
G2.cmtyListInit()
G .minimize(gamma=1)
G2.minimize(gamma=1)
print G .entropy
print G2.entropy
print util.mutual_information_python(G, G2), \
      util.mutual_information_c(G, G2)
assert approxequal(util.mutual_information_c(G, G2),
                   util.mutual_information_python(G,G2))

#M = MultiResolution(low=.01, high=100, G=G)
#M.do()
#M.print_()


# Random tests
for i in range(100):
    G.cmtyCreate(randomize=True)
    G.minimize(gamma)
    assert 0 == G.check()
    # energy
    assert G.energy(gamma) \
           == sum(G.energy_cmty(gamma, c) for c in range(G.Ncmty))
    s = G.entropy
    mi = util.mutual_information(G, G)
    print mi, s
    assert approxequal(mi, s)
    assert approxequal(util.mutual_information_c(G, G2),
                       util.mutual_information_python(G,G2))


# Overlap N tests
print "\n\n\nOverlap N tests"
#G.trials(gamma=gamma, trials=5, minimizer='overlapMinimize')
G2 = G.copy()
print G.N, G2.N, G.Ncmty, G2.Ncmty
for c0 in range(G.Ncmty):
    p = util.H2_python(G, G2, c0, 0)
    c = util.H2_c(G, G2, c0, 0)
    if p==float('inf') and c==float('inf'): continue
    assert round(p-c,10)==0, "%r %r"%(p, c)
    print p, c

nmio_python = util.mutual_information_overlap(G, G2)
# extremely unpythoic way of doing this
util.HX_Ynorm = util.HX_Ynorm_python
util.H2 = util.H2_python

nmio_c = util.mutual_information_overlap(G, G2)

print nmio_python, nmio_c


G.make_sparse('auto')
assert approxequal(G.energy(.1), G.energy_sparse(.1))
assert approxequal(G.energy(10), G.energy_sparse(10))
assert approxequal(G.energy(.5), G.energy_sparse(.5))
