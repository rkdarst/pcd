# Richard Darst, July 2011

import copy

from math import exp, log

import models
import util

approxequal = lambda x,y: abs(x-y) < 5e-11*(x+y)
gamma = 1.0

# Initialized to one community per list
G = models.random_graph(size=5)

# Do we have the right number of initial communites - one per particle?
assert G.q == 25
assert G.q == len(set(G.cmty))
G.minimize(gamma)
assert G.q == len(set(G.cmty))
G.cmtyCreate(randomize=False)

# Test energies
assert G.energy(gamma) == sum(G.energy_cmty(gamma, c) for c in range(G.Ncmty))
G.minimize(gamma)
assert G.energy(gamma) == sum(G.energy_cmty(gamma, c) for c in range(G.Ncmty))
G.cmtyCreate(randomize=False)

# Test entropies
# Evenly distribution
s = G.entropy
calculatedS = -(25 * (1/25.) * log(1/25., 2))
#print s, calculatedS, s-calculatedS
assert round(s-calculatedS, 10)==0, "Etropy calculation: %s"%((s, calculatedS))
# Spiked distribution, all in one community, entropy should be zero.

G.cmty[:] = 0
G.cmtyListInit()
assert G.entropy == 0

# Reset to normal
G.cmtyCreate(randomize=False)

# Test mutual information with itself
mi = util.mutual_information(G, G)
assert mi == s
G.minimize(gamma=gamma)
assert approxequal(util.mutual_information(G, G), G.entropy), (util.mutual_information(G, G), G.entropy)


# Test a case where mutual information should be zero.
G.cmtyCreate(randomize=False)
G2 = G.copy()
print id(G), id(G2), id(G.cmty), id(G2.cmty)
for i in range(5):
#    G.cmty[5*(i):5*(i+1)] = i
    G2.cmty[i::5] = i
    pass
print G.cmty
print G2.cmty
G.cmtyListInit()
G2.cmtyListInit()
print G.entropy, G2.entropy
print util.mutual_information(G, G2)

#M = MultiResolution(low=.01, high=100, G=G)
#M.do()
#M.print_()
