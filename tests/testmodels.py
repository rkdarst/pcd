# Richard Darst, July 2011

import copy

from math import exp, log

import models
import util

# Initialized to one community per list
G = models.random_graph(size=5)
G.cmtyListCheck()

G2 = G.copy()
assert id(G) != id(G2)
assert id(G) != id(G2)


#
# Test related to cmty list data structures
#
G.cmtyListCheck()

# Test that cmtyListCheck does find errors:
G.cmty[:] = 0
try:                     G.cmtyListCheck()
except AssertionError:   pass
else:        raise AssertionError("cmtyListCheck did not detect an error "
                                  "when there should be one.")
G.cmtyListInit()
G.cmtyListCheck()

# Test various ways of manipulating the lists
G.minimize(1.0)
G.cmtyListCheck()

G.minimize(0)
G.cmtyListCheck()

G.remapCommunities()
G.cmtyListCheck()
