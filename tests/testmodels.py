# Richard Darst, July 2011

import copy
from math import exp, log
import sys
import types

import graphs
import models
import util

# Initialized to one community per list
G = graphs.random_graph(size=5)
G.cmtyListCheck()

# Test copying
G2 = G.copy()
assert id(G) != id(G2)
assert id(G) != id(G2)

# Test that various attributes do NOT have the same id:
for name in dir(G):
    if name[:2] == '__' \
       or isinstance(getattr(G, name), types.MethodType)\
       or name in ('_layout', 'q', ):
        continue
    print name, id(getattr(G, name)), id(getattr(G2, name)), \
          type(getattr(G, name))
    v1 = getattr(G , name)
    v2 = getattr(G2, name)
    assert id(v1) != id(v2), name
# Continue testing copying.  Make sure that 
G .cmtyListCheck()
G2.cmtyListCheck()
G .cmtyCreate(randomize=True)
G2.cmtyCreate(randomize=True)
assert (G.cmty != G2.cmty).any()
assert (G.cmtyl != G2.cmtyl).any()
#from fitz import interactnow
#sys.exit(1)
G .cmtyListCheck()
G2.cmtyListCheck()
G .minimize(1.0)
G2.minimize(.01) # use different gammas, otherwise it is too likely to
                 # be identical results
G .cmtyListCheck()
G2.cmtyListCheck()
assert (G.cmty != G2.cmty).any()
assert (G.cmtyl != G2.cmtyl).any()


#
# Test related to cmty list data structures
#
G.cmtyListCheck()

# Test that cmtyListCheck does find errors:
G.cmty[:] = 0
print "***Ignore errors below, this is a test of errorchecking***"
try:                    G.cmtyListCheck()
except AssertionError:
    print "***Ignore errors above, this is a test of errocchecking***"
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
