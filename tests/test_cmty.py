# Richard Darst, October 2013

import os
from os.path import join, dirname

import networkx

import pcd.cmty

datadir = join(dirname(pcd.__file__), 'data')
fname = join(datadir, 'test_cmty.txt')

# Test iterator
cmtys = pcd.cmty.CommunityListIterator(fname)

#print cmtys.__dict__

assert cmtys.label == 'test-communities'
assert len(tuple(cmtys.iteritems())) == 2 # Do a full iteration
assert len(cmtys) == 2

pcd.cmty._test_interface(cmtys)

#print cmtys.__dict__

cmtys_full = cmtys.to_dict()

assert cmtys_full.label == 'test-communities'

assert cmtys_full.cmtynames() == set(('one', 'two'))
#from fitz import interactnow

assert len(cmtys_full) == 2

pcd.cmty._test_interface(cmtys_full)
