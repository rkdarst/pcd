# Richard Darst, October 2013

import os
from os.path import join, dirname

import networkx
import pcd.cmty

datadir = join(dirname(pcd.__file__), 'data')
fname = join(datadir, 'test_cmty.txt')
fname2 = join(datadir, 'test_cmty2.txt')

# Test full object
cmtynodes = dict(
    a=set((0, 1, 2, 3)),
    b=set((3, 4, 5, 6)),
    c=set((7, 8, 9)),
    )
cmtys = pcd.cmty.Communities(cmtynodes)
assert cmtys.is_cover() == True
assert cmtys.is_non_overlapping() == False
assert cmtys.is_partition() == False
cmtys.load_networkx(networkx.complete_graph(10))

# Test iterator
cmtys = pcd.cmty.CommunityListIterator(fname)
#print cmtys.__dict__
assert cmtys.label == 'test-communities'
assert len(tuple(cmtys.iteritems())) == 2 # Do a full iteration
assert len(cmtys) == 2
pcd.cmty._test_interface(cmtys)

# Test conversion to full communities.
cmtys_full = cmtys.to_full()
assert cmtys_full.label == 'test-communities'
assert cmtys_full.cmtynames() == set(('one', 'two'))
#from fitz import interactnow
assert len(cmtys_full) == 2
pcd.cmty._test_interface(cmtys_full)

# Test label loading
cmtys = pcd.cmty.CommunityFile(fname2)
assert cmtys.label == 'test-communities2'


cmtys = pcd.cmty.Communities(cmtynodes)
def filter(c, ns):
    yield str(c)+'.1', ns
    yield str(c)+'.2', ns
cmtys_filter = pcd.cmty.CommunityFilter(cmtys, filter)
cmtys_filter_full = cmtys_filter.to_full()
cmtys_filter_full.q == cmtys_filter.q
pcd.cmty._test_interface(cmtys_filter)
assert cmtys_filter_full.nodes_spanned() == set(range(10))
assert cmtys_filter_full.cmtysizes_sum() == 22
assert cmtys_filter_full.overlap() == 2.2

cmtys = pcd.cmty.Communities(cmtynodes)
def filter(c, ns):
    if len(ns) == 4:
        yield c, ns
cmtys_filter = pcd.cmty.CommunityFilter(cmtys, filter)
assert len(cmtys_filter) == 2
