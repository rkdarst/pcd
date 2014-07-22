# Richard Darst, October 2013

import os
from os.path import join, dirname
import unittest

import networkx
import pcd
import cmty

cmtynodes = dict(
    a=set((0, 1, 2, 3)),
    b=set((3, 4, 5, 6)),
    c=set((7, 8, 9)),
    )

def test_full():
    # Test full object
    cmtys = cmty.Communities(cmtynodes)
    assert cmtys.is_cover() == True
    assert cmtys.is_non_overlapping() == False
    assert cmtys.is_partition() == False
    cmtys.load_networkx(networkx.complete_graph(10))


    # Test object with missing some nodes (non-consecutive).  This is
    # mainly a test of the cmtyintmap and nodeintmap.
    cn2 = {
        0:set((0, 1, 2, 3)),
        1:set((3, 4, 5, )),
        5:set((7, 8, )),    #missing 6 and 9
        }
    cm2 = cmty.Communities(cn2)
    cmty._test_interface(cm2)


def test_iter():
    datadir = join(dirname(pcd.__file__), 'data')
    fname = join(datadir, 'test_cmty.txt')
    fname2 = join(datadir, 'test_cmty2.txt')

    # Test iterator
    cmtys = cmty.CommunityListIterator(fname)
    #print cmtys.__dict__
    assert cmtys.label == 'test-communities'
    assert len(tuple(cmtys.iteritems())) == 2 # Do a full iteration
    assert len(cmtys) == 2
    cmty._test_interface(cmtys)


    # Test conversion to full communities.
    cmtys_full = cmtys.to_full()
    assert cmtys_full.label == 'test-communities'
    assert cmtys_full.cmtynames() == set(('one', 'two'))
    #from fitz import interactnow
    assert len(cmtys_full) == 2
    cmty._test_interface(cmtys_full)


    # Test label loading from the .names file
    cmtys = cmty.CommunityFile(fname2)
    assert cmtys.label == 'test-communities2'


def test_filter():
    # Test some filters
    cmtys = cmty.Communities(cmtynodes)
    def filter(c, ns):
        yield str(c)+'.1', ns
        yield str(c)+'.2', ns
    cmtys_filter = cmty.CommunityFilter(cmtys, filter)
    cmtys_filter_full = cmtys_filter.to_full()
    cmtys_filter_full.q == cmtys_filter.q
    cmty._test_interface(cmtys_filter)
    assert cmtys_filter_full.nodes_spanned() == set(range(10))
    assert cmtys_filter_full.cmtysizes_sum() == 22
    assert cmtys_filter_full.overlap() == 2.2

    cmtys = cmty.Communities(cmtynodes)
    def filter(c, ns):
        if len(ns) == 4:
            yield c, ns
    cmtys_filter = cmty.CommunityFilter(cmtys, filter)
    assert len(cmtys_filter) == 2

def test_union():
    # Test unions
    cn1 = dict(A=set((0, 1, 3)), B=set((1, 2, 3, 4)))
    cn2 = dict(A=set((0, 1, 3)), B=set((5, 6, 7, 8, 9)))
    c1 = cmty.Communities(cn1)
    c2 = cmty.Communities(cn2)
    cU = cmty.CommunityUnion((c1, c2))
    assert len(cU) == 3
    assert cU.nodes == set(range(10))
    assert set(frozenset(x) for x in cU.itervalues()) \
           == set((frozenset((0,1,3)),frozenset((1,2,3,4)),frozenset((5,6,7,8,9))))
    #print list(cU.iteritems())
    assert len(set(cU.iterkeys())) == 3
    assert len(list(cU.itervalues())) == 3
    cmty._test_interface(cU)

    # Same test, but without dup_ok.
    cU = cmty.CommunityUnion((cn1, cn2), dup_ok = True)
    assert len(cU) == 4
    assert len(set(cU.iterkeys())) == 4
    assert cU.nodes == set(range(10))
    #print list(cU.iteritems())
    cmty._test_interface(cU)

def test_cmty_graph():
    # Test cmty_graph:
    g = networkx.complete_graph(7)
    g.remove_edge(3, 5)
    g.remove_edge(1, 2)
    cmtys = cmty.Communities({0:set((0,1,2)), 1:set((3,4)), 'a':set((5,6))})
    cmty_graph = cmtys.cmty_graph(g)
    assert cmty_graph.node[0]['size'] == 3           # number of nodes
    assert cmty_graph.node[0]['weight'] == 2  # edges within cmty
    assert cmty_graph.node['a']['size'] == 2
    assert cmty_graph.node['a']['weight'] == 1
    assert cmty_graph[0][1]['weight'] == 6    # edges between
    assert cmty_graph[1]['a']['weight'] == 3
    assert cmty_graph['a'][0]['weight'] == 6

    # Test cmty_graph on directed graph:
    g = networkx.complete_graph(7, create_using=networkx.DiGraph())
    g.remove_edge(3, 5)
    g.remove_edge(1, 2)
    cmtys = cmty.Communities({0:set((0,1,2)), 1:set((3,4)), 'a':set((5,6))})
    cmty_graph = cmtys.cmty_graph(g)
    assert cmty_graph.node[0]['size'] == 3       # number of nodes
    assert cmty_graph.node[0]['weight'] == 5     # edges within cmty
    assert cmty_graph.node['a']['size'] == 2
    assert cmty_graph.node['a']['weight'] == 2
    assert cmty_graph[0][1]['weight'] == 6     # edges between
    assert cmty_graph[1][0]['weight'] == 6
    assert cmty_graph[1]['a']['weight'] == 3
    assert cmty_graph['a'][1]['weight'] == 4
    assert cmty_graph['a'][0]['weight'] == 6
    assert cmty_graph[0]['a']['weight'] == 6



g = networkx.Graph([
    (0,1),
    (1,2),
    (2,3),
    (3,4),
    (4,5),
    (5,6),
    ])
cmtys = cmty.Communities(
    {0: set((0,1,2)),
     1: set((2,3,4)),
     2: set((5,6)),
     })


g2 = networkx.Graph([
    (0,1), (0,2), (1,2),
    (3,4), (4,5), (5,6), (6,3),
    (2,3),
    ])
cmtys2 = cmty.Communities(
    {'a': set((0,1,2)),
     'b': set((3,4,5,6)),
     'c': set((2,3,4)),
     })


class TestSubgraph(unittest.TestCase):
    def test_single(self):
        "Subgraph of a single community"
        self.assertEqual(set(cmtys.subgraph(g, 1).nodes()),
                         set((2,3,4)))
    def test_shells_cmty(self):
        "Subgraph of community and neighboring communities"
        sg = cmtys.subgraph(g, 0, shells_cmty=1)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3,4)))
        sg = cmtys.subgraph(g, 1, shells_cmty=1)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3,4,5,6)))
        # multiple neighbors
        sg = cmtys.subgraph(g, 0, shells_cmty=2)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3,4,5,6)))
    def test_shells_neighbors(self):
        "Subgraph of community and neighboring nodes"
        sg = cmtys.subgraph(g, 0, shells_node=1)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3)))
        sg = cmtys.subgraph(g, 1, shells_node=1)
        self.assertEqual(set(sg.nodes()), set((1,2,3,4,5)))
        # more than one neighbor
        sg = cmtys.subgraph(g, 0, shells_node=2)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3,4)))

    def test_shells_overlaps(self):
        "Subgraph of community and neighboring nodes"
        sg = cmtys.subgraph(g, 0, shells_overlap=1)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3,4)))
        sg = cmtys.subgraph(g, 1, shells_overlap=1)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3,4)))
        # more than one neighbor
        sg = cmtys.subgraph(g, 0, shells_overlap=2)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3,4)))
    def test_shells_mixed(self):
        "Subgraph of community and neighboring nodes"
        # cmtys is applied first, the nodes
        sg = cmtys.subgraph(g, 0, shells_cmty=1, shells_node=1)
        self.assertEqual(set(sg.nodes()), set((0,1,2,3,4,5)))
        sg = cmtys.subgraph(g, 2, shells_cmty=1, shells_node=1)
        self.assertEqual(set(sg.nodes()), set((1,2,3,4,5,6)))

class TestCmtygraph(unittest.TestCase):
    def test_basic(self):
        cg = cmtys.cmty_graph(g)
        self.assertEqual(len(cg), 3)
        self.assertEqual(set(cg.nodes()), set((0,1,2)))
        self.assertEqual(set(frozenset(x) for x in cg.edges()),
                         set(frozenset(x) for x in ((0,1),(1,2))))
    def test_weights(self):
        cg = cmtys.cmty_graph(g)
        self.assertEqual(cg.node[0]['weight'], 2)
        self.assertEqual(cg.node[1]['weight'], 2)
        self.assertEqual(cg.node[2]['weight'], 1)
        self.assertEqual(cg.edge[0][1]['weight'], 2)
        self.assertEqual(cg.edge[1][2]['weight'], 1)

    def test_basic2(self):
        cg = cmtys2.cmty_graph(g2)
        self.assertEqual(len(cg), 3)
        self.assertEqual(set(cg.nodes()), set('abc'))

        self.assertEqual(set(frozenset(x) for x in cg.edges()),
                         set(frozenset(x) for x in
                             (('a','b'),('b','c'),('a','c'))))
    def test_weights2(self):
        cg = cmtys2.cmty_graph(g2)
        self.assertEqual(cg.node['a']['weight'], 3)
        self.assertEqual(cg.node['b']['weight'], 4)
        self.assertEqual(cg.node['c']['weight'], 2)
        self.assertEqual(cg.edge['a']['b']['weight'], 1)
        # The following test is 5 becasue overlaps count twice.
        self.assertEqual(cg.edge['b']['c']['weight'], 5)
        self.assertEqual(cg.edge['a']['c']['weight'], 3)


# Test measures

gb = networkx.union(
    networkx.complete_graph(4),
    networkx.relabel_nodes(networkx.complete_graph(3), {0:4,1:5,2:6}),
    )
gb.add_edge(0, 4)
cb = cmty.Communities({0:set((0,1,2,3)), 1:set((4,5,6))})


from nose.tools import *

class TestMeasures(unittest.TestCase):
    def test_embeddedness(self):
        assert_equal(cb.cmty_embeddedness(gb),
                     {0:(3*4/float(3*4+1)),  1:2*3/float(2*3+1)})

        g = networkx.complete_graph(5)
        c = cmty.Communities({0:set(range(5))})
        assert_equal(c.cmty_embeddedness(g)[0],
                     1.0)

        g = networkx.path_graph(5)
        c = cmty.Communities({0:set(range(5))})
        assert_equal(c.cmty_embeddedness(g)[0],
                     1.0)

        # Singleton comm
        g = networkx.complete_graph(1)
        c = cmty.Communities({0:set(range(1))})
        assert_equal(c.cmty_embeddedness(g)[0],
                     0.0)

        # Empty comm
        g = networkx.complete_graph(0)
        c = cmty.Communities({0:set(range(0))})
        assert_equal(c.cmty_embeddedness(g)[0],
                     0.0)

        # Total embeddednesses
        ce = (4*3/(4*3+1.) * 4 + 3*2/(3*2+1.) * 3) / (4+3)
        assert_equal(cb.tot_embeddedness(gb), ce)

    def test_sld(self):
        assert_equal(cb.cmty_scaledlinkdensities(gb),
                     {0:4,  1:3})

        g = networkx.complete_graph(5)
        c = cmty.Communities({0:set(range(5))})
        assert_equal(c.cmty_scaledlinkdensities(g)[0],   5)

        g = networkx.path_graph(5)
        c = cmty.Communities({0:set(range(5))})
        assert_equal(c.cmty_scaledlinkdensities(g)[0],   2)

        # Singleton comm
        g = networkx.complete_graph(1)
        c = cmty.Communities({0:set(range(1))})
        assert_equal(c.cmty_scaledlinkdensities(g)[0],   0.0)

        # Empty comm
        g = networkx.complete_graph(0)
        c = cmty.Communities({0:set(range(0))})
        assert_equal(c.cmty_scaledlinkdensities(g)[0],   0.0)

        # Total SLD
        ce = (4 * 4 + 3 * 3) / float(4+3)
        assert_equal(cb.tot_scaledlinkdensity(gb), ce)
