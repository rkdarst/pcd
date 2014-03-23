
import unittest

import networkx
import pcd.cmty

g = networkx.Graph([
    (0,1),
    (1,2),
    (2,3),
    (3,4),
    (4,5),
    (5,6),
    ])
cmtys = pcd.cmty.Communities(
    {0: set((0,1,2)),
     1: set((2,3,4)),
     2: set((5,6)),
     })


g2 = networkx.Graph([
    (0,1), (0,2), (1,2),
    (3,4), (4,5), (5,6), (6,3),
    (2,3),
    ])
cmtys2 = pcd.cmty.Communities(
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


if __name__ == "__main__":
    unittest.main()
