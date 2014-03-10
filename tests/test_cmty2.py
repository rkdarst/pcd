
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


if __name__ == "__main__":
    unittest.main()
