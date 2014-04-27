# Richard Darst, 2013 April

import networkx

import pcd
import pcd.cmty

g = networkx.Graph()
g.add_edge(0, 1)
g.add_edge(2, 3)
g.add_edge(0, 2)
g.add_edge(1, 3)
g.add_edge(0, 3)
g.add_edge(1, 2)

g.add_edge(10, 11)
g.add_edge(12, 13)
g.add_edge(10, 12)
g.add_edge(11, 13)
g.add_edge(10, 13)
g.add_edge(11, 12)

#g.add_edge(0, 10)
#g.add_edge(3, 13)

G = pcd.Graph.fromNetworkX(g)

C = pcd.cmty.Communities(cmtynodes={0:set((0,1,2,3)), 1:set((10,11,12,13))}, nodes=set(g.nodes()))
C.load_pcd(G)

print G.modularity()




import pcd.graphs
import pcd.support.algorithms

g = pcd.graphs.karate_club(weighted=False)
gw = pcd.graphs.karate_club(weighted=True)
assert networkx.is_isomorphic(g, gw)


#G = pcd.Graph.from_networkx(g)
#G.enableModularity()

#a = pcd.support.algorithms.PCDmodSA(g, quiet=False, betafactor=1.001)
#print a.modularity
## 0.41978961209730375 - what my code returns
## 0.4197896           - Most precise one I see in literature
#assert abs(a.modularity - 0.4197896) < 1e-6
#assert sorted(a.cmtys.cmtysizes().values()) == [5, 6, 11, 12]



g_dol = pcd.graphs.dolphins()
for a,b,d in g_dol.edges_iter(data=True):
    d.clear()

a_dol = pcd.support.algorithms.PCDmodSA(g_dol, quiet=False, betafactor=1.0001)
print a_dol.modularity

from fitz import interactnow
