# Richard Darst, September 2011

from util import color_graph

gamma=1
from graphs import random_graph
G = random_graph(size=20)



G.minimize(gamma=gamma)
g2 = G.supernode_networkx()

color_graph(g2)

#G.viz()
