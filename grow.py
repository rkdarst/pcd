
import random

class Grower(object):
    def __init__(self, g, p, seed=None):
        if g is None:
            g = self.graph(g)
        self.g = g
        self.nodes = sorted(g.nodes())
        # In each growing step, prob. to attach to a random node in
        # the network.  Prob 1-p to attach to a neighbor of one of the
        # previous nodes in the addition sequence.
        self.p = 1
        self.rng = random.Random(seed) # seed not given: random seed.

    def add(self, n0=None):
        if n0 is None:
            n0 = len(self.nodes)  # one greater than greatest node
        # first neighbor
        n1 = self.rng.choice(self.nodes)
        self.g.add_edge(n0, n1)
        neighs = set()
        neighs.update(self.g.neighbors(n1))

        if self.rng.uniform(0, 1) < p:
            # p chance of making second link to a completly random node:
            while n2 is None or n2 in 
            n2 = random.choice(self.nodes)

        else:
            # 1-p chance of next link to a neighbor of node n1
            n1_neighs = 
