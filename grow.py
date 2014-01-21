"""Growing graph models.
"""


import networkx

from math import exp, pow

import random

import pcd.util

class GrowFitness(object):
    """Growing graph model, with power-law attachment"""
    def __init__(self, g, p, beta, kappa, m=2, seed=None):
        # In each growing step, prob. to attach to a random node in
        # the network.  Prob 1-p to attach to a neighbor of one of the
        # previous nodes in the addition sequence.
        self.p = p
        self.beta = beta
        self.kappa = kappa
        self.m = m
        self.rng = random.Random(seed) # seed not given: random seed.

        # Set up initial graph
        if g is None:
            g = networkx.Graph()
            self.add_three_nodes(g)
        self.g = g

        # Set up fitness lists
        self.fitnesses = dict((n, g.node[n]['fitness']) for n in g.nodes_iter() )

        self.chooser = pcd.util.WeightedChoice(self.fitnesses.iteritems())

    def add_three_nodes(self, g):
        g.add_edge(0, 1)
        g.add_edge(1, 2)
        g.add_edge(2, 0)
        g.node[0]['fitness'] = exp(-self.beta*pow(self.rng.uniform(0, 1),1./(self.kappa+1)))
        g.node[1]['fitness'] = exp(-self.beta*pow(self.rng.uniform(0, 1),1./(self.kappa+1)))
        g.node[2]['fitness'] = exp(-self.beta*pow(self.rng.uniform(0, 1),1./(self.kappa+1)))

    def add(self, n0=None):
        g = self.g
        if n0 is None:
            n0 = len(self.g)  # one greater than greatest node

        g.add_node(n0)

        # Generate and set fitnesses
        fitness = exp(-self.beta*pow(self.rng.uniform(0, 1), 1./(self.kappa+1)))
        self.fitnesses[n0] = fitness
        g.node[n0]['fitness'] = fitness

        self.chooser.add(n0, fitness)

        # Current neighbors of linked nodes
        neighs = set()
        # Existing links.  Do *not* make new links to these.
        links_exclude = set((n0,))


        # Choose the first node to link to.
        #chooser = pcd.util.WeightedChoice(self.fitnesses.iteritems())
        while True:
            n1 = self.chooser.choice()
            if n1 != n0:
                break
        g.add_edge(n0, n1)
        # Update neighbors:
        neighs.update(g.neighbors(n1))
        links_exclude.add(n1)

        for _ in range(self.m - 1):

            # NOTE: the Communities.c file seems to be backwards here,
            # with a probability of 1-p instead of p.  The real
            # situation needs to be worked out...
            if self.rng.uniform(0, 1) < self.p:
                # p chance of making second link to a completly random node:
                assert len(self.chooser) == len(g)
                while True:
                    n_next = self.chooser.choice()
                    if n_next not in links_exclude:
                        break
                assert not g.has_edge(n0, n_next)
                g.add_edge(n0, n_next)
                neighs.update(g.neighbors(n_next))
                links_exclude.add(n_next)

            else:
                # 1-p chance of next link to a neighbor of node n1
                neigh_fitnesses = [ (n, self.fitnesses[n]) for n in neighs ]
                neighbor_chooser = pcd.util.WeightedChoice(neigh_fitnesses)
                # Choose the next node
                while True:
                    n_next = neighbor_chooser.choice()
                    if n_next not in links_exclude:
                        break
                assert not g.has_edge(n0, n_next)
                g.add_edge(n0, n_next)
                neighs.update(g.neighbors(n_next))
                links_exclude.add(n_next)


    def grow(self, N):
        """Add enough nodes to graph to reach N nodes."""
        while len(self.g) < N:
            #print len(self.g)
            self.add()


def growsf_gb(N, p, beta, kappa):
    """Make a gb grown graph."""
    grower = GrowFitness(p=p, beta=beta, kappa=kappa, m=2, g=None)

    grower.grow(N)
    assert len(grower.g) == N
    assert grower.g.number_of_edges() == 2*N - 3
    return grower.g


class GrowBA(object):
    def __init__(self, m=1, g=None):
        self.m = m
        if g is None:
            g = networkx.Graph()
        self.g = g
        assert set(g.nodes()) == set(range(len(g)))
        #self.chooser = chooser = pcd.util.WeightedChoice((n, g.degree(n)) for
        #                                                 n in range(len(g)))
    def add(self, n0=None):
        g = self.g
        if n0 is None:
            n0 = len(self.g)  # one greater than greatest node
        assert n0 not in self.g.node

        if len(g) == 0:
            g.add_node(n0)
            return
        elif len(g) == 1:
            g.add_edge(n0, next(iter(g.nodes())))
            g.add_node(n0)
            return

        #chooser = self.chooser
        links = set()
        chooser = pcd.util.WeightedChoice((n, g.degree(n)**2)
                                          for n in g.nodes())

        g.add_node(n0)

        for i in range(min(self.m, len(g))):
            while True:
                n1 = chooser.choice()
                if n1 not in links:
                    break
            g.add_edge(n0, n1)
            links.add(n1)
        ## Add to chooser
        #chooser.add(n0, m)
        #for n1 in links:
        #    chooser[n1] += 1
        #chooser.norm += m
        #chooser.check()

    @classmethod
    def create(cls, N, m=1):
        grower = cls(m=m)
        for _ in range(N):
            grower.add()
        return grower.g

if __name__ == "__main__":
    print len(GrowBA.create(10000, m=2))
    #print len(networkx.barabasi_albert_graph(n=10000, m=2))
