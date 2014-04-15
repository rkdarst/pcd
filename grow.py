"""Growing graph models.
"""


import networkx

from math import exp, pow

import random

import pcd.util

class GrowingGraph(object):
    def grow(self, N):
        """Add enough nodes to graph to reach N nodes."""
        while len(self.g) < N:
            #print len(self.g)
            self.add()

    @classmethod
    def get(cls, **kwargs):
        N = kwargs.pop('N')
        grower = cls(**kwargs)
        grower.grow(N)
        return grower.g


class GrowFitness(GrowingGraph):
    """Growing graph model, with power-law attachment

    neighbor_mode: how to add edges beyond the second
        all: neighbors of all previous nodes attached to
        last: neighbors of last node attached to
        first: neighbors of first node attached to.
        all_nonrandom: like all, but exclude random nodes (from p)
        last_nonrandom: like last, but exclude random nodes (from p)
    """
    neighbor_mode = 'all'
    track_edges = False
    def __init__(self, g, p, beta, kappa, m=2, seed=None,
                 **kwargs):
        # In each growing step, prob. to attach to a random node in
        # the network.  Prob 1-p to attach to a neighbor of one of the
        # previous nodes in the addition sequence.
        self.p = p
        self.beta = beta
        self.kappa = kappa
        self.m = m
        self.rng = random.Random(seed) # seed not given: random seed.
        for k, v in kwargs.iteritems():
            assert hasattr(self, k), "GrowFitness doesn't have attr %s"%k
            setattr(self, k, v)

        # Set up initial graph
        if g is None:
            g = networkx.complete_graph(self.m+1)
            self.assign_initial_fitnesses(g)
        self.g = g

        # Set up fitness lists
        self.fitnesses = dict((n, g.node[n]['fitness']) for n in g.nodes_iter() )

        self.chooser = pcd.util.WeightedChoice(self.fitnesses.iteritems(), rng=self.rng)

    def assign_initial_fitnesses(self, g):
        for n in g.nodes_iter():
            g.node[n]['fitness'] = \
                  exp(-self.beta*pow(self.rng.uniform(0, 1),1./(self.kappa+1)))

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

        # Choose the first node to link to.
        #chooser = pcd.util.WeightedChoice(self.fitnesses.iteritems())
        while True:
            n1 = self.chooser.choice()
            if n1 != n0:
                break
        g.add_edge(n0, n1)
        # Current neighbors of linked nodes
        neighs = set(g.adj[n1])
        # Existing links.  Do *not* make new links to these.
        links_exclude = set((n0,))
        links_exclude.add(n1)

        # Track edge creation, if requested.
        if self.track_edges:
            g.node[n0]['edge_order'] = [n1]

        for _ in range(self.m - 1):

            # NOTE: the Communities.c file seems to be backwards here,
            # with a probability of 1-p instead of p.  The real
            # situation needs to be worked out...
            random_last = False
            if self.rng.uniform(0, 1) < self.p:
                # p chance of making second link to a completly random node:
                random_last = True
                assert len(self.chooser) == len(g)
                while True:
                    n_next = self.chooser.choice()
                    if n_next not in links_exclude:
                        break

            else:
                # 1-p chance of next link to a neighbor of node n1
                assert neighs - links_exclude != 0, "We have no remaining neighoring nodes to connect to"
                #print n0, len(g), len(neighs), len(links_exclude)
                #if len(neighs) < .05 * len(g):
                if len(neighs) < 100:
                    neigh_fitnesses = [ (n, self.fitnesses[n]) for n in neighs ]
                    neighbor_chooser = pcd.util.WeightedChoice(neigh_fitnesses, rng=self.rng)
                    # Choose the next node
                    while True:
                        n_next = neighbor_chooser.choice()
                        if n_next not in links_exclude:
                            break
                else:
                    #print 'global'
                    while True:
                        n_next = self.chooser.choice()
                        if n_next in neighs and n_next not in links_exclude:
                            break

            assert not g.has_edge(n0, n_next)
            g.add_edge(n0, n_next)
            links_exclude.add(n_next)
            # Update our neighbor lists depending on how we do the
            # addition:
            if self.neighbor_mode == 'all':
                neighs.update(g.adj[n_next])
            elif self.neighbor_mode == 'all_nonrandom':
                if not random_last:
                    neighs.update(g.adj[n_next])
            elif self.neighbor_mode == 'last':
                neighs = set(g.adj[n_next])
            elif self.neighbor_mode == 'last_nonrandom':
                if not random_last:
                    neighs = set(g.adj[n_next])
            elif self.neighbor_mode == 'first':
                pass
            else:
                raise ValueError('Unknown neighbor_mode %s'%self.neighbor_mode)
            # Track edge growth if desired:
            if self.track_edges:
                g.node[n0]['edge_order'].append(n_next)


class HolmeGraph(GrowingGraph):
    def __init__(self, m, mt=None, Pt=None, m0=3, seed=None):
        """

        m: int
            additional edges for each successive node
        Pt: float
            Probability of triad formation.
        mt: float
            Avg number of clustered edges per run.  Pt is then set to mt/(m-1)
        m0: int
             original number of nodes.
        """
        self.g = networkx.Graph()
        self.m = m
        m0 = max(m0, m+1)
        self.m0 = m0
        self.rng = random.Random(seed) # seed not given: random seed.
        # We can specify eather Pt or Mt
        if Pt is None:
            Pt = mt/float(m-1)
        else:
            assert mt is None
        assert Pt <= 1.0, "Pt=%s must be <= 1.0"%Pt
        self.Pt = Pt
        # Add initial nodes
        for n in range(m0):
            self.g.add_node(n)  # One free unit for each node.
        # Preferential attachment selector.  Selects random nodes
        # proportional to degree.
        self.pa_selector = list(range(m0))
    def add(self):
        g = self.g
        n = len(g)
        pa_selector = self.pa_selector
        rng = self.rng

        g.add_node(n)
        pa_selector.append(n)  # One free unit for each node.

        def pa_select():
            if len(g) <= len(edges_made):
                return None
            while True:
                n1 = rng.choice(pa_selector)
                if n1 in edges_made:
                    continue
                break
            return n1

        edges_made = set((n, ))
        assert len(g) > self.m, "Graph to small and not enough edges to connect to."
        for e in range(self.m):
            # Random edge
            n1 = pa_select()
            if n1 is None:
                print "  can't make more edges at %d"%n
                break
            assert not g.has_edge(n, n1)
            assert n != n1
            g.add_edge(n, n1)
            edges_made.add(n1)
            pa_selector.extend((n, n1))

            if rng.uniform(0,1) <= self.Pt:
                # Next, add clustered edge
                second_neighbors = set(g.neighbors(n1))
                available = second_neighbors - edges_made
                # PA or triad closure step?
                if available:
                    # Triad closure
                    n2 = rng.sample(available, 1)[0]
                else:
                    # No second neighbor available.  Do a random PA step
                    # instead.
                    n2 = pa_select()
                    if n2 is None:
                        print "    can't make more edges at %d, %d"%(n, n1)
                        continue
                assert not g.has_edge(n, n2)
                assert n != n2
                assert n1 != n2

                g.add_edge(n, n2)
                edges_made.add(n2)
                pa_selector.extend((n, n2))


class SquareClosure(GrowingGraph):
    def __init__(self, ns, seed=None, m0=3):
        self.g = networkx.Graph()
        self.ns = ns
        self.rng = random.Random(seed) # seed not given: random seed.
        for n in range(m0):
            self.g.add_node(n)
            for n1 in range(n):
                self.g.add_edge(n, n1)
    def add(self):
        g = self.g
        rng = self.rng
        n = len(g)
        g.add_node(n)

        # Random base of new attachment
        n1 = random.randint(0, len(g)-2)  # exclude endpoint and n
        #print 'graph status:', len(g), g.number_of_edges(), n, n1#, g.nodes()
        g.add_edge(n, n1)

        inner_nodes = set((n, n1))
        shell_nodes = inner_nodes
        for dist, n_at_dist in enumerate(self.ns):
            dist += 1  # dist starts from 1, not zero
            nodes_available = set.union(*(set(g.neighbors(_))
                                          for _ in shell_nodes))
            #print '  nodes_available:', nodes_available
            #print '  inner_nodes:', inner_nodes
            shell_nodes = nodes_available - inner_nodes
            #print '  shell_nodes:', shell_nodes
            #print ' ', dist, n_at_dist, len(inner_nodes), len(nodes_available)
            if len(shell_nodes) == 0:
                # We can never add any more from now on.
                break

            n_at_dist = min(n_at_dist, len(shell_nodes))
            if n_at_dist == 0:
                inner_nodes |= shell_nodes
                continue

            new_edges = rng.sample(shell_nodes, n_at_dist)
            #print '  new edges:', new_edges, len(shell_nodes)
            for n2 in new_edges:
                #print '    linking to:', n2
                if n2 == n or n2 == n1: raise
                assert not g.has_edge(n, n2)
                g.add_edge(n, n2)

            inner_nodes |= nodes_available



def growsf_gb(N, p, beta, kappa, m=2, **kwargs):
    """Make a gb grown graph."""
    grower = GrowFitness(p=p, beta=beta, kappa=kappa, m=m, g=None,
                         **kwargs)

    grower.grow(N)
    assert len(grower.g) == N
    assert grower.g.number_of_edges() == (N*m - ((m+1)**2 -(m+1))/2)
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
    #print len(GrowBA.create(10000, m=2))
    #print len(networkx.barabasi_albert_graph(n=10000, m=2))
    growsf_gb(N=100000, p=0, beta=20, kappa=6, m=5)
