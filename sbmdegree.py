# Richard Darst, September 2013

"""Degree-corrected SBM graphs, as in Karer and Newman"""

import networkx

import sbmgraphs
import nxutil
import numpy

def graph(q, n, distribution, lambda_=1):
    N = q*n
    U = range(N)
    g = networkx.Graph()
    g.add_nodes_from(U)

    node_degrees = dict((n, distribution.rv()) for n in U)
    # Pick our communities from the universe.
    communities = sbmgraphs.makeCommunities(U, q, n)
    community_degrees = dict((q_, sum(node_degrees[n] for n in ns))
                             for q_, ns in enumerate(communities))

    # Set the planted communities:
    nxutil.cmtyInit(g)
    for c, nodes in enumerate(communities):
        nxutil.cmtyAddFromList(g, c, nodes)

    # planted configuration
    omega_planted = numpy.zeros((q, q))
    for q_, ns in enumerate(communities):
        omega_planted[q_,q_] = sum(node_degrees[n] for n in ns)
    # random configuration
    E = sum(node_degrees.itervalues()) / 2.
    omega_random = numpy.zeros((q, q))
    for q1 in range(q):
        for q2 in range(q):
            omega_random[q1, q2] = (community_degrees[q1]*community_degrees[q2]
                                    / float(2. * E))

    # Mix the random and planted partitions.
    omega = lambda_*omega_planted + (1-lambda_)*omega_random

    from pcd.util import WeightedChoice
    node_choosers = [ WeightedChoice(  (n, node_degrees[n]) for n in communities[q_]  )
                      for q_ in range(q) ]
    # Add all the edges
    redundant_edges = 0
    self_edges = 0
    for q1 in range(q):
        for q2 in range(q1, q):
            expected = omega[q1, q2]
            if q1 == q2:
                expected = .5 * omega[q1, q2]
            for i in range(int(round(expected))):
                n1 = node_choosers[q1].choice()
                n2 = node_choosers[q2].choice()
                if n1 == n2:
                    self_edges += 1
                    continue
                if g.has_edge(n1, n2):
                    redundant_edges += 1
                    continue
                g.add_edge(n1, n2)

    # checks:
    print sum(g.degree().values()) + redundant_edges*2 + self_edges*2, sum(node_degrees.values())


    #from fitz import interactnow
    return g

def _test1():
    from pcd.support.powerlaw import PowerLaw
    dist = PowerLaw(-1, xmin=10, xmax=30)
    g = graph(2, 1000, dist, lambda_=0)


if __name__ == "__main__":
    _test1()
