
from nose.tools import assert_true, assert_equal, assert_greater_equal, assert_less

from pcd.support.growsf_gb import *


from pcd.support.testutil import assert_isomorphic

import networkx
G = networkx.Graph


def test_sole():
    # For small graphs we can exactly specify what the outcome should be:
    # alpha=0, delta=0
    assert_isomorphic(sole(T=3, alpha=0, delta=0),
                      G({0:(1,2), 1:(0,2)}))

    assert_isomorphic(sole(T=4, alpha=0, delta=0),
                      G({0:(1,2), 1:(0,2), 3:(1,2)}))

    # T=5 is harder but the graph should match ONE of these!
    for g in [G({0:(1,2), 1:(0,2), 3:(1,2), 4:(1,2),   }),
              G({0:(1,2), 1:(0,2), 3:(1,2), 4:(0,2,3), }),
              ]:
        if networkx.is_isomorphic(sole(T=5, alpha=0, delta=0), g):
            'matches x'
            break
    else:
        raise AssertionError("One did not match.")

    # With alpha=1 it should be a clique
    # alpha=1, delta=0
    for T in [3, 4, 7, 10, 15]:
        assert_isomorphic(sole(T=T, alpha=1, delta=0),
                          networkx.complete_graph(T),
                          "Clique of size %d"%10)

    # With alpha=0, delta=0 we should have at least two edges for each new node
    # alpha=0, delta=0
    for T in [3, 4, 7, 10, 25]:
        g = sole(T=T, alpha=0, delta=0)
        assert_greater_equal(g.number_of_edges(), 2*T-3)

    # Test duplication and keeping only one of the links.
    # alpha=1, delta=1
    for T in [3, 4, 5, 10, 15, 25]:
        g = sole(T=T, alpha=1, delta=1)
        assert_equal(g.number_of_edges(), T)

    # Test Keeping only one of the links.  No link to original.  This graph is disconnected.
    # alpha=0, delta=1
    for T in [3, 4, 5, 10, 15, 25]:
        g = sole(T=T, alpha=0, delta=1)
        assert_equal(g.number_of_edges(), 3)


def test_marsili():
    # These should generate clique graphs.
    for N in [5, 10, 20]:
        g = marsili(N=N, T=2000, lambda_=0, xi=0)
        #assert_isomorphic(g, networkx.complete_graph(N))
        assert_greater_equal(g.number_of_edges(), N*(N-1)/2 - 2)

    # This should also be clique graphs, but forcing to only triadic closure.
    for N in [5, 10, 20]:
        print N
        while g.number_of_nodes() < N:
            g = marsili(N=N, T=2000, lambda_=0, xi=100000)
        #assert_isomorphic(g, networkx.complete_graph(N))
        #assert_equal(g.number_of_edges(), N*(N-1)/2 - 2)
        assert_greater_equal(g.number_of_edges(), N*(N-1)/2 - 2)
        #E = g.number_of_edges()
        #assert_true(E in ((N*(N-1))/2,      (N*(N-1))/2-1,
        #                  ((N-1)*(N-2))/2,  ((N-1)*(N-2))/2-1,  ))

    # This should be an empty graph.
    for N in [5, 10, 20]:
        print N
        g = marsili(N=N, T=2000, lambda_=10000, xi=0)
        assert_less(g.number_of_edges(), 3, msg="Almost empty graph.")
