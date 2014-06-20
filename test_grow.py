from nose.tools import assert_true, assert_equal, assert_greater_equal, assert_less, assert_greater

from pcd.grow import *

from pcd.support.testutil import assert_isomorphic

G = networkx.Graph

import networkx


def test_holme():

    for N in (5, 10, 15):
        for m in (1, 2, 3, 4):
            m0 = max(3, m)    # must reproduce logic of the model.
            g = HolmeGraph.get(N=N, m=m, Pt=1, m0=m0)
            assert_equal(g.number_of_edges(),   (N-m0)*m)

    for N in (5, 10, 15):
        for m in (1, 2, 3, 4):
            m0 = max(3, m)    # must reproduce logic of the model.
            g = HolmeGraph.get(N=N, m=m, Pt=.5, m0=m0)
            assert_equal(g.number_of_edges(),   (N-m0)*m)


    # Should return itself
    g0 = networkx.complete_graph(5)
    g = HolmeGraph.get(N=5, m=3, Pt=1, g0=g0)
    assert_equal(networkx.triangles(g), dict((i,4*3/2) for i in range(5)))

    # Test number of triangles created for new nodes.
    g0 = networkx.complete_graph(5)
    g = HolmeGraph.get(N=6, m=3, Pt=1, g0=g0)
    assert_equal(networkx.triangles(g)[5], 3)

    g0 = networkx.complete_graph(6)
    g = HolmeGraph.get(N=7, m=3, Pt=1, g0=g0)
    assert_equal(networkx.triangles(g)[6], 3)

    # Test number of triangles created for new nodes.
    def _make():
        g = HolmeGraph.get(N=6, m=3, m0=5, Pt=0)
        return networkx.triangles(g)
    sizes = [_make() for _ in range(10)]
    assert_true(any(_[5]==0 for _ in sizes))



def test_triad():
    for N in [5, 10, 20]:
        for m in [2, 3]:
            g = TriadGraph.get(N=N, m=m, p=0)
            assert_equal(g.number_of_edges(), 3+(N-3)*m)

    for N in [5, 10, 20]:
        for m in [2, 3]:
            g = TriadGraph.get(N=N, m=m, mmean=m+1, p=0)
            assert_greater(g.number_of_edges(), 3+(N-3)*m)

    g = TriadGraph.get(N=40, m=3, p=1)
    assert_greater(networkx.average_clustering(g), .5)

    g = TriadGraph.get(N=40, m=3, p=.0)
    assert_less(networkx.average_clustering(g), .5)
