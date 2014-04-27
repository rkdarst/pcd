import pcd
import pcd.graphs
import pcd.cmty

def approxeq(a, b, tol=1e-6):
    return abs(a-b)/float(abs(a)+abs(b)) < tol



def dotest_mod(g, known_mod=None):
    G = pcd.Graph.fromNetworkX(g)
    G.enableModularity(None)

    G.verbosity = -1
    if known_mod:
        #G.anneal(1.0, Escale=100, betafactor=1.00001, p_binary=0.05, p_collective=0)
        G.trials(1.0, 10, minimizer='anneal',
                 Escale=100, betafactor=1.001, p_binary=0.05, p_collective=0)
    else:
        G.anneal(1.0)

    cmtys = pcd.cmty.Communities.from_pcd(G)

    #print cmtys.q
    #print repr(cmtys.Q(g))
    print repr(cmtys._Q_cmty(g))
    #print repr(cmtys._Q_pcd(g))
    #print repr(G.modularity())
    print repr(cmtys._Q_pcd_energy(g))

    assert approxeq(cmtys.Q(g), cmtys._Q_cmty(g))
    assert approxeq(cmtys.Q(g), cmtys._Q_pcd(g))
    assert approxeq(cmtys.Q(g), G.modularity())
    #assert approxeq(cmtys.Q(g), cmtys._Q_pcd_energy(g))

    if known_mod:
        assert approxeq(cmtys.Q(g), known_mod, .001), (cmtys.Q(g), known_mod)

def test():
    g = pcd.graphs.karate_club()
    dotest_mod(g)
    # known_mod test is disabled for now, since annealing is
    # non-deterministic and it doesn't work every time.
    # karate club Q is Q=0.4197896 (to 6 decimal places)
    #test_mod(g, known_mod=0.419789)

    #g = pcd.graphs.dolphins()
    #test_mod(g)
