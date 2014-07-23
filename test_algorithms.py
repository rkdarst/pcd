
import pcd.support.algorithms as algs
import pcd.graphs
import pcd.cmty

import networkx

from nose.plugins.attrib import attr
from nose.tools import assert_set_equal, assert_true

# Methods to NOT test.
methods_exclude = set(('CDMethod',
                       'BeliefPropogationZ',
                       'BeliefPropogationq',
                       'PCDmodSA',   # works but slow.
                       ))

# These methods don't pass the two-clique test
methods_dontwork = set(('SnapAGMfit',
                        'NullCD',
                        'SnapBigClam',
                        'BeliefPropagationZ',
                        'LouvainWeighted',
                        'FlowSpectrum',
                        'LinkCommunities',
                        'Infomap_dir',
                        'InfomapSingle_dir',
                        ))


methods = sorted(
    name for name, cda in vars(algs).iteritems()
    if isinstance(cda, type)
    and issubclass(cda, algs.CDMethod)
    and name[0]!='_')

methods = [x for x in methods if
           'APM' not in x]

method_args = {
    'BeliefPropogationZ':dict(q=2),
    }

default_args = dict(verbosity=0)

# karate club
g_karate = pcd.graphs.karate_club()

@attr('slow')
def test_algorithms():
    """Just test that methods run successfully."""

    for method in methods:
        if method in methods_exclude:
            continue
        g = g_karate.copy()
        yield do_run, method, g

def do_run(method, g):
    print method
    cda = getattr(algs, method)
    args = default_args.copy()
    args.update(method_args.get(method, {}))
    r = cda(g, **args)

    assert_true(hasattr(r, 'cmtys'))
    assert_true(hasattr(r, 'results'))


nclq = 10  # clique size
g1 = networkx.complete_graph(nclq)
g2 = networkx.relabel_nodes(g1, dict(zip(range(0,nclq),range(nclq,nclq*2))), copy=True)
g_2clique = networkx.union(g1, g2)
g_2clique.add_edge(0, nclq)
g_2clique_cmtys = pcd.cmty.Communities({0:set(range(0,nclq)), 1:set(range(nclq,nclq*2))})


@attr('slow')
def test_2clique():
    """Test that algorithms run and give reasonable results."""

    for method in methods:
        if method in methods_exclude:
            continue
        if method in methods_dontwork:
            continue
        g = g_2clique.copy()
        yield do_method, method, g_2clique, g_2clique_cmtys

def do_method(method, g, cmtys_true):
    for i in range(10):
        try:
            return do_method_trial(method, g, cmtys_true)
        except AssertionError:
            continue
    raise

def do_method_trial(method, g, cmtys_true):
    print method
    cda = getattr(algs, method)
    args = default_args.copy()
    args.update(method_args.get(method, {}))
    r = cda(g, **args)

    cmtys = r.cmtys

    cmtys_dict   = set(frozenset(_) for _ in cmtys.itervalues())
    cmtys_dict_0 = set(frozenset(_) for _ in cmtys_true.itervalues())

    assert_set_equal(cmtys_dict, cmtys_dict_0)
