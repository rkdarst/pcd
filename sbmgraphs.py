# Richard Darst, July 2012

import collections
import math
import random
from scipy.stats import binom

import networkx

from pcd import nxutil

product = lambda l: reduce(lambda x,y: x*y, l, 1)


#def partition_nodes(U, q, n):
#    """Partition a set of nodes
#
#    Given a set of nodes, return min(N//Nsets, nodesPerSet) equal
#    partitions of the nodes.  This function does allow overlaps.
#
#    Algorithm: Divide U into U//q (approximately) equal size
#    communities first, spanning the system.  Then, choose
#    """
#    U = list(U)
#    random.shuffle(U)
#    N = len(U)//q
#    N = min(N, n)
#    sets = [ U[N*i:N*(i+1)] for i in range(q) ]
#    print "  ", [len(s) for s in sets]
#    for s in sets:
#        while len(s) < n:
#            node = random.choice(U)
#            if node in s:
#                continue
#            s.append(node)
#    print "  ", [len(s) for s in sets]
#    return sets

def makeCommunities(U, q, n):
    """Take a universe of nodes U and partition it into q communities of
    n, allowing overlaps."""
    U = set(U)
    Ulist = list(U)
    random.shuffle(Ulist)
    communities = [ ]

    # Make initial communities about equally sized, picking from all
    n_initial = int(math.ceil(len(U)/float(q)))
    n_initial = len(U)/float(q)
    if isinstance(n, (int, float)):
        assert n_initial <= n, "let(U)//q < n."
    else:
        raise NotImplementedError
        assert n_initial <= min(n) # so our initial community assignments work
        assert sum(n) > U

    for i in range(q):
        nodes = Ulist[int(round(n_initial*i)):int(round(n_initial*(i+1)))]
        communities.append(set(nodes))
    print "initial:", [ len(nodes) for nodes in communities ]
    assert 0 == len(U - set.union(*communities)) # every node is at least once
    assert len(U) == sum(len(n) for n in communities) # no nodes are duplicated
    assert sum(len(_) for _ in communities) == len(U)

    # Now pick random nodes among all other nodes to top-off each
    # community t, 'n' nodes.
    for nodes in communities:
        nodes |= set(random.sample(U-nodes, n-len(nodes)))
    print "final community sizes:", [ len(nodes) for nodes in communities ]

    return communities


def random_sample(nodes, n):
    """Return a random sample of size n."""
    nodes = list(nodes)
    return set(random.sample(nodes, n))

def random_sample_q(nodes, q, n):
    """Return list of q random samples of size n."""
    nodes = list(nodes)
    return [ set(random.sample(nodes, n)) for _ in range(q) ]

def add_edges(g, nodes, p, weight=None):
    """Add edges between """
    nodes = list(nodes)
    for i, x in enumerate(nodes):
        for j, y in enumerate(nodes[i+1:]):
            if weight:
                #g.add_edge(x, y, weight=weight)
                g.add_edge(a, b, weight=None)
                g.edge[a][b].setdefault('weights', []).append(p)
            elif random.uniform(0, 1) < p:
                g.add_edge(x, y)
def add_edges_cm(g, nodes, p, weighted=None):
    if weighted: raise NotImplementedError("add_edges_cm + weighted")
    nodes = list(nodes)
    nNodes = len(nodes)
    #nEdges = int(round(.5 * len(nodes) * (len(nodes)-1) * p))
    nEdges = binom(.5*len(nodes)*(len(nodes)-1), p).rvs()
    vertices = [ ]
    for n in nodes:
        vertices.extend([n] * int(round(2*nEdges/float(nNodes))))
    #if len(vertices)
    random.shuffle(vertices)
    #for a, b in zip(vertices[0::2], vertices[1::2]):
    #    if g.has_edge(a, b):
    #        ...
    #    g.add_edge(a, b)
    tries = 0
    while True:
        if len(vertices) < 2:
            break
        if len(vertices) == 0:
            break
        a = vertices.pop(0)
        b = vertices.pop(0)
        if g.has_edge(a, b):
            if tries > 10:
                break
            vertices.insert(random.randint(0, len(vertices)), a)
            vertices.insert(random.randint(0, len(vertices)), b)
            #if tries > 0: print vertices
            tries += 1
            continue
        g.add_edge(a, b)
        tries = 0
    print p, len(nodes), len(g.edges())
    #raw_input()
def add_edges_exact(g, nodes, p, weighted=None,
                    links=None):
    if not links:
        links = set(frozenset((a, b))
                    for a in nodes for b in nodes
                    if a != b)
        assert len(links) == len(nodes) * (len(nodes)-1) / 2
    else:
        assert not nodes
    if weighted == 'exact':
        #for x in links:
        #    print x
        #    a, b = x
        e = 0
        for a,b in links:
            g.add_edge(a, b, weight=None)
            g.edge[a][b].setdefault('weights', []).append(p)
            assert g.edge[b][a]['weights'] == g.edge[a][b]['weights']
            #print g.edge[a][b]['weight']
            e += 1
        print p, len(nodes), len(g.edges()), e
        #raw_input()
        return
    links = list(links)
    #nNodes = len(nodes)
    #nEdges = int(round(len(links) * p))
    nEdges = binom(len(links), p).rvs()
    edges = random.sample(links, nEdges)
    e = 0
    for a,b in edges:
        g.add_edge(a, b)
        e += 1
    #print p, len(nodes), len(g.edges()), e

def sbm_incomplete(U, q, n, p, weighted=None,
                   edges_constructor=None):
    """Create a SBM graph from a universe of nodes, possibly
    excluding some of the nodes."""
    U = set(U)
    g = networkx.Graph()

    # Pick our communities from the universe.  Note that the size of
    # the universe decreases.
    # is n constant, or does it change for all communities?
    if isinstance(n, (int, float)):
        communities = [ random.sample(U, n) for _ in range(q) ]
    else:
        communities = [ random.sample(U, int(round(_))) for _ in n ]
    U = set().union(*communities)
    print sum(len(x) for x in communities), [ len(x) for x in communities ]
    print len(U)

    # Add all nodes from the (new) universe, then add all of their edges.
    g.add_nodes_from(U)
    # is p constant for all communities, or per-community?
    if isinstance(p, (int, float)):
        for nodes in communities:
            edges_constructor(g, nodes, p, weighted=weighted)
    else:
        for nodes, p in zip(communities, p):
            edges_constructor(g, nodes, p, weighted=weighted)

    if weighted:
        for a, b in g.edges_iter():
            newweight = 1. - product(1.-x for x in g.edge[a][b]['weights'])
            g.edge[a][b]['weight'] = newweight

    # Set the planted communities:
    nxutil.cmtyInit(g)
    for c, nodes in enumerate(communities):
        nxutil.cmtyAddFromList(g, c, nodes)
    return g

def sbm(U, q, n, p, weighted=None, edges_constructor=add_edges_exact):
    U = set(U)
    g = networkx.Graph()

    # Pick our communities from the universe.
    communities = makeCommunities(U, q, n)

    # Add all nodes from the (new) universe, then add all of their edges.
    g.add_nodes_from(U)
    for nodes in communities:
        #add_edges(g, nodes, p, weighted=weighted)
        #add_edges_cm(g, nodes, p, weighted=weighted)
        #add_edges_exact(g, nodes, p, weighted=weighted)
        edges_constructor(g, nodes, p, weighted=weighted)

    if weighted:
        for a, b in g.edges_iter():
            newweight = 1. - product(1.-x for x in g.edge[a][b]['weights'])
            g.edge[a][b]['weight'] = newweight

    # Set the planted communities:
    nxutil.cmtyInit(g)
    for c, nodes in enumerate(communities):
        nxutil.cmtyAddFromList(g, c, nodes)
    return g

def compose_graphs(gs, weighted=None):
    g = networkx.Graph()
    for g_ in gs: g.add_nodes_from(g_.nodes_iter())
    for g_ in gs: g.add_edges_from(g_.edges_iter())
    if weighted:
        for a, b in g.edges_iter():
            #weights = [ _g.edge[a][b]['weight']
            #            for _g in gs if _g.has_edge(a, b) ]
            weights = [ ]
            [ weights.extend(_g.edge[a][b]['weights'])
              for _g in gs if _g.has_edge(a, b) ]
            newweight = 1. - product(1.-x for x in weights)
            #if len(weights) > 1:
            #    print weights, newweight
            g[a][b]['weight'] = newweight
    #sys.exit()
    #print gs
    for _g in gs:
        nxutil.graphproperties(_g )
    #raw_input('>')
    return g

def add_edges_out(g, p, g_layers, weighted=False,
                  edges_constructor=add_edges_exact):
    """Adds edges between any links _not_ in the same community."""
    nodes = set(g.nodes())
    if weighted: raise NotImplementedError("weighted+add_pout")
    links = set()
    _iterCmtys = nxutil._iterCmtys
    for n1, n1data in g.nodes_iter(data=True):
        for n2, n2data in g.nodes_iter(data=True):
            if n2 <= n1: continue
            # If there is no overlap
            if [1 for g_ in g_layers if
                  set(_iterCmtys(g_.node[n1]))&set(_iterCmtys(g_.node[n2]))]:
                continue
            links.add((n1, n2))
    edges_constructor(g, nodes=None, p=p, weighted=weighted,
                      links=links)

def multiLayerSBM(N_U, config, weighted=None, incomplete0=False,
                  edges_constructor=add_edges_exact,
                  pOut=None):
    """Config is a list of (q, n, p) tuples."""
    U = range(N_U)
    subgraphs = [ ]
    for level, (q, n, p) in enumerate(config):
        if level==0 and incomplete0:
            _g = sbm_incomplete(U, q=q, n=n, p=p, weighted=weighted,
                                edges_constructor=edges_constructor)
            U = tuple(_g.edges())
        else:
            _g = sbm(U, q=q, n=n, p=p, weighted=weighted,
                     edges_constructor=edges_constructor)
        subgraphs.append(_g)
    g = compose_graphs(subgraphs, weighted=weighted)
    if pOut:
        add_edges_out(g, p=pOut, g_layers=subgraphs, weighted=weighted,
                      edges_constructor=edges_constructor)
    #for a,b in  g.edges_iter(): assert _g.has_edge(a, b)
    #for a,b in _g.edges_iter(): assert  g.has_edge(a, b)
    #for a,b in _g.edges_iter():
    #    assert g.edge[a][b]['weight'] == _g.edge[a][b]['weight']
    return g, subgraphs

def overlap_communities(g, n_min=1):
    """Return graph with communities set to all pairwise overlap of
    given graph communities.

    n_min = minimum size of overlap to be returned."""
    g_new = g.copy()
    nxutil.cmtyInit(g_new)

    cmtys = nxutil.communities(g)
    cmtys_list = list(cmtys)

    overlap_cmtys = { }
    overlap_cmty_i = 0
    print cmtys_list
    for i, c1 in enumerate(cmtys_list):
        c1nodes = cmtys[c1]
        for j, c2 in enumerate(cmtys_list[i+1:]):
            c2nodes = cmtys[c2]
            newCmtyNodes = c1nodes & c2nodes
            print c1, c2, len(c1nodes), len(c2nodes), len(newCmtyNodes)
            if len(newCmtyNodes) > n_min:
                overlap_cmtys[overlap_cmty_i] = newCmtyNodes
                overlap_cmty_i += 1
                nxutil.cmtyAddFromList(g_new, overlap_cmty_i, newCmtyNodes)
    return g_new


def make_overlap_test(N_U, q1, n1, p1, q2, n2, p2, pU):
    """

    N_U: number of nodes in the universe.
    q1, q2: sizes of level1 and level2 communities.
    n1, n2: number of nodes in level1 and level2 communities.
    p1, p2: edge densities
    pU: background connection probability.
    """
    U = set(range(N_U))
    g = networkx.Graph()

    #n1 = 25
    #q1 = 40
    #p1 = .95
    level1 = [ random.sample(U, n1) for _ in range(q1) ]
    U = set().union(*level1)
    print sum(len(x) for x in level1), [ len(x) for x in level1 ]
    print len(U)

    g.add_nodes_from(U)

    for points in level1:
        add_edges(g, points, p1)

    #level2 = None
    #n2 = 100
    #q2 = 10
    #p2 = .65
    #level2 = [ random.sample(U, size2) for _ in range(q2) ]
    level2 = partition_nodes(U, q2, n2)

    for points in level2:
        add_edges(g, points, p2)

    add_edges(g, list(g.nodes()), pU)

    g1 = g.copy()
    nxutil.cmtyInit(g1)
    for c, nodes in enumerate(level1):
        nxutil.cmtyAddFromList(g1, c, nodes)


    g2 = g.copy()
    nxutil.cmtyInit(g2)
    for c, nodes in enumerate(level2):
        nxutil.cmtyAddFromList(g2, c, nodes)

    return g, g1, g2

