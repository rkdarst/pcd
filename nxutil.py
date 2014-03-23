# Richard Darst, July 2012

"""Module to handle networkx graphs and their interfaces.

Networkx graphs have two ways of representing the community:

- node attribute 'cmty' which indicates the respective community of
  the graph.
- node attribute 'cmtys' which is an iterable listing muultiple
  communities for the node.
"""

import collections
import numpy
import random

def edges_between(g, nbunch1, nbunch2):
    """Edges betwen nbunches.

    If nbunch1==nbunch2, return double the number of edges."""
    nbunch2 = set(nbunch2)
    edges = set()
    for n1 in g.nbunch_iter(nbunch1):
        for neighbor in g.neighbors_iter(n1):
            if neighbor in nbunch2:
                edges.add((n1, neighbor))
    return edges
def n_edges_between(g, nbunch1, nbunch2):
    """Number of edges between two node sets.

    nbunches should be SETS for efficiency.

    If nbunch1==nbunch2, return double the number of edges."""
    n = 0
    #nbunch2 = set(nbunch2)
    if len(nbunch1) > len(nbunch2):
        nbunch1, nbunch2 = nbunch2, nbunch1
    for n1 in nbunch1:
        for neighbor in g.adj[n1]:
            if neighbor in nbunch2:
                n += 1
    return n

def _iterCmtys(d):
    """Iterate over the communities in a node' data dict.

    This abstracts out whether or not there is one community or
    multiple communities for the node."""
    if 'cmty' in d:
        #raise
        if 'cmtys' in d:
            raise ValueError('Multiple "cmty" and "cmtys" both exist.')
        yield d['cmty']
    else:
        for c in d['cmtys']:
            yield c

def setCmtyAssignments(G, g):
    """Loads communities from a networkx graph g into a pcd graph G.

    Make sure you clear the G  first.  G.cmtyListClear()"""
    G.oneToOne = 0
    for node, data in g.nodes_iter(data=True):
        #print node, data['cmtys']
        for c in _iterCmtys(data):
            G.cmtyListAddOverlap(c, G._nodeIndex[node])
    return G

def cmtyInit(g):
    """Clear/initialize all communities in a networkx graph."""
    for n, d in g.nodes_iter(data=True):
        d.pop('cmty', 0)
        d['cmtys'] = set()
def cmtyClear(g):
    for n, d in g.nodes_iter(data=True):
        d.pop('cmty', 0)
        d.pop('cmtys', 0)
def cmtyAddFromList(g, cmtyID, nodes):
    """For a networkx graph, add all 'nodes' to 'cmtyID'."""
    #print cmtyID, nodes
    #for n in g.nbunch_iter(nodes):
    for n in nodes:
        g.node[n]['cmtys'].add(cmtyID)


def graphproperties(g):
    import collections
    cmty_d = collections.defaultdict(set)
    for node, data in g.nodes_iter(data=True):
        for c in _iterCmtys(data):
            cmty_d[c].add(node)
    print "number of communities:", len(cmty_d)
    print "mean cmty size:", numpy.mean([len(v) for v in cmty_d.values()])
    print "std cmty size:", numpy.std([len(v) for v in cmty_d.values()])
    U = set(g.nodes_iter())
    U_c = set()
    #print cmty_d
    #print [ v for v in cmty_d.itervalues() ]
    [ U_c.update(v) for v in cmty_d.itervalues() ]
    print "total nodes:", len(U)
    print "total nodes in communities:", len(U_c)


def communities(g):
    raise Exception
    """Given a graph, return a mapping c -> nodes_in_c."""
    cmtys = collections.defaultdict(set)
    for node, data in g.nodes_iter(data=True):
        for c in _iterCmtys(data):
            cmtys[c].add(node)
    return cmtys
cmtynodes = communities
def nodecmtys(g):
    raise Exception
    """Given a graph, return mapping n -> its communities"""
    cmtys = collections.defaultdict(set)
    for node, data in g.nodes_iter(data=True):
        for c in _iterCmtys(data):
            cmtys[node].add(c)
    return cmtys



def graphcolor(g, colormap=None, distance=1):
    """Return a coloring of a graph.

    Arguments:

    g: undirected networkx.Graph
        graph to be colored
    colormap: str, matplotlib colormap name.
        If given, should be a string listing a colormap name.  Instead
        of returning node->int map, return node->color tuple map.
    distance: int, default 1
        Nodes have unique colors for the n-th nearest neighbors,
        instead of only nearest neighbors.

    Returns:

    dict node->int
        indexes for colors (or color tuples).
    """
    colors = set()
    coloring = { }
    # Go through and make a greedy coloring.  Coloring-satisfying
    # integer for each node.
    for n in g:
        # 1-st nearest neighbors
        neighbors = set(g.neighbors(n))
        # n-th nearest neighbors
        for _ in range(distance-1):
            neighbors = set(y for neigh in neighbors
                            for y in g.neighbors(neigh))
        #
        used_colors = set(coloring[neigh] for neigh in neighbors
                           if neigh in coloring)
        avail_colors = colors - used_colors
        if avail_colors:
            color = random.choice(list(avail_colors))
        else:
            color = len(colors)
            colors.add(color)
        coloring[n] = color
    # This step goes through and re-randomizes choices, given the
    # minimal set of colors we found already.
    for n in g:
        used_colors = set(coloring[neigh] for neigh in g.neighbors(n)
                           if neigh in coloring)
        avail_colors = colors - used_colors
        color = random.choice(list(avail_colors))
        coloring[n] = color

    # If wanted, do a matplotlib colormap
    if colormap:
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors

        colormap = cm.get_cmap(colormap)
        normmap = mcolors.Normalize(vmin=0, vmax=len(colors))
        coloring = dict((n, colormap(normmap(color)))
                          for n, color in coloring.iteritems())
    return coloring

