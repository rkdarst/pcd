# Richard Darst, August 2011

"""Functions to return standard graphs.

There is an important distinction between the Graph class in
models.py, and NetworkX graphs.  Most 
"""

import numpy
import os
import random
from xml.etree import ElementTree

import networkx
from networkx.generators.classic import grid_graph
import networkx.generators.random_graphs

try:
    from networkx.relabel import relabel_nodes as nx_relabel_nodes
except ImportError:
    from networkx.convert import relabel_nodes as nx_relabel_nodes

import models
import util

def random_graph(graph=None, size=10, cluster=True, layout=None):
    """Return a random graph (models.py Graph class).

    This was my original test graph function.  It makes a 2D
    non-periodic lattice, adds some random noise to it, and returns it.
    """
    if graph is None:
        g = grid_graph(dim=[size,size])
    g = nx_relabel_nodes(
        g, dict((n,i) for i,n in enumerate(g.nodes())))
    if layout:
        layout = nx.drawing.nx_pydot.pydot_layout(g)
#    layout = nx.drawing.layout.shell_layout(g)
    if cluster:
        clusterGraph(g, p=.05)

    # Set weights and make the Graph object.
    for a,b,d in g.edges(data=True):
        d['weight'] = -10

    G = models.Graph.fromNetworkX(g, layout=layout, defaultweight=1)
    return G

def relabel_nodes(g, copy=False):
    #for node in g.nodes():
    #    g[node]['name'] = node
    g = networkx.convert.relabel_nodes(
        g,
        dict((n,i) for i,n in enumerate(g.nodes())))
    return g

def clusterGraph(g, p=0.01):
    """Add random noise to a graph."""
    nodes = random.sample(g.nodes(), int(p*len(g)))
    for i in nodes:
        for n in g.neighbors(i):
            for n2 in g.neighbors(n):
                g.add_edge(i, n2)


def polopa_tribes(weightAllied=-1, weightHostile=2):
    """Highland Polopa Tribes, NetworkX representation.

    highPolopaTribes.xml comes from http://www.casos.cs.cmu.edu/computational_tools/datasets/external/gama/index2.html
    """

    fname = os.path.join(os.path.dirname(__file__),
                         "data/highlandPolopaTribes.xml")
    r = ElementTree.parse(fname).getroot()
    nodes = [ ]
    for node in r[0][0][0]:
        name = node.attrib['id']
        nodes.append(name)

    allied = [ ]
    hostile = [ ]
    def _parse(tree, lst):
        for link in tree:
            if link.tag == "properties":
                continue
            source = link.attrib['source']
            dest = link.attrib['target']
            lst.append((source, dest))
    _parse(r[0][1][0], lst=allied)
    _parse(r[0][1][1], lst=hostile)

    # Check it was parsed right
    assert len(nodes) == 16
    # confirm it is bi-directional
    for edges in (allied, hostile):
        for edge in edges:
            assert (edge[1], edge[0]) in edges

    g = networkx.Graph()
    for node in sorted(nodes,
                       key=lambda n: len([1 for x in allied if (n in x)])):
        g.add_node(node)
    for edge in allied:
        g.add_edge(*edge, weight=weightAllied, color='green')
    for edge in hostile:
        g.add_edge(*edge, weight=weightHostile, color='red')
    #print G.nodes()
    #print G.edges()
    #print G.edge['GAVEV']['OVE']
    #from fitz import interactnow
    return g

def polopa_tribes_peter(fname="highland_one.txt"):
    fname = os.path.join(os.path.dirname(__file__),
                         "data/highland_peter/"+fname)
    f = open(fname)
    f.readline() ; f.readline() ; f.readline()

    lines = f.readlines()
    lines = [ l.strip().split() for l in lines ]
    array = [ [int(weight) for weight in row] for row in lines ]

    array = numpy.asarray(array)
    array = -array
    g = util.networkx_from_matrix(array, ignore_values=(0, 1))
    return g

def dolphins(weightFriend=-1):
    """Dolphins test graph, NetworkX representation.

    http://www-personal.umich.edu/~mejn/netdata/
    http://www-personal.umich.edu/~mejn/netdata/dolphins.zip

    Dolphin social network: an undirected social network of frequent
    associations between 62 dolphins in a community living off
    Doubtful Sound, New Zealand. Please cite D. Lusseau, K. Schneider,
    O. J. Boisseau, P. Haase, E. Slooten, and S. M. Dawson, Behavioral
    Ecology and Sociobiology 54, 396-405 (2003). Thanks to David
    Lusseau for permission to post these data on this web site.
    """
    fname = os.path.join(os.path.dirname(__file__),
                         "data/dolphins.gml")
    g = networkx.read_gml(fname)
    g = networkx.Graph(g)
    for a,b in g.edges_iter():
        #d['weight'] = -1
        g.edge[a][b]['weight'] = weightFriend
    return g

def nussinov_256node(weight=-1):
    """Load the 256-node hierarchical graph from PRE 80, 016109 (2009)

    NetworkX representation.
    """
    fname = os.path.join(os.path.dirname(__file__),
                         "data/256node/example1.3.gml")
    g = networkx.read_gml(fname)
    g = networkx.Graph(g)
    for a,b in g.edges_iter():
        #d['weight'] = -1
        g.edge[a][b]['weight'] = weight
    return g


def hierarchical_graph(l, probs, random=random):
    """A hierarchical graph

    Example usage.  This call should approximately reproduce the graph
    in Figure 1a from PRE 80, 016109 (2009):
    hierarchical_graph(
      [[16, 17], [18, 18], [5, 19, 18], [16, 19, 22, 12], [11, 17, 16, 17, 15]]
      [.1, .3, .9]
    )
    """
    if isinstance(l, int):
        g = networkx.generators.random_graphs.erdos_renyi_graph(l, probs[0],
                                                  seed=random.randint(0, 1e9))
        return g

    subgraphs = [ hierarchical_graph(n, probs[1:], random=random) for n in l ]
    # relabel nodes
    for i, subg in enumerate(subgraphs):
        nodes = subg.nodes()
        mapping = dict((n, (i,n)) for n in nodes)
        g = networkx.convert.relabel_nodes(subg, mapping)
        subgraphs[i] = g # in-place
    g = subgraphs[0].__class__(name=str(l)+"  "+str(probs))
    for subg in subgraphs:
        print "sg:",
        print "  ", subg.number_of_nodes(), subg.number_of_edges()
        print "  ", sorted(subg.nodes())
        g.add_nodes_from(subg.nodes_iter(data=True))
        g.add_edges_from(subg.edges_iter(data=True))
    # Add links between the different subgraphs
    for i0, sg0 in enumerate(subgraphs):
        for i1, sg1 in enumerate(subgraphs[i0+1:]):
            for n0 in sg0.nodes_iter():
                for n1 in sg1.nodes_iter():
                    print "random adding edge:", n0, n1
                    if random.uniform(0, 1) < probs[0]:
                        print "-> yes"
                        g.add_edge(n0, n1)
                        pass
    return g

def _test_hierarchical_graph():
    # Return a test hierarchical_graph
    rnd = random.Random(20)
    g = hierarchical_graph([10, 5],
                           [.1, .5, .9],
                           random=rnd)
    print
    print "final:"
    print "  ", g.number_of_nodes(), g.number_of_edges()
    print "   nodes:", sorted(g.nodes())


if __name__ == "__main__":
    print polopa_tribes()
