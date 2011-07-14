# Richard Darst, July 2011

import ctypes
import itertools
import numpy
import random

import networkx as nx
#import 

from cmodels import cGraph_fields
import cmodels

class _cobj(object):
    def __init__(self, cModel):
        self._struct = cModel()
        self._struct_p = ctypes.pointer(self._struct)

    def _allocArray(self, name, **args):
        """Allocate an array in a way usable by both python and C.
        
        `name` is the name of the array to allocate: it will be
        self.`name` and self.SD.`name`.  `**args` are arguments passed
        to `numpy.zeros` to allocate the arrays.
        
        If you want to initilize the array to something, you have to
        do it afterwards, like this:
        self.lattsite[:] = S12_EMPTYSITE
        """
        # This code can be hard to understand, since we can't use the
        # name, but it is equivalent to this for the name `lattsite`:
        ##  self.__dict__["lattsite"] = numpy.zeros(**args)
        ##  self.SD.lattsite = self.lattsite.ctypes.data
        # Get dtype if not already a argument
        if 'dtype' not in args:
            args['dtype'] = getattr(self._struct, name)._type_
        self.__dict__[name] = array = numpy.zeros(**args)
        setattr(self._struct, name,
                array.ctypes.data_as(ctypes.POINTER(
                    getattr(self._struct, name)._type_)))
        #setattr(self._struct, name, array.ctypes.data)
    def __getattr__(self, attrname):
        """Wrapper to proxy attribute gets to the C SimData struct
        """
        if attrname in cGraph_fields:
            return getattr(self._struct, attrname)
        raise AttributeError("No Such Attribute: %s"%attrname)
    def __setattr__(self, attrname, value):
        """Wrapper to proxy attribute sets to the C SimData struct
        """
        if attrname in cGraph_fields:
            setattr(self._struct, attrname, value)
        else:
            self.__dict__[attrname] = value


class Graph(_cobj, object):
    def __init__(self, graph, layout=None):
        _cobj.__init__(self, cmodels.cGraph)

        self._graph  = graph
        self._layout = layout

        self.n = len(graph.nodes())
        self._allocArray("cmty", shape=self.n)
        self._allocArray("cmtyN", shape=self.n)
        self.cmty[:] = range(self.n)
        self.cmtyN[:] = 1
        self._allocArray("interactions", shape=(self.n, self.n))

        self._allocArray("cmtyii", shape=[self.n]*2)
        self._allocArray("cmtyi", shape=self.n, dtype=ctypes.c_void_p)
        for i, row in enumerate(self.cmtyii):
            print hex(row.ctypes.data)
#            print row.ctypes.data_as(self._struct.cmtyi._type_)
            self.cmtyi[i] = row.ctypes.data#_as(self._struct.cmtyi._type_)


    def viz(self):
        print "vizing"
        import matplotlib.pyplot as plt
        #from fitz import interactnow

        colorMap = list(range(self.n))
        r = random.Random()
        r.seed('color mapping')
        r.shuffle(colorMap)

        nx.draw(self._graph,
                         node_color=[colorMap[x] for x in self.cmty],
                         labels = dict(zip(self._graph.nodes(), self.cmty)),
                         pos=self._layout)
        #plt.margins(tight=True)
        plt.show()

    def minimize(self, gamma):
        print "beginning minimization"

        for round_ in itertools.count():
            print "Minimizing, round", round_
            changes = cmodels.minimize(self._struct_p, gamma)
            print self.cmty
            print set(self.cmty), len(set(self.cmty))
            if changes == 0:
                break
            if round_ > 50:
                break

    def remapCommunities(self):
        """Collapse all communities into the lowest number needed to
        represent everything."""
        print "remapping communities"
        #oldCommunities = self.cmty[:]
        mapping = { }
        for i in range(len(self.cmty)):
            old = self.cmty[i]
            if old not in mapping:
                mapping[old] = len(mapping)
            self.cmty[i] = mapping[old]

    def test(self):
        return cmodels.test(self._struct_p)
        


def random_graph(graph=None, size=10):
    print "making random graph"
    from networkx.generators.classic import grid_graph
    from networkx.convert import relabel_nodes
    if graph is None:
        g = grid_graph(dim=[size,size])
    g = relabel_nodes(g, dict((n,i) for i,n in enumerate(g.nodes()) ))
    layout = nx.drawing.nx_pydot.pydot_layout(g)
#    layout = nx.drawing.layout.shell_layout(g)

    clusterGraph(g, p=.05)
    G = Graph(g, layout=layout)
    G.interactions[:] = 10
    for i in g.nodes():
        for j in g.neighbors(i):
            G.interactions[i,j] = -100
            G.interactions[j,i] = -100
    return G
        
def clusterGraph(g, p=0.01):
    print "re-clustering"
    nodes = random.sample(g.nodes(), int(p*len(g)))
    for i in nodes:
        for n in g.neighbors(i):
            for n2 in g.neighbors(n):
                g.add_edge(i, n2)



if __name__ == "__main__":
    import networkx.generators.classic
    G = random_graph(size=20)
    #G.test()
    G.minimize(gamma=1)
    G.remapCommunities()
    G.viz()
