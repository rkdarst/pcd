# Richard Darst, July 2011

import collections
import copy
import ctypes
import itertools
import math
from math import exp, log, floor, ceil
import numpy
import cPickle as pickle
import random
import sys

import networkx as nx


from cmodels import cGraph_fields
import cmodels
import util

log2 = lambda x: math.log(x, 2)

class _cobj(object):
    def _allocStruct(self, cModel):
        self._struct = cModel()
        self._struct_p = ctypes.pointer(self._struct)

    def _allocArray(self, name, array=None, **args):
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
        if array is None:
            array = numpy.zeros(**args)
        self.__dict__[name] = array
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
    def __init__(self, graph, layout=None, randomize=True):
        self._fillStruct(len(graph.nodes()))
        self.n     = len(graph.nodes())
        self.Ncmty = len(graph.nodes())

        self._graph  = graph
        self._layout = layout

        self.cmtyCreate(randomize=randomize)
    def _fillStruct(self, n=None):
        _cobj._allocStruct(self, cmodels.cGraph)

        self._allocArray("cmtyll", shape=[n]*2)
        self._allocArray("cmtyl", shape=n, dtype=ctypes.c_void_p)
        for i, row in enumerate(self.cmtyll):
#            print hex(row.ctypes.data)
#            print row.ctypes.data_as(self._struct.cmtyi._type_)
            self.cmtyl[i] = row.ctypes.data#_as(self._struct.cmtyi._type_)

        self._allocArray("cmty",  shape=n)
        self._allocArray("cmtyN", shape=n)
        self._allocArray("interactions", shape=(n, n))

    #def __getinitargs__(self):
    #    return self._graph, self._layout
    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_struct']
        del state['_struct_p']
        #for name, type_ in self._struct._fields_:
        #    if name in state  and isinstance(state[name], numpy.ndarray):
        #        setattr(self, name, state[name])
        #        del state[name]
        return state
    def __setstate__(self, state):
        self.__init__(state['_graph'], state['_layout'])
        #self._fillStruct(n=state['n'])
        for name, type_ in self._struct._fields_:
            #print name
            if name in state and isinstance(state[name], numpy.ndarray):
                self._allocArray(name, array=state[name])
                del state[name]
        self.__dict__.update(state)
    def copy(self):
        """Return a copy of the object (not sharing data arrays).
        """
        # Right now we turn this into a string and back, in order to
        # make actual deep copies of all the arrays.
        return pickle.loads(pickle.dumps(self))

    def cmtyCreate(self, cmtys=None, randomize=True):
        """Create initial communities for a model.

        cmtys=<array of length nNodes>: initialize for this initial
        array

        randomize=<bool>: if initial cmtys array is not given, should
        the innitial array be range(nNodes), or should it be randomized?

        By default, arrays are initialized to one particle per
        community, in random order.
        """
        if cmtys:
            pass
        else:
            cmtys = range(self.n)
            if randomize:
                random.shuffle(cmtys)
        self.cmty[:] = cmtys
        self.cmtyN[:] = 1
        self.cmtyListInit()

    def cmtyListInit(self):
        """Initialize community list.

        The community list (self.cmtyll[community_index, ...]) is a
        data structure optimized for the community lookups.  This
        function fills self.cmtyll from self.cmty[i].
        """
        print "cmtyListInit"
        cmodels.cmtyListInit(self._struct_p)
    def cmtyListCheck(self):
        print "cmtyListCheck"
        errors = cmodels.cmtyListCheck(self._struct_p)
        assert errors == 0, "We have %d errors"%errors

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
                node_size=150,
                labels = dict(zip(self._graph.nodes(), self.cmty)),
                pos=self._layout)
        #plt.margins(tight=True)
        plt.show()

    def energy(self, gamma):
        """Return energy due to just community c"""
        return cmodels.energy(self._struct_p, gamma)
    def energy_cmty(self, gamma, c):
        """Return energy due to just community c"""
        return cmodels.energy_cmty(self._struct_p, gamma, c)
    @property
    def q(self):
        "number of clusters"
        q = len([1 for n in self.cmtyN if n>0])
        assert q == len(set(self.cmty))
        return q
    @property
    def entropy(self):
        N = self.n
        H = sum( (n/float(N))*log2(n/float(N))  for n in self.cmtyN
                 if n!=0 )
        return - H

    def minimize(self, gamma):
        print "beginning minimization (n=%s, gamma=%s)"%(self.n, gamma)

        changesCombining = None
        for round_ in itertools.count():
            changes = cmodels.minimize(self._struct_p, gamma)
            #print self.cmty
            #print set(self.cmty),
            print "  (r%2s) cmtys, changes: %4d %4d"%(round_, self.q, changes)

            if changes == 0 and changesCombining == 0:
                break

            # If we got no changes, then try combining communities
            # before breaking out.
            if changes == 0:
                changesCombining = cmodels.combine_cmtys(self._struct_p, gamma)
                print "  (r%2s) cmtys, changes: %4d %4d"%(
                                           round_, self.q, changesCombining), \
                      "(<-- combining communities)"
            
            # If we have no changes in regular and combinations
            if changes == 0 and changesCombining == 0:
                break
            if round_ > 50:
                print "  Exceeding maximum number of rounds."
                break
        #print set(self.cmty),

    def remapCommunities(self):
        """Collapse all communities into the lowest number needed to
        represent everything."""
        print "remapping communities"
        def find_empty_cmty():
            """Start at zero, and find the first (lowest-index) empty
            community.
            """
            for m in range(self.Ncmty):
                if self.cmtyN[m] == 0:
                    assert len([ 1 for i in range(self.n)
                                 if self.cmty[i] == m  ]) == 0
                    return m

        for oldCmty in range(self.Ncmty-1, -1, -1):
            #print "oldCmty:", oldCmty
            # Ignore when we already have a community
            if self.cmtyN[oldCmty] == 0:
                continue
            # Where do we store this?
            newCmty = find_empty_cmty()
            #print "newCmty:", newCmty
            if newCmty is None:
                break
            assert newCmty != oldCmty
            if newCmty > oldCmty:
                #print "No more communities empty"
                break
            # We have established we want to move oldCmty to newCmty,
            # now do it.
            self.cmtyll[newCmty, :] = self.cmtyll[oldCmty, :]
            for i in range(self.n):
                if self.cmty[i] == oldCmty:
                    self.cmty[i] = newCmty
            self.cmtyN[newCmty] = self.cmtyN[oldCmty]
            self.cmtyN[oldCmty] = 0

        self.cmtyListCheck()

    def _test(self):
        return cmodels.test(self._struct_p)

class MultiResolutionCorrelation(object):
    def __init__(self, gamma, Gs):
        pairs = [ ]
        for i, G0 in enumerate(Gs):
            for j, G1 in enumerate(Gs[i+1:]):
                pairs.append((G0, G1))

        Is = [ util.mutual_information(G0,G1) for G0,G1 in pairs ]

        VI = numpy.mean([G0.entropy + G1.entropy - 2*mi
                         for ((G0,G1), mi) in zip(pairs, Is)])
        In = numpy.mean([2*mi / (G0.entropy + G1.entropy)
                         for ((G0,G1), mi) in zip(pairs, Is)])

        self.gamma = gamma
        self.q = min(G.q for G in Gs)
        self.E = min(G.energy(gamma) for G in Gs)
        self.entropy = numpy.mean([G.entropy for G in Gs])

        self.I = numpy.mean(Is)
        self.VI = VI
        self.In = In

class MultiResolution(object):
    def __init__(self, low, high):
        self.indexLow  = int(floor(util.logTimeIndex(low )))
        self.indexHigh = int(ceil (util.logTimeIndex(high)))

        self._data = { } #collections.defaultdict(dict)
    def do(self, Gs, trials=10):
        for index in range(self.indexLow, self.indexHigh+1):
            gamma = util.logTime(index)
            #minima = [ ]
            minGs = [ ]
            for G in Gs:
                thisGs = [ ]
                for i in range(trials):
                    G.cmtyCreate() # randomizes it
                    G.minimize(gamma)
                    #minima.append((G.energy(gamma), G.q))
                    thisGs.append(copy.copy(G))
                minG = min(thisGs, key=lambda G: G.energy)
                minGs.append(minG)

            self._data[index] = MultiResolutionCorrelation(gamma, Gs)
    def print_(self):
        import matplotlib.pyplot as pyplot
        MRCs = [ MRC for i, MRC in sorted(self._data.iteritems()) ]

        gammas    = [mrc.gamma     for mrc in MRCs]
        qs        = [mrc.q         for mrc in MRCs]
        Es        = [mrc.E         for mrc in MRCs]
        entropies = [mrc.entropy   for mrc in MRCs]

        Is        = [mrc.I         for mrc in MRCs]
        VIs       = [mrc.VI        for mrc in MRCs]
        Ins       = [mrc.In        for mrc in MRCs]


        pyplot.semilogx(gammas, qs)
        pyplot.xlabel('\gamma')
        pyplot.ylabel('q')
        #pyplot.show()
        pyplot.ion()
        from fitz import interactnow

def random_graph(graph=None, size=10, cluster=True, layout=None):
    print "making random graph"
    from networkx.generators.classic import grid_graph
    from networkx.convert import relabel_nodes
    if graph is None:
        g = grid_graph(dim=[size,size])
    g = relabel_nodes(g, dict((n,i) for i,n in enumerate(g.nodes()) ))
    if layout:
        layout = nx.drawing.nx_pydot.pydot_layout(g)
#    layout = nx.drawing.layout.shell_layout(g)

    if cluster:
        clusterGraph(g, p=.05)
    G = Graph(g, layout=layout)
    G.interactions[:] = 1
    for i in g.nodes():
        for j in g.neighbors(i):
            G.interactions[i,j] = -10
            G.interactions[j,i] = -10
    return G
        
def clusterGraph(g, p=0.01):
    print "re-clustering"
    nodes = random.sample(g.nodes(), int(p*len(g)))
    for i in nodes:
        for n in g.neighbors(i):
            for n2 in g.neighbors(n):
                g.add_edge(i, n2)



if __name__ == "__main__":
    command = None
    if len(sys.argv) > 1:
        command = sys.argv[1]

    if command == "test":
        import networkx.generators.classic
        G = random_graph(size=5)
        G.cmtyListCheck()
        #G.test()
        G.minimize(gamma=1.0)
        G.remapCommunities()
        G.remapCommunities()
        G.cmtyListCheck()
        G.viz()
        print G.q

        print G.cmty
        print G.cmtyN
        #import cPickle as pickle
        #s = pickle.dumps(G)
        #G2 = pickle.loads(s)
        #G2.viz()
    if command == "multi":
        Gs = [ random_graph(size=20) for _ in range(10) ]
        M = MultiResolution(low=.1, high=10)
        M.do(Gs=Gs, trials=1)
        M.print_()
