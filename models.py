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
import time
import numbers

import networkx

import anneal
import cmodels
import util

log2 = lambda x: math.log(x, 2)

NO_CMTY = -1


class Graph(anneal._GraphAnneal, cmodels._cobj, object):
    """Core Potts-based community detection object.

    This returns an object which can be used to do a community
    detection on a graph.

    There is an important distinction between this object, and the
    underlying NetworkX graph object.  The NetworkX object is a
    Pythonic representation of the graph, but can not be used for
    minimization itself.  The Graph object is the thing which can be
    used for minimizing.  The Graph object can be created from a
    NetworkX object easily.  A typical use is to set up a NetworkX
    graph from loading a file, or your other simulation, and then
    using that to set up Graph object.

    """
    verbosity = 2
    _use_overlap = False

    def __init__(self, N, randomize=True, rmatrix=False, overlap=False,
                 sparse=False):
        """Initialize Graph object.

        For more easier minimization, look at the various class
        methods Graph.from* .

        N: number of nodes in this graph

        randomize: initialize communities randomly (default True).

        rmatrix: if True, initialize a rmatrix (matrix of repulsive
        terms).  Really, this just makes
          E = imatrix[n,m] + gamma*rmatrix[n,m]
        Without rmatrix,
          E = imatrix[n,m] if value<0 and gamma*imatrix[n,m] if value>0

        overlap: If true, will result in each .minimize() call
        attempting to find overlapping communities at the end of each
        minimization cycle.
        """
        self._fillStruct(N, sparse=sparse)
        self.N = N
        self.oneToOne = 1
        if sparse:
            self.hasFull = 0
        else:
            self.hasFull = 1
        self.hasPrimaryCmty = 1
        if rmatrix:
            self._allocRmatrix()
        if overlap:
            self.setOverlap(True)
        #seenList = cmodels.LList(N)
        #self.seenList = ctypes.pointer(seenList)
        #self.__dict__['seenList'] = seenList
        self._struct.seenList = cmodels.C.SetInit()
        cmodels.C.hashCreate(self._struct_p)
        self.cmtyCreate(randomize=randomize)
    def __del__(self):
        cmodels.C.hashDestroy(self._struct_p)
        cmodels.C.SetDestroy(self._struct.seenList)
    def _fillStruct(self, N=None, sparse=False):
        """Fill C structure."""
        cmodels._cobj._allocStruct(self, cmodels.cGraph)

        #self._allocArray("cmtyll", shape=[N]*2)
        #self._allocArray("cmtyl", shape=N, dtype=ctypes.c_void_p)
        #self._allocArrayPointers(self.cmtyl, self.cmtyll)
        self._allocArray("cmty",  shape=N)
        self._allocArray("cmtyN", shape=N)
        self._allocArray("tmp",  shape=N)
        self._allocArray("randomOrder", shape=N)
        self._allocArray("randomOrder2", shape=N)
        self.randomOrder[:]  = numpy.arange(N)
        self.randomOrder2[:] = numpy.arange(N)
        if not sparse:
            self._allocArray("imatrix", shape=(N, N))
        self._allocArray('cmtyListHash', shape=N, dtype=cmodels.c_void_p,
                        dtype2=cmodels.c_void_p)
    def _allocRmatrix(self):
        """Allocate a matrix for repulsive interactions"""
        self._allocArray("rmatrix", shape=(self.N, self.N))

    def __getstate__(self):
        """Get state for pickling"""
        state = self.__dict__.copy()
        del state['_struct']
        del state['_struct_p']
        del state['cmtyListHash']
        state['__extra_attrs'] = { }
        for name in ('N', 'Ncmty', 'oneToOne', 'hasSparse', 'hasFull',
                     'hasPrimaryCmty',
                     'simatrixDefault', 'srmatrixDefault', 'simatrixLen'):
            state['__extra_attrs'][name] = getattr(self, name)
        #for name, type_ in self._struct._fields_:
        #    if name in state  and isinstance(state[name], numpy.ndarray):
        #        setattr(self, name, state[name])
        #        del state[name]
        state['cmtystate'] = self.getcmtystate()
        return state
    def __setstate__(self, state):
        """Set state when unpickling"""
        sparse = False
        if not state['__extra_attrs'].get('hasFull', 1):
            sparse = True
        self.__init__(state['__extra_attrs']['N'], sparse=sparse)
        #self._fillStruct(n=state['n'])
        for name, type_ in self._struct._fields_:
            #print name
            if name in state and isinstance(state[name], numpy.ndarray):
                self._allocArray(name, array=state[name])
                del state[name]
        #self._allocArrayPointers(self.cmtyl, self.cmtyll)
        for k, v in state['__extra_attrs'].items():
            setattr(self, k, v)
        del state['__extra_attrs']
        cmtystate = state.pop('cmtystate')
        self.__dict__.update(state)
        self._allocArray('cmtyListHash', shape=self.N, dtype=cmodels.c_void_p,
                        dtype2=cmodels.c_void_p)
        cmodels.C.hashCreate(self._struct_p)
        cmodels.C.hashInit(self._struct_p)
        self._struct.seenList = cmodels.C.SetInit()
        self.setcmtystate(cmtystate)
    def copy(self):
        """Return a copy of the object (not sharing data arrays).
        """
        # Right now we turn this into a string and back, in order to
        # make actual deep copies of all the arrays.  FIXME.
        #new = pickle.loads(pickle.dumps(self, -1))
        # Better way, sharing some data structures.  This means we
        # just have a shallow copy, so if you change mutable objects
        # it will change on all objects!
        new = self.__new__(self.__class__)
        state = self.__getstate__()
        for name, type_ in self._struct._fields_:
            if name in state and isinstance(state[name], numpy.ndarray) \
                   and name not in ('imatrix', 'rmatrix',
                                    'simatrix', 'srmatrix',
                                    'simatrixN', 'simatrixId', 'simatrixIdl'):
                state[name] = state[name].copy()
        new.__setstate__(state)


        return new
    def getcmtystate(self):
        """Return an object representing the communities.

        This can be pickeled and re-loaded into this object later.

        See also .cmtyDict().
        """
        state = { }
        state['version'] = 0
        state['oneToOne'] = self.oneToOne
        state['hasPrimaryCmty'] = self.hasPrimaryCmty
        state['cmtyContents'] = { }
        for c in self.cmtys():
            state['cmtyContents'][c] = tuple(self.cmtyContents(c))
        state['cmtyList'] = tuple(self.cmty)
        return state
    def setcmtystate(self, state):
        """Re-load state of communities from state dict.
        """
        self.cmtyListClear()
        if state['version'] == 0:
            seenNodes = set()
            actuallyOneToOne = True
            for c, nodes in sorted(state['cmtyContents'].items()):
                for n in nodes:
                    #print "clAO:", c, n, self.cmtyN[c]
                    self.cmtyListAddOverlap(c, n)
                # Detect if we have a repeat of any nodes.  If so, we
                # can not be one to one.
                if actuallyOneToOne and (seenNodes & set(nodes)):
                    actuallyOneToOne = False
                seenNodes |= set(nodes)
            # If stored state shows we should be oneToOne, and
            # detected state shows we are not, there is a problem.
            if state['oneToOne'] and not actuallyOneToOne:
                raise Exception("oneToOne value does not match.")
            self.oneToOne = state['oneToOne']
            self.hasPrimaryCmty = state['hasPrimaryCmty']
            self.cmty[:] = state['cmtyList']
        else:
            raise Exception("Unknown version of state to load communities.")
        #self.hashInit() #cmtyListAddOverlap keeps hashes properly working now.
        #self.check() # can't check non-oneToOne (yet)
    def hash(self):
        """Hash of state of self.

        This is not a real Python hash, since this is a mutable object."""
        return hash((
            self.oneToOne,
            self.hasPrimaryCmty,
            self.hasFull, self.hasSparse,
            tuple(self.cmty.flat),
            tuple(self.cmtyN.flat),
            tuple(self.cmtys()), # possibly redundant
            tuple(tuple(sorted(self.cmtyContents(c))) for c in self.cmtys()),
            ))
    def has_rmatrix(self):
        """True if system has self.rmatrix defined
        """
        if not isinstance(self.rmatrix, numpy.ndarray) and not self.rmatrix:
            return False
        return True
    def setOverlap(self, overlap):
        self._use_overlap = bool(overlap)

    #
    # Constructor methods
    #
    @classmethod
    def fromNetworkX(cls, graph, defaultweight=1,
                     coords=None, randomize=True,
                     diagonalweight=None):
        """Create a Graph glass from a NetworkX graph object.

        The NetworkX graph object is expected to have all edges have a
        attribute `weight` which is set to the interaction between the
        nodes.  Note that negative weight is attractive, positive
        weight is repulsive.

        `defaultweight` is the weight for all not explicitely defined
        pairs of nodes.  Note that positive weight is repulsive and
        negative weight is attractive.  `diagonalweight` is the weight
        of all diagonals on the matirx.

        `layout` is a mapping nodeid->(x,y) which is used in the
        .viz() method.

        `randomize` is passed to the constructor, randomizes node
        initial assignments.
        """
        G = cls(N=len(graph), randomize=randomize)
        G._graph = graph
        G.coords = coords

        G._nodeIndex = { }
        G._nodeLabel = { }

        # Set up node indexes (since NetworkX graph object nodes are
        # not always going to be integeras in range(0, self.N)).
        for i, name in enumerate(sorted(graph.nodes())):
            G._nodeIndex[name] = i
            G._nodeLabel[i] = name
            graph.node[name]['index'] = i

        # Default weighting
        G.imatrix[:] = defaultweight
        # Default diagonal weighting
        if diagonalweight is not None:
            for i in range(G.N):
                G.imatrix[i, i] = diagonalweight
        # All explicit weighting from the graph
        for n0 in graph.nodes():
            for n1 in graph.neighbors(n0):
                i0 = graph.node[n0]['index']
                i1 = graph.node[n1]['index']
                G.imatrix[i0,i1] = graph[n0][n1]['weight']
        # Check that the matrix is symmetric (FIXME someday: directed graphs)
        if not numpy.all(G.imatrix == G.imatrix.T):
            print "Note: interactions matrix is not symmetric"
        return G

    @classmethod
    def from_imatrix(cls, imatrix, coords=None, rmatrix=None):
        """Create a graph structure from a matrix of interactions.

        imatrix: the interaction matrix to use.

        coords: optional mapping for nodes->coordinates.  Can be
        either an array or mapping.

        rmatrix: If given, initialize an rmatrix with this.  See
        __init__ documentation for what this means.
        """
        # This should be improved sometime.
        assert len(imatrix.shape) == 2
        assert imatrix.shape[0] == imatrix.shape[1]

        G = cls(N=imatrix.shape[0])
        G.imatrix[:] = imatrix

        if rmatrix is not None:
            assert len(rmatrix.shape) == 2
            assert rmatrix.shape[0] == rmatrix.shape[1]
            assert imatrix.shape[0] == rmatrix.shape[0]
            G._allocRmatrix()
            G.rmatrix[:] = rmatrix

        G.coords = coords
        return G
    @classmethod
    def from_coords_and_efunc(cls, coords, efunc, boxsize=None,
                              refunc=None):
        """Create graph structure from coordinates+energy function.

        This helper class method is used to create an imatrix from a
        list of coordinates ( [[x0,y0,...], [x1,y1,...], ...]) and a
        efunc, which maps distance -> energy_of_interaction at that
        distance.  The efunc must be able to be applied to an array.

        If `boxsize` is given, it is used as a periodic boundary
        length: either a number or a [Lx, Ly, ...] array.

        If `refunc` is given, this is an efunc to be used to create
        the repulsive matrix.  See __init__ documentation for what
        this means.
        """
        G = cls(N=coords.shape[0])
        G.coords = coords
        G._makematrix_fromCoordsAndEfunc(coords, efunc, boxsize=boxsize)
        if refunc is not None:
            G._allocRmatrix()
            G._makematrix_fromCoordsAndEfunc(coords, refunc, boxsize=boxsize,
                                            matrix="rmatrix")
        return G

    def _makematrix_fromCoordsAndEfunc(self, coords, efunc, boxsize=None,
                                      matrix="imatrix"):
        orig_settings = numpy.seterr()
        numpy.seterr(all="ignore")
        matrix = getattr(self, matrix)
        for n1 in range(len(coords)):
            delta = coords[n1] - coords
            if not boxsize is None:
                delta -=  numpy.round(delta/boxsize)*boxsize
            delta = delta**2
            dist = numpy.sum(delta, axis=1)
            dist = numpy.sqrt(dist)
            matrix[n1, :] = efunc(dist)
        numpy.seterr(**orig_settings)
    @classmethod
    def from_mapping(cls, nodes, interactions,
                     maxinteractions=None,
                     defaultweight=1):
        """

        interactions: Callable interactions[node1][node2]
        """
        G = cls(N=len(nodes))
        G.imatrix[:] = defaultweight

        G._nodeIndex = { }
        G._nodeLabel = { }
        for i, label in enumerate(nodes):
            G._nodeIndex[label] = i
            G._nodeLabel[i] = label

        if maxinteractions is None:
            maxinteractions = max([len(connections) for connections
                                   in interactions.values()])

        for n1, connections in interactions.items():
            for n2, interaction in connections.items():
                #print interaction
                G.imatrix[G._nodeIndex[n1] , G._nodeIndex[n2]] = interaction

        return G
    @classmethod
    def from_sparseiter(cls, nodes, weights, default, maxconn,
                        imatrixDefault=0):

        nodeIndex = { }
        nodeLabel = { }
        for i, label in enumerate(nodes):
            nodeIndex[label] = i
            nodeLabel[i] = label
        #if nnodes is None:
        #    maxconn = len(nodeIndex)
        nnodes = len(nodeIndex)

        G = cls(N=nnodes, sparse=True)
        G._nodeIndex = nodeIndex
        G._nodeLabel = nodeLabel

        G.simatrixDefault = imatrixDefault
        G.srmatrixDefault = default
        G._alloc_sparse(simatrixLen=maxconn)
        G.simatrix[:] = numpy.nan

        simatrixN = G.simatrixN
        simatrix = G.simatrix
        simatrixId = G.simatrixId
        for node1, node2, weight in weights:
            if node1 == node2:
                continue
            i = nodeIndex[node1]
            j = nodeIndex[node2]
            idx = simatrixN[i]
            assert idx < maxconn
            simatrix[i,idx] = weight
            simatrixId[i,idx] = j
            simatrixN[i] += 1
        return G


    def _alloc_sparse(self, simatrixLen):
        self.hasSparse = 1
        self.simatrixLen = simatrixLen
        self._allocArray('simatrix', shape=(self.N, simatrixLen))
        if isinstance(self.rmatrix, numpy.ndarray):
            self._allocArray('srmatrix', shape=(self.N, simatrixLen))
        self._allocArray('simatrixN', shape=self.N)
        self._allocArray('simatrixId', shape=(self.N, simatrixLen))
        self._allocArray("simatrixIdl", shape=self.N, dtype=ctypes.c_void_p)
        self._allocArrayPointers(self.simatrixIdl, self.simatrixId)

    def make_sparse(self, default, cutoff=None,
                    imatrixDefault=0, imatrixCutoff=None):
        if self.hasSparse:
            return
        assert self.hasFull # have full to create sparse from it.
        assert not isinstance(self.rmatrix, numpy.ndarray), \
                                       "sparse does not support rmatrix yet"
        if default == 'auto':
            default = self.imatrix.max()
        if cutoff is None:
            cutoff = default

        self.simatrixDefault = imatrixDefault
        self.srmatrixDefault = default

        imatrix = self.imatrix
        #simatrixLen = 0
        #for row in imatrix:
        #    x = numpy.sum(row < minval)
        #    simatrixLen = max(x, simatrixLen)
        simatrixLen = (imatrix<cutoff).sum(axis=1).max()
        print "Making graph sparse: %d nodes, %d reduced nodes"%(
            self.N, simatrixLen)
        print "  imatrix max/min:", numpy.min(imatrix), numpy.max(imatrix)
        self._alloc_sparse(simatrixLen=simatrixLen)

        for i,j in zip(*numpy.where(imatrix < cutoff)):
            if i == j:
                continue
            idx = self.simatrixN[i]
            #assert idx < self.simatrixLen
            self.simatrix[i,idx] = imatrix[i, j]
            self.simatrixId[i,idx] = j
            self.simatrixN[i] += 1
    def _check_sparse(self, imatrixDefault):
        #if numpy.max(self.seenList.data[:self.seenList.maxcount]) > 10000:
        #    from fitz import interactnow
        #    assert 0
        #
        for i, row in enumerate(self.simatrix):
            for _j, val in enumerate(row[0:self.simatrixN[i]]):
                j = self.simatrixId[i,_j]
                if not val == self.imatrix[i, j]:
                    from fitz import interactnow
                assert val == self.imatrix[i, j], "%d %d %d %d"%(i,j,_j,val)
        # Compare everything in imatrix, make sure it is in simatrix
        # if needed.
        for i, row in enumerate(self.imatrix):
            for j, val in enumerate(row):
                # We don't store on-diagonal values
                if i == j:
                    continue
                # If it isn't the default, it should be in in the matrices
                if val != imatrixDefault:
                    assert j     in self.simatrixId[i,:self.simatrixN[i]], \
                           "%d %d %d"%(i,j,val)
                # And if it is the default, it shouldn't be in.
                else:
                    assert j not in self.simatrixId[i,:self.simatrixN[i]], \
                           "%d %d %d"%(i,j,val)
        # Is it sorted?
        for i in range(self.N):
            assert list(sorted(self.simatrixId[i,:self.simatrixN[i]])) == \
                   list(self.simatrixId[i,:self.simatrixN[i]])

    #
    # Conversion methods
    #
    def make_networkx(self, ignore_values=()):
        """Return a networkX representation of the interaction matrix.

        This could be useful where you have created a graph from an
        imatrix or some other constructor, and you want a networkx
        graph out of it in order to do various graph
        algorithms/transforms on it."""
        g = util.networkx_from_matrix(self.imatrix,
                                      ignore_values=ignore_values)
        return g
    def supernode_networkx(self):
        """Make a networkx graph representing attrraction between cmtys.

        This returns a graph with nodes for each community and edges
        for every pair of communities with at least one attractive
        interaction.
        """
        g = networkx.Graph()
        g.add_nodes_from(self.cmtys())
        for c in self.cmtys():
            for n in self.cmtyContents(c):
                for m, interaction in enumerate(self.imatrix[n]):
                    if interaction < 0:
                        g.add_edge(c, self.cmty[m])
        return g
    def supernodeGraph(self, gamma, multiplier=100):
        """Great a sub-graph of super-nodes composed of communities.

        Each community is turned into a super-node, and energy of
        interaction between communities is calculated.

        WARNING: When self.rmatrix exists, this method returns a graph
        showing that every community is connected to every other
        community.
        """
        self.remap()
        N = self.q
        G = self.__class__(N=N)
        for c1 in range(N):
            for c2 in range(c1+1, N):
                # Find the energy of interaction between c1 and c2
                E = self.energy_cmty_cmty(gamma, c1, c2)
                G.imatrix[c1, c2] = E*multiplier
                G.imatrix[c2, c1] = E*multiplier
        return G
    def loadFromSupernodeGraph(self, G):
        """Reload our community assignments from a supernode Graph G.
        """
        # For each community in the sub-graph
        mapping = { }
        for c in range(G.N):
            for n in self.cmtyContents(c):
                mapping[n] = G.cmty[c]
        #print mapping
        print "Subgraph changes in q: %d -> %d"%(self.q, G.q)
        assert max(mapping) == len(mapping)-1
        #print self.cmty
        for n in range(self.N):
            #print "x", n, self.cmty[n], "->", mapping[n]
            self.cmty[n] = mapping[n]
        #print self.cmty
        self.cmtyListInit()
        #self.check()



    #
    # Utility methods
    #
    def check(self):
        """Do an internal consistency check."""
        print "cmtyListCheck"
        errors = cmodels.cmtyListCheck(self._struct_p)
        assert errors == 0, "We have %d errors"%errors
        #errors = self._check_sparse(imatrixDefault=500)
        #assert errors == 0, "We have %d errors (2)"%errors
        if self.oneToOne:
            assert self.q == len(set(self.cmty))
        # Verify that oneToOne is correct
        # (no node is in more than one community):
        if self.oneToOne:
            for c in self.cmtys():
                for n in self.cmtyContents(c):
                    assert self.cmty[n] == c
                    if self.cmty[n] != c:
                        errors += 1
        return errors
    def _test(self):
        return cmodels.test(self._struct_p)
    def _gen_random_order(self, randomize=True):
        """Generate random orders to use for iterations."""
        if not randomize:
            # Fake mode: actually sort them
            numpy.sort(self.randomOrder)
            numpy.sort(self.randomOrder2)
        else:
            numpy.random.shuffle(self.randomOrder)
            numpy.random.shuffle(self.randomOrder2)
    def _check_overlaps(self):
        """Check some theoretical properties of overlaps.

        1) once you apply overlaps, none of the communities are
        identical.
        """
        if self.oneToOne:
            print ("Warning: self.oneToOne is true, do you want to call "
                   "this method?")
        # Make a set of all community contents:
        cmtys = set()
        for c in self.cmtys():
            nodes = frozenset(self.cmtyContents(c))
            if nodes in cmtys:
                # We have a overlap - we hope that any
                print "Overlapping community:", c, nodes
        # Now make sure that there are no strict subsets.
        for s1 in cmtys:
            for s2 in cmtys:
                if s1 < s2:
                    print "Set1 is strict subset of set2 "
                    print s1
                    print s2
    def _interaction_summary(self):
        pass
    def is_lowest_gamma(self):
        """Return True if this graph reduced fully.

        This method tests that every community has *only* repulsive
        interactions with every other community."""
        # rmatrix can always be reduced further
        if self.rmatrix is not None:
            return False
        if self.q == 1:
            return True


    #
    # Methods that deal with communities and community lists
    #
    def cmtyCreate(self, cmtys=None, randomize=None):
        """Create initial communities for a model.

        Used to make initial random community assignments.  By
        default, arrays are initialized to one particle per community,
        in random order.

        `cmtys`=<array of length nNodes>: initialize with this initial
        array.  Default: None

        `randomize`=<bool>: Randomize the initial community
        assignments.  Default: True for cmtys=None, False otherwise.
        """
        if cmtys is None:
            cmtys = numpy.arange(self.N)
            randomize = True
        else:
            randomize = False
        if randomize:
            numpy.random.shuffle(cmtys)
        self.cmty[:] = cmtys
        #self.Ncmty = numpy.max(cmtys)+1 done in cmtyListInit()
        #self.cmtyN[:] = 1   # done in cmtyListInit()
        #self.oneToOne = 1   # done in C
        self.cmtyListInit()
    def cmtyListInit(self):
        """Initialize the efficient community lists.

        The community list (self.cmtyll[community_index, ...]) is a
        data structure optimized for the community lookups.  This
        function fills self.cmtyll from self.cmty[i].
        """
        if self.verbosity >= 3:
            print "cmtyListInit"
        cmodels.cmtyListInit(self._struct_p)
    def cmtyListClear(self):
        """Clear all particles from all communities"""
        #raise NotImplementedError("Update this for hash lists.")
        #cmodels.hashInit(self._struct_p)
        self.cmty[:] = NO_CMTY
        self.cmtyN[:] = 0
        self.Ncmty = 0
        self.hashClear()
    def cmtyListAdd(self, c, n):
        """Add node n to community c."""
        cmodels.cmtyListAdd(self._struct_p, c, n)
    def cmtyListAddOverlap(self, c, n):
        """Add node n to community c (overlaps allowed)."""
        cmodels.cmtyListAddOverlap(self._struct_p, c, n)
    def cmtyListRemove(self, c, n):
        """Remove node n from community c."""
        cmodels.cmtyListRemove(self._struct_p, c, n)
    def cmtyListRemoveOverlap(self, c, n):
        """Remove node n from community c (overlaps allowed)."""
        cmodels.cmtyListRemoveOverlap(self._struct_p, c, n)
    def cmtyMove(self, n, cold, cnew):
        """Move node n from cmty cold to cmty cnew.

        This assumes that particle is currently in cold, and is *not*
        currently in cnew."""
        return cmodels.cmtyMove(self._struct_p, n, cold, cnew)
    def cmtyMoveSafe(self, n, cold, cnew):
        """Move node n from cmty cold to cmty cnew.

        Does additional errorchecking to get around limitations of cmtyMove."""
        return cmodels.cmtyMoveSafe(self._struct_p, n, cold, cnew)
    def cmtySet(self, n, c):
        """Set the community of a particle."""
        raise NotImplementedError("cmtySet not implemented yet.")
    def cmtyContents(self, c):
        """Array of all nodes in community c."""
        cN = cmodels.c_int()
        contents = numpy.zeros(self.cmtyN[c], dtype=ctypes.c_int)
        cmodels.cmtyGetContents(self._struct_p, c,
                                contents.ctypes.data_as(cmodels.c_int_p),
                                ctypes.pointer(cN))
        return contents
    def cmtyContains(self, c, n):
        """True if cmty c contains node n."""
        return cmodels.isInCmty(self._struct_p, c, n)
    def cmtyIntersect(self, c1, c2):
        """Number of nodes in intersect of c1 and c2."""
        return cmodels.cmtyIntersect(self._struct_p, c1, self._struct_p, c2)
    def cmtyUnion(self, c1, c2):
        """Number of nodes in union of c1 and c2."""
        return cmodels.cmtyUnion(self._struct_p, c1, self._struct_p, c2)
    def cmtyIsSubset(self, csmall, cbig):
        """Is cmty csmall a subset of cmty cbig?"""
        return cmodels.cmtyIsSubset(self._struct_p, csmall, cbig)
    def cmtyDict(self):
        """Return dict of community state.

        Keys: community indexes
        Values: tuples of nodes in that community

        See also .getcmtystate()
        """
        return dict((c, tuple(self.cmtyContents(c))) for c in self.cmtys())
    def cmtys(self):
        """Return an iterator over all non-empty community IDs.

        This is a useful method as it only returns communities that
        are non-empty, unlike range(G.cmtyN)
        """
        return ( c for c in range(self.Ncmty) if self.cmtyN[c]>0 )
    def cmtyConnectivitySingle(self, c, nodes=None):
        if not nodes:
            nodes = tuple(self.cmtyContents(c))
        nodes = set(nodes)
        connections = [ ]
        for n in nodes:
            conn = 0
            connInCmty = 0
            for i in range(self.simatrixN[n]):
                n2 = self.simatrixId[n,i]
                interaction = self.simatrix[n,i]
                if interaction >= 0: continue
                conn += 1
                if n2 in nodes:
                    connInCmty += 1
            connections.append((conn, connInCmty))
        #print connections
        return \
           sum(i for c,i in connections)/float(sum(c for c,i in connections)),\
           sum(1 for c,i in connections if i>.5*c)/float(len(connections))
    def cmtyConnectivity(self, c=None, nodes=None):
        if nodes is not None:
            return numpy.mean([self.cmtyConnectivitySingle(c=None, nodes=ns)
                               for ns in nodes],
                              axis=0)
        if c is None:
            c = self.cmtys()
        return numpy.mean([self.cmtyConnectivitySingle(c=c_) for c_ in c],
                          axis=0)
    def hashInit(self):
        return cmodels.hashInit(self._struct_p)
    def hashClear(self):
        return cmodels.hashClear(self._struct_p)


    #
    # Methods that deal with energies or other thermodynamic measures
    # or statistics of the current state of the graph.
    #
    def energy(self, gamma):
        """Return total energy of the graph."""
        return cmodels.energy(self._struct_p, gamma)
    def energy_sparse(self, gamma):
        """Return total energy of the graph."""
        return cmodels.energy_sparse(self._struct_p, gamma)
    def energy_cmty(self, gamma, c):
        """Return energy due to community c."""
        return cmodels.energy_cmty(self._struct_p, gamma, c)
    def energy_cmty_n(self, gamma, c, n):
        """Return energy due to particle n if it was in community c."""
        return cmodels.energy_cmty_n(self._struct_p, gamma, c, n)
    def energy_cmty_cmty(self, gamma, c1, c2):
        """Total entery of interaction between two communities.

        This can be used for making the supernode graph, but is not
        used in normal energy calculations as this is not included in
        normal energy calculations."""
        return cmodels.energy_cmty_cmty(self._struct_p, gamma, c1, c2)
    def energy_n(self, gamma, n):
        """Energy of particle n in its own community."""
        return cmodels.energy_n(self._struct_p, gamma, n)
    @property
    def q_python(self):
        """Number of communities in the graph.

        The attribute Ncmty is an upper bound for the number of
        communites, but some of those communities may be empty.  This
        attribute returns the number of non-empty communities."""
        q = len([1 for c in range(self.Ncmty) if self.cmtyN[c] > 0])
        #q = sum(self.cmtyN > 0)
        #q = len(set(self.cmty))  #FIXME this is fast, but a bit wrong.
        #assert q == len(set(self.cmty))
        return q
    @property
    def q_c(self):
        """Number of communities in the graph.

        The attribute Ncmty is an upper bound for the number of
        communites, but some of those communities may be empty.  This
        attribute returns the number of non-empty communities.

        Optimized implementation in C."""
        return cmodels.q(self._struct_p)
    q = q_c
    @property
    def entropy_python(self):
        """Return the entropy of this graph.
        """
        N = float(self.N)
        H = sum( (n/N)*log2(n/N)  for n in self.cmtyN
                 if n!=0 )
        return -H
    @property
    def entropy_c(self):
        """Return the entropy of this graph.  Optimized impl. in C.
        """
        return cmodels.entropy(self._struct_p)
    entropy = entropy_c
    def n_counts(self):
        """Return countsogram of community sizes."""
        #community sizes can go from 1 - N
        n_counts = {}
        for c in self.cmtys():   # iterate over non-empty cmtys only
            size = self.cmtyN[c]
            n_counts[size] = n_counts.get(size, 0) + 1
        return n_counts
    def n_mean(self):
        """Average num nodes per community"""
        return sum(self.cmtyN[c] for c in self.cmtys())/float(self.q)



    #
    # Methods that deal with minimization/optimization
    #
    def trials(self, gamma, trials, initial='random',
               minimizer='greedy', **kwargs):
        """Minimize system using .minimize() and `trials` trials.

        This will minimize `trials` number of times, and set the final
        graph to be the state that had the best energy of all the
        trials.

        trials: Number of trials to minimize.

        minimize: Minimizer function to use.  Common ones will be
        self.minimize, self.overlapMinimize, self.anneal.  If a
        string, use getattr(self, name).  Default: 'minimize'.

        initial: Initial state for each trial.  Default is the string
        'random', which calls self.cmtyCreate() each trial to
        randomize.  The string 'current' will get current state and
        reload it each time.  Otherwise, it should be something from
        self.getcmtystate(), and will be set with self.setcmtystate()
        before each minimization round.

        **kwargs: Keyword arguments passed to minimizer function.
        Minimizer is called as minimize(gamma, **kwargs).
        """
        if isinstance(minimizer, str):
            minimizer = getattr(self, minimizer)
        if self.verbosity > 0:
            print "Trials: %d %s runs (gamma=%f)"%(trials,
                                                   minimizer.func_name,gamma)
        minE = float('inf')
        minCmtyState = self.getcmtystate()
        if initial == 'current':
            initial = self.getcmtystate()
        nChanges = [ ]

        for i in range(trials):
            if initial == 'random':
                self.cmtyCreate() # randomizes it
            else:
                self.setcmtystate(initial)
            changes = minimizer(gamma, **kwargs)
            nChanges.append(changes)
            thisE = self.energy(gamma)
            if thisE < minE:
                minE = thisE
                minCmtyState = self.getcmtystate()
            if self.verbosity >= 1.5:
                print "Trial %d, minimum energy %f"%(i, thisE)

        self.setcmtystate(minCmtyState)
        self.nChanges = numpy.mean(nChanges, axis=0)
        if self.verbosity > 0:
            print "Trials done: %d %s runs (gamma=%f, minE=%f)"%(
                                     trials, minimizer.func_name, gamma, minE),\
                                     self.nChanges
        return self.nChanges
    # Keep backwards compatibility for now.
    minimize_trials = trials

    def alternate(self, funcs, mode="restart"):
        """Alternately call funcA and funcB.
        """
        if self.verbosity >= 0:
            print "beginning alternating", " ".join(f.func_name for f in funcs)
        lastChanges = [ None ] * len(funcs)
        rounds = [ 0 ] * len(funcs)

        for round_ in itertools.count():
            # For each of the `funcs` passed,
            for i, func in enumerate(funcs):
                changes = func()
                self.remap(check=False)
                if changes > 0:
                    rounds[i] += 1
                lastChanges[i] = changes
                if self.verbosity >= 2:
                    print "  (r%2s) %s: cmtys, changes: %4d %4d"%(
                        round_, func.func_name, self.q, changes)
                # At the first function that makes any changes, start over
                if mode == "restart" and changes > 0:
                    break
            if None not in lastChanges and sum(lastChanges) == 0:
                break

            if round_ > 250:
                print "  Exceeding maximum number of rounds."
                break
        return tuple(rounds) + (sum(lastChanges), )


    def greedy(self, gamma):
        """Minimize the communities at a certain gamma.

        This function requires a good starting configuration.  To do
        that, use self.cmtyCreate() first.
        """
        if self.verbosity >= 0:
            print "beginning minimization (n=%s, gamma=%s)"%(self.N, gamma)
        changes = 0
        changesCombining = None
        roundsMoving = 0
        roundsCombining = 0

        for round_ in itertools.count():
            self._gen_random_order()
            changesMoving = self._greedy(gamma=gamma)
            if changesMoving > 0:
                roundsMoving += 1
            #if round_%2 == 0:
            #    self.remapCommunities(check=False)
            changes += changesMoving
            if self.verbosity >= 2:
                print "  (r%2s) greedy: cmtys, changes: %4d %4d"%(round_,self.q,
                                                          changesMoving)

            if changesMoving == 0 and changesCombining == 0:
                break

            # If we got no changes, then try combining communities
            # before breaking out.
            if changesMoving == 0:
                self.remap(check=False)
                changesCombining = self.combine(gamma=gamma)
                #changesCombining = self.combine_cmtys_supernodes(gamma=gamma)
                if changesCombining > 0:
                    roundsCombining += 1
                changes += changesCombining
                if self.verbosity >= 2:
                    print "  (r%2s) greedy: cmtys, changes: %4d %4d"%(
                                          round_, self.q, changesCombining), \
                                                  "(<-- combining communities)"

                self.remap(check=False)
                roundsCombining += 1
            # If we have no changes in regular and combinations, escape
            if changesMoving == 0 and changesCombining == 0:
                break
            if round_ > 250:
                print "  Exceeding maximum number of rounds."
                break
        #print set(self.cmty),
        if self._use_overlap:
            self.greedyOverlap(gamma)
        return roundsMoving, roundsCombining, changes
    minimize = greedy
    def _greedy(self, gamma):
        return cmodels.greedy(self._struct_p, gamma)
    def ovGreedy(self, gamma):
        """Attempting to add particles to overlapping communities.
        """
        if not self.hasSparse:
            self.make_sparse(default='auto')
        if self.verbosity > 0:
            print "Overlap minimize (gamma=%f)"%gamma
        changes = 0
        changes_add = 0
        changes_remove = None
        roundsAdding = 0
        roundsRemoving = 0

        for round_ in itertools.count():
            changes_add = self._overlapAdd(gamma)
            if changes_add > 0:
                roundsAdding += 1
            changes += changes_add
            if self.verbosity >= 2:
                print "  (r%2s) overlap: adding  : %4d changes, %4d cmtys"%(
                    round_, changes_add, self.q)
            if changes_add == 0 and changes_remove == 0:
                break
            if changes_add == 0:
                changes_remove = self._overlapRemove(gamma)
                if changes_remove > 0:
                    roundsRemoving += 1
                changes += changes_remove
                if self.verbosity >= 2:
                    print "  (r%2s) overlap: removing: %4d changes, %4d cmtys"%(
                        round_, changes_remove, self.q)
            if changes_add == 0 and changes_remove == 0:
                break
            if round_ > 250:
                print "  Exceeding maximum number of rounds"
                break
        return roundsAdding, roundsRemoving, changes
    def _overlapAdd(self, gamma):
        return cmodels.overlapAdd(self._struct_p, gamma)
    def _overlapRemove(self, gamma):
        return cmodels.overlapRemove(self._struct_p, gamma)
    def ovgreedy(self, gamma):
        def A(): return cmodels.overlapAdd(self._struct_p, gamma)
        def B(): return cmodels.overlapRemove(self._struct_p, gamma)
        def C(): return cmodels.combine_sparse_overlap(self._struct_p, gamma)
        def D(): return self.combineSubsets()
        return self.alternate(funcs=(A, B, C, D))
    def ovfree(self, gamma):
        def A(): return cmodels.overlapAdd(self._struct_p, gamma)
        def B(): return cmodels.overlapRemove2(self._struct_p, gamma)
        def C(): return cmodels.combine_sparse_overlap(self._struct_p, gamma)
        def D(): return self.combineSubsets()
        return self.alternate(funcs=(A, B, C, D))

    def combine(self, gamma):
        """Attempt to combine communities if energy decreases."""
        return cmodels.combine(self._struct_p, gamma)
    def combine_supernodes(self, gamma, maxdepth=1):
        """Attempt to combine communities if energy decreases.

        Uses a more general recursive-Graph method."""
        depth = getattr(self, "depth", 0)
        if depth >= maxdepth:
            return 0
        # Can't minimize if there are not enough communities.
        if self.q == 1:
            return 0
        print "==== minimizing subgraph, depth:", depth
        subG = self.supernodeGraph(gamma=gamma, multiplier=1000)
        subG.depth = depth + 1
        # If entire interaction matrix is positive, combining is futile.
        if not numpy.any(subG.imatrix < 0):
            return 0
        changes = subG.greedy(gamma=gamma)
        # No need to try to reload subgraph if no changes.
        if changes == 0:
            return 0
        # A bit of error checking (note: we should never even get here
        # since this is conditinod out above).
        if not numpy.any(subG.imatrix < 0):
            assert changes == 0
        # Reload
        self.loadFromSupernodeGraph(subG)
        print "==== done minimizing subgraph, depth:", depth
#        raw_input(str(changes)+'>')
        return changes
    def combineSubsets(self):
        print "combineSubsets:"
        changes = 0
        for c0 in range(self.Ncmty):
            #print "c0:", c0
            if self.cmtyN[c0] == 0: continue
            for c1 in range(self.Ncmty-1, c0, -1):
                if self.cmtyN[c1] == 0: continue
                if self.cmtyIsSubset(c1, c0):  # is c1 a subset of c0?
                    # if so, delete c1, move these items into c0:
                    print "moving c1=%d (cN=%d) into c0=%s (cN=%d)"%(
                        c1, self.cmtyN[c1], c0, self.cmtyN[c0])
                    #print "c0:", len(cset0)
                    #print "c1:", len(cset1)
                    for n in tuple(self.cmtyContents(c1)):
                        #print "moving n:", n
                        self.cmtyMoveSafe(n, c1, c0)
                    changes += 1
        return changes

    def remap_python(self, check=True):
        """Collapse all communities into the lowest number needed to
        represent everything.

        For example, if we have communities 0, 1, 5, 6 with nonzero
        numbers of particles, this will change us to have communities
        0, 1, 2, 3 only."""
        print "        remapping communities"
        def find_empty_cmty():
            """Start at zero, and find the first (lowest-index) empty
            community.
            """
            for m in range(self.Ncmty):
                if self.cmtyN[m] == 0:
                    # Commented out because this significantly slows things.
                    #assert len([ 1 for i in range(self.N)
                    #             if self.cmty[i] == m  ]) == 0
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
            #for i in range(self.N):
            #    if self.cmty[i] == oldCmty:
            #        self.cmty[i] = newCmty
            self.cmty[self.cmty==oldCmty] = newCmty
            self.cmtyN[newCmty] = self.cmtyN[oldCmty]
            self.cmtyN[oldCmty] = 0
            self.Ncmty -= 1;
        if check:
            self.check()
    def remap_c(self, check=False):
        if self.verbosity >= 2:
            print "        remapping communities (c):",
        changes = cmodels.remap(self._struct_p)
        if self.verbosity >= 2:
            print changes, "changes"
        if check:
            self.check()
    remap = remap_c



    #
    # Methods that deal with visualization
    #
    def get_colormapper(self, colormap_name='gist_rainbow',
                        useGraphColoring=None):
        """Return a colormapper object: colormapper(cmty_id) -> integer_color

        This is a helper method which will return an object which can
        be called with communities and return colors.

        useGraphColoring: if True, color graph cleverly.  If False,
        don't color graph with minimum number of colors.  This is not
        efficient for large graph, default is True for self.N <= 1000,
        False otherwise.
        """
        if useGraphColoring is None:
            if self.N <= 1000:  useGraphColoring = True
            else:               useGraphColoring = False
        if useGraphColoring:
            g = self.supernode_networkx()
            return util.ColorMapper(g, colormap_name=colormap_name)
        else:
            import matplotlib.cm as cm
            import matplotlib.colors as mcolors

            colormap = cm.get_cmap(colormap_name)
            normmap = mcolors.Normalize(vmin=0, vmax=self.Ncmty)
            return list(colormap(normmap(c)) for c in range(self.Ncmty))

    def viz(self, show=True, fname=None):
        """Visualize the detected communities.

        This uses matplotlib to visualize the graphs.  You must have a
        networkx representation of the graph stored at `self._graph`.

        If `show` is True (default), use pyplot.show() to make an
        visual display.  If `fname` is a string, save a copy of the
        graph to that file.
        """
        import matplotlib.pyplot as plt
        #from fitz import interactnow

        colorMap = list(range(self.Ncmty))
        r = random.Random()
        r.seed('color mapping')
        r.shuffle(colorMap)

        g = self._graph.copy()
        for a,b,d in g.edges_iter(data=True):
            weight = d['weight']
            width = 1
            #print a,b,weight
            if self.cmty[g.node[a]['index']] == self.cmty[g.node[b]['index']]:
                newweight = 15
                width = 5
            elif weight < 0:
                newweight = 6
                width = 5
            else:
                newweight = 3
            g.edge[a][b]['weight'] = newweight
            g.edge[a][b]['width'] = width
            #del g.edge[a][b]['weight']

        cmtys = [ self.cmty[d['index']] for (n,d) in g.nodes(data=True) ]
        networkx.draw(g,
                node_color=[colorMap[c] for c in cmtys],
                edge_color=[d.get('color', 'green') for a,b,d in g.edges(data=True)],
                width=[d.get('width', 1) for a,b,d in g.edges(data=True)],
                node_size=900,
                labels = dict((n, "%s (%s)"%(n, c))
                              for (n, c) in zip(g.nodes(), cmtys)),
                #pos=self.coords
                      )
        #plt.margins(tight=True)
        if show:
            plt.show()
        if fname:
            plt.savefig(fname, bbox_inches='tight')
        return g

    def savefig(self, fname, coords=None,
                radii=None, base_radius=0.5, boxsize=None,
                energies=None, energy_radius_ratio=.25,
                nodes='circles',
                hulls=None,
                cmtyColormap=None,
                dpiScale=1,
                explodescale=None,
                **kwargs):
        """Save a copy of layout to `fname`.

        `coords` is a mapping of node index to (x,y) position.  If
        `coords` is not given, use self.coords if it exists.

        `nodes` can be either 'circles' or 'squares', and switches
        between drawing circles for nodes, or drawing a square on the
        background indicating the community of that site (mostly
        useful for square lattices).

        `radii` is a list of relative sizes with the same length as
        `coords`. Smallest particles should be size 1.  `base_radius`
        is a mulitplier for all radii.

        `energies`, if given, will place a small circle in the center
        of each node to indicate its relative energy.  `energies`
        should be an n-length array indicating energy of each
        particle.  `energy_radius_ratio` (default=.25) is a
        multiplying factor indicating how large the energy circle
        should be relative to the node radius.

        `boxsize`, if given, is the (square) box length and used to
        add periodic images of particles when needed.  If boxsize is
        not given, use self.boxsize if it exists.

        `hulls`, if true, will draw convex hulls of communities on the
        background.

        `dpiScale`, a float, is a scale factor for total pixel count
        of the final image.  This makes images with greater resolution.

        `explodescale`: if given, this is an amount to 'explode' the
        graph to visualize overlapping communities.
        """
        if self.verbosity > 0:
            print "Savefig:", fname
        import matplotlib.figure
        import matplotlib.backends.backend_agg
        from matplotlib.patches import Rectangle, Polygon
        from matplotlib.patches import Circle
        #from matplotlib.patches import CirclePolygon as Circle
        import matplotlib.cm as cm
        import matplotlib.colors as colors

        if cmtyColormap is None:
            cmtyColormap = self.get_colormapper()

        if coords is None and hasattr(self, "coords"):
            coords = self.coords
        if coords is None:
            raise Exception("If coords is not given, self.coords must exist.")
        if boxsize is None and hasattr(self, "boxsize"):
            boxsize = self.boxsize
        if boxsize is not None:
            coords = coords % boxsize

        #f = matplotlib.backends.backend_agg.Figure()
        fig = matplotlib.figure.Figure()
        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
        ax = fig.add_subplot(111, aspect='equal')

        if radii is not None:
            radii=radii*base_radius
        elif hasattr(self, 'radii'):
            radii=self.radii*base_radius
        else:
            radii=numpy.ones(self.N)*base_radius

        if explodescale is not None:
            escale = explodescale
        else:
            escale = 1

        if explodescale is None:
            for n in range(self.N):
                if nodes == 'circles':
                    p = Circle(coords[n], radius=radii[n], axes=ax,
                               #color=colormap.to_rgba(n),
                               color=cmtyColormap[self.cmty[n]],
                               **kwargs)
                elif nodes == 'squares':
                    p = Rectangle(coords[n]-(base_radius,base_radius),
                                  width=2*base_radius, height=2*base_radius,
                                  axes=ax,
                                  color=cmtyColormap[self.cmty[n]],
                                  **kwargs)
                ax.add_patch(p)
                if boxsize is not None and nodes == 'circles':
                    greaterthanl=( coords[n] + radii[n] > boxsize  )
                    lessthanzero=( coords[n] - radii[n] < 0        )
                    if greaterthanl.any() or lessthanzero.any():
                        p = Circle(
                           coords[n]-boxsize*greaterthanl+boxsize*lessthanzero,
                            radius=radii[n],
                            color=cmtyColormap[self.cmty[n]],
                            **kwargs)
                        ax.add_patch(p)

        else:
            # Exploded community map
            for c in range(self.Ncmty):
                if self.cmtyN[c] == 0: continue
                #coords = coords % boxsize
                cmtynodes = list(self.cmtyContents(c))
                cmtycoords = coords[cmtynodes]
                mean_position = util.mean_position(
                    cmtycoords[[ i for i,n in enumerate(cmtynodes)
                                if self.cmty[n]==c ]],
                    boxsize)
                cmtycoords = util.wrap_dists(cmtycoords-mean_position,boxsize)\
                             + mean_position

                cmtycoords += (escale-1) * mean_position
                cmtycoords %= escale*boxsize

                #pos_addition = (escale-1) * mean_position
                for i, n in enumerate(cmtynodes):
                    if self.cmty[n] == c:
                        alpha = 1
                    else:
                        alpha = .25
                    p = Circle(cmtycoords[i],
                               radius=radii[n], axes=ax,
                               #color=colormap.to_rgba(n),
                               color=cmtyColormap[c],
                               alpha=alpha,
                               linewidth=0,
                               **kwargs)
                    ax.add_patch(p)


        if energies is not None:
            e_colormap = cm.get_cmap('RdYlBu')
            e_norm = colors.Normalize()
            e_norm.autoscale(energies)
            for n in range(self.N):
                cir = Circle(coords[n], radius=radii[n]*energy_radius_ratio,
                             axes=ax,
                             color=e_colormap(e_norm(energies[n])),
                             zorder=1,
                             **kwargs)
                ax.add_patch(cir)

        if hulls is not None:
            #from support.convex_hull import convex_hull
            #from support.convexhull import convexHull
            from support.convexhull import convex_hull
            for c in self.cmtys():
                ns = self.cmtyContents(c)
                points = [ coords[n] for n in ns ]
                points = numpy.asarray(points)
                #if len(points) > 5:
                #    print points
                #    points = convex_hull(points.T)
                points = convex_hull(tuple(p) for p in points)
                # Detect if any points are too long
                if boxsize is not None:
                    pts1 = points[:]
                    pts2 = points[1:] + points[0:1]
                    d = numpy.subtract(pts1, pts2)
                    # We don't sum along axis=-1, instead we do it per-axis
                    d /= boxsize # inplace
                    #print d
                    if numpy.abs(d).max() > .5:
                        continue
                # Do actual plotting
                p = Polygon(points, alpha=.25, color=cmtyColormap[c],
                            zorder=-2,
                            axes=ax)
                ax.add_patch(p)


        ax.autoscale_view(tight=True)
        if boxsize is not None:
            if isinstance(boxsize, numbers.Number):
                ax.set_xlim(0,boxsize*escale)
                ax.set_ylim(0,boxsize*escale)
            else:
                ax.set_xlim(0,boxsize[0]*escale)
                ax.set_ylim(0,boxsize[1]*escale)
        canvas.print_figure(fname, dpi=fig.get_dpi()*dpiScale*escale,
                            bbox_inches='tight')
        if self.verbosity > 0:
            print "Done saving"







if __name__ == "__main__":
    command = None
    if len(sys.argv) > 1:
        command = sys.argv[1]

    if command == "test":

        random.seed(169)
        G = random_graph(size=20, layout=True)
        random.seed(time.time())
        G.cmtyCreate()
        G.check()
        #G.test()
        G.greedy(gamma=1.0)
        #G.remapCommunities()
        #G.remapCommunities()
        #G.check()
        G.viz()

    if command == "multi":
        Gs = [ random_graph(size=20) for _ in range(10) ]
        M = MultiResolution(low=.01, high=100)
        M.do(Gs=Gs, trials=3)
        M.print_()
