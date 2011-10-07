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

import cmodels
import util

log2 = lambda x: math.log(x, 2)

NO_CMTY = -1


class Graph(cmodels._cobj, object):
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

    def __init__(self, N, randomize=True, rmatrix=False, overlap=False):
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
        self._fillStruct(N)
        self.N = N
        self.cmtyCreate(randomize=randomize)
        self.oneToOne = 1
        if rmatrix:
            self._allocRmatrix()
        if overlap:
            self.setOverlap(True)
    def _fillStruct(self, N=None):
        """Fill C structure."""
        cmodels._cobj._allocStruct(self, cmodels.cGraph)

        self._allocArray("cmtyll", shape=[N]*2)
        self._allocArray("cmtyl", shape=N, dtype=ctypes.c_void_p)
        self._allocArrayPointers(self.cmtyl, self.cmtyll)
        self._allocArray("cmty",  shape=N)
        self._allocArray("cmtyN", shape=N)
        self._allocArray("randomOrder", shape=N)
        self._allocArray("randomOrder2", shape=N)
        self.randomOrder[:]  = numpy.arange(N)
        self.randomOrder2[:] = numpy.arange(N)
        self._allocArray("imatrix", shape=(N, N))
    def _allocRmatrix(self):
        """Allocate a matrix for repulsive interactions"""
        self._allocArray("rmatrix", shape=(self.N, self.N))

    def __getstate__(self):
        """Get state for pickling"""
        state = self.__dict__.copy()
        del state['_struct']
        del state['_struct_p']
        state['__extra_attrs'] = {'N':self.N, 'oneToOne':self.oneToOne }
        #for name, type_ in self._struct._fields_:
        #    if name in state  and isinstance(state[name], numpy.ndarray):
        #        setattr(self, name, state[name])
        #        del state[name]
        return state
    def __setstate__(self, state):
        """Set state when unpickling"""
        self.__init__(state['__extra_attrs']['N'])
        #self._fillStruct(n=state['n'])
        for name, type_ in self._struct._fields_:
            #print name
            if name in state and isinstance(state[name], numpy.ndarray):
                self._allocArray(name, array=state[name])
                del state[name]
        self._allocArrayPointers(self.cmtyl, self.cmtyll)
        for k, v in state['__extra_attrs'].items():
            setattr(self, k, v)
        del state['__extra_attrs']
        self.__dict__.update(state)
    def copy(self):
        """Return a copy of the object (not sharing data arrays).
        """
        # Right now we turn this into a string and back, in order to
        # make actual deep copies of all the arrays.  FIXME.
        new = pickle.loads(pickle.dumps(self, -1))
        #new._allocArray('imatrix', array=self.imatrix)
        return new
    def getcmtystate(self):
        """Return an object representing the communities.

        This can be pickeled and re-loaded into this object later.

        See also .cmtyDict().
        """
        state = { }
        state['version'] = 0
        state['oneToOne'] = self.oneToOne
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
            self.cmty[:] = state['cmtyList']
        else:
            raise Exception("Unknown version of state to load communities.")
        #self.check() # can't check non-oneToOne (yet)
    def hash(self):
        """Hash of state of self.

        This is not a real Python hash, since this is a mutable object."""
        return hash((
            self.oneToOne,
            tuple(self.cmty.flat),
            tuple(self.cmtyN.flat),
            tuple(tuple(self.cmtyContents(c)) for c in self.cmtys()),
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
                     layout=None, randomize=True,
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
        G._layout = layout

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
        self.remapCommunities()
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
        assert self.q == len(set(self.cmty))
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
        self.cmtyListInit()
    def cmtyListInit(self):
        """Initialize the efficient community lists.

        The community list (self.cmtyll[community_index, ...]) is a
        data structure optimized for the community lookups.  This
        function fills self.cmtyll from self.cmty[i].
        """
        print "cmtyListInit"
        cmodels.cmtyListInit(self._struct_p)
    def cmtyListClear(self):
        """Clear all particles from all communities"""
        self.cmty[:] = NO_CMTY
        self.cmtyN[:] = 0
        self.Ncmty = 0
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
    def cmtySet(self, n, c):
        """Set the community of a particle."""
        raise NotImplementedError("cmtySet not implemented yet.")
    def cmtyContents(self, c):
        """Array of all nodes in community c."""
        return self.cmtyll[c, :self.cmtyN[c]]
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



    #
    # Methods that deal with energies or other thermodynamic measures
    # or statistics of the current state of the graph.
    #
    def energy(self, gamma):
        """Return total energy of the graph."""
        return cmodels.energy(self._struct_p, gamma)
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
        community_sizes=[ self.cmtyN[c] for c in range(self.Ncmty)
                          if self.cmtyN[c] > 0 ]
        for size in community_sizes:
            if size in n_counts.keys():
                n_counts[size]=n_counts[size]+1
            else:
                n_counts[size]=1
        return n_counts



    #
    # Methods that deal with minimization/optimization
    #
    def minimize_trials(self, gamma, trials):
        """Minimize system using .minimize() and `trials` trials.

        This will minimize `trials` number of times, and set the final
        graph to be the state that had the best energy of all the
        trials.
        """
        if self.verbosity > 0:
            print "Minimizing with %d trials (gamma=%f)"%(trials,gamma)
        minE = float('inf')
        minCmtyState = self.getcmtystate()

        for i in range(trials):
            self.cmtyCreate() # randomizes it
            changes = self.minimize(gamma)
            thisE = self.energy(gamma)
            if thisE < minE:
                minE = thisE
                minCmtyState = self.getcmtystate()

        self.setcmtystate(minCmtyState)
        if self.verbosity > 0:
            print "... Done minimizing with %d trials (gamma=%f)"%(
                                                                  trials,gamma)
        return changes

    def minimize(self, gamma):
        """Minimize the communities at a certain gamma.

        This function requires a good starting configuration.  To do
        that, use self.cmtyCreate() first.
        """
        print "beginning minimization (n=%s, gamma=%s)"%(self.N, gamma)
        changes = 0

        changesCombining = None
        for round_ in itertools.count():
            self._gen_random_order()
            changesMoving = self._minimize(gamma=gamma)
            changes += changesMoving
            if self.verbosity > 1:
                print "  (r%2s) cmtys, changes: %4d %4d"%(round_, self.q,
                                                          changesMoving)

            if changesMoving == 0 and changesCombining == 0:
                break

            # If we got no changes, then try combining communities
            # before breaking out.
            if changesMoving == 0:
                changesCombining = self.combine_cmtys(gamma=gamma)
                #changesCombining = self.combine_cmtys_supernodes(gamma=gamma)
                changes += changesCombining
                if self.verbosity > 1:
                    print "  (r%2s) cmtys, changes: %4d %4d"%(
                                          round_, self.q, changesCombining), \
                                                  "(<-- combining communities)"

                self.remapCommunities(check=False)
            # If we have no changes in regular and combinations, escape
            if changesMoving == 0 and changesCombining == 0:
                break
            if round_ > 250:
                print "  Exceeding maximum number of rounds."
                break
        #print set(self.cmty),
        if self._use_overlap:
            self.overlapMinimize(gamma)
        return changes
    def _minimize(self, gamma):
        return cmodels.minimize(self._struct_p, gamma)
    def overlapMinimize(self, gamma):
        """Attempting to add particles to overlapping communities.
        """
        if self.verbosity > 0:
            print "Minimizing by trying overlaps."
        changes = 0
        changes_add = 0
        changes_remove = None
        for round_ in itertools.count():
            changes_add = self._overlapMinimize_add(gamma)
            changes += changes_add
            print "  (r%2s) overlap: adding  : %4d"%(round_, changes_add)
            if changes_add == 0 and changes_remove == 0:
                break
            if changes_add == 0:
                changes_remove = self._overlapMinimize_remove(gamma)
                changes += changes_remove
                print "  (r%2s) overlap: removing: %4d"%(
                                                     round_, changes_remove)
            if changes_add == 0 and changes_remove == 0:
                break
            if round_ > 250:
                print "  Exceeding maximum number of rounds"
                break
        return changes
    def _overlapMinimize_add(self, gamma):
        return cmodels.overlapMinimize_add(self._struct_p, gamma)
    def _overlapMinimize_remove(self, gamma):
        return cmodels.overlapMinimize_remove(self._struct_p, gamma)

    def combine_cmtys(self, gamma):
        """Attempt to combine communities if energy decreases."""
        return cmodels.combine_cmtys(self._struct_p, gamma)
    def combine_cmtys_supernodes(self, gamma, maxdepth=1):
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
        changes = subG.minimize(gamma=gamma)
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

    def remapCommunities_python(self, check=True):
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
    def remapCommunities_c(self, check=False):
        if self.verbosity > 0:
            print "        remapping communities (c):",
        changes = cmodels.remap_cmtys(self._struct_p)
        if self.verbosity > 0:
            print changes, "changes"
        if check:
            self.check()
    remapCommunities = remapCommunities_c



    #
    # Methods that deal with visualization
    #
    def get_colormapper(self, colormap_name='gist_rainbow'):
        """Return a colormapper object: colormapper(cmty_id) -> integer_color

        This is a helper method which will return an object which can
        be called with communities and return colors.
        """
        g = self.supernode_networkx()
        return util.ColorMapper(g, colormap_name=colormap_name)
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
                pos=self._layout)
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
        """
        if self.verbosity > 0:
            print "Savefig:", fname
        import matplotlib.figure
        import matplotlib.backends.backend_agg
        from matplotlib.patches import Circle, Rectangle, Polygon
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

        #f = matplotlib.backends.backend_agg.Figure()
        fig = matplotlib.figure.Figure()
        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
        ax = fig.add_subplot(111, aspect='equal')

        if type(radii)==type(None):
            radii=numpy.ones(self.N)*base_radius
        else:
            radii=radii*base_radius

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
                ax.set_xlim(0,boxsize)
                ax.set_ylim(0,boxsize)
            else:
                ax.set_xlim(0,boxsize[0])
                ax.set_ylim(0,boxsize[1])
        canvas.print_figure(fname, bbox_inches='tight')
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
        G.minimize(gamma=1.0)
        #G.remapCommunities()
        #G.remapCommunities()
        #G.check()
        G.viz()

    if command == "multi":
        Gs = [ random_graph(size=20) for _ in range(10) ]
        M = MultiResolution(low=.01, high=100)
        M.do(Gs=Gs, trials=3)
        M.print_()
