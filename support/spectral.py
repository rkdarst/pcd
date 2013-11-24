# Richard Darst, August 2013

import cmath
import collections
import math
import numpy
#from numpy import linalg
import random
import scipy.linalg
import scipy.sparse
from pcd import cmty

import networkx

class Spectrum(object):
    check_empty = True
    sparse = False
    verbosity = 1
    def __init__(self, **kwargs):
        for k, v in kwargs.iteritems():
            if not hasattr(self, k): raise ValueError("Unknown option '%s'"%k)
            setattr(self, k, v)

    def _init(self):
        """Internal initialization.

        Create sorted eigenvalues/eigenvectors.  We always sort
        eigenvectors/eigenvalues by real value, ascending.  So 0 is
        the lowest, and -1 is the greatest."""
        ev_sorted = sorted(((ev, i) for i,ev in enumerate(self.ev)), key=lambda (ev,i): ev.real)
        self.ev_sorted = tuple(x for x,i in ev_sorted)
        self.ev_sorted_vecindex = dict((j, i)
                                        for j, (ev,i) in enumerate(ev_sorted))
        #from fitz import interactnow
    def _calc_ev(self, M, k=None, herm=False):
        """Calculate eigenvalues/eigenvectors.

        if herm is True, use a hermitian EV calculator.

        Returns (ev, evec)"""
        if self.verbosity >= 1:
            if k:
                print "computing %d EVs (matrix size %s)..."%(k, M.shape,)
            else:
                print "computing EVs (matrix size %s)..."%(M.shape,)
        import os
        utime0 = os.times()[0] # utime
        if k is None:
            if herm:
                # Assume hermitian matrix
                ev, evec = scipy.linalg.eigh(M)
            else:
                ev, evec = scipy.linalg.eig(M)
        else:
            # Calculate only 'k' eigenvectors.  Only negative indexes
            # make sense in this case.  LR=largest real part.
            def eigs_sparse(M, k):
                #M = M.astype(float)
                #return scipy.linalg.eig(M)
                try:
                    return scipy.sparse.linalg.eigs(M, k=k, which='LR')
                except AttributeError:
                    # Older scipy: should be removed once possible
                    import scipy.sparse.linalg.eigen.arpack as arpack
                    return arpack.eigen(M, k=k, which='LR')
            if herm:
                ev, evec = scipy.sparse.linalg.eigsh(M, k, which='LR')
            else:
                #ev, evec = scipy.sparse.linalg.eigs(M, k, which='LR')
                ev, evec = eigs_sparse(M, k)
        self.time_ev = os.times()[0]-utime0
        if self.verbosity >= 1:
            print "done, took %f seconds."%(self.time_ev)
        return ev, evec

    def get_ev_index(self, i):
        """Return the i'th eigenvalue (ascending)"""
        return self.ev_sorted[i]
    def get_evec_index(self, i, tol=1e-15):
        """Return the i'th eigenvector (corresponding to i'th evalue ascending)"""
        i = i%len(self.ev_sorted)
        evec = self.evec[:, self.ev_sorted_vecindex[i]]
        for i, x in enumerate(evec):
            if abs(x) < tol:
                evec[i] = 0.0
        return evec
    def get_evec_from_ev(self, ev):
        for i in range(len(self.ev)):
            if self.ev[i] == ev:
                return self.evec[:, i]
        raise ValueError("Eigenvalue %s not found"%ev)
    def get_distinct_ev(self, rounding=5):
        pass
    def get_real_ev(self, tol=1e-6):
        """Get an iterator over all eigenvalues, sorted ascending"""
        for ev in self.ev:
            if abs(ev.imag) < tol:
                yield ev

    def cmtys(self, q=None):
        """Detect communities.

        q: if given, assume this number of communities.  Otherwise,
        use self.q().

        self.check_empty: if true (default), raise a warning if we detect some empty communities."""
        from pcd import cmty
        if q is None and not hasattr(self, 'q'):
            raise ValueError('Must supply q to detect.')
        if q is None:
            q = self.q()
        # q==1 is one giant community.  We have no vectors to classify.
        if q == 1:
            nodes = self._all_nodes()
            if len(nodes) != self.N:
                print "We seem to not spanned the graph.  Do all nodes have an incoming edge?"
            cmtynodes={0:nodes}
            if hasattr(self, 'self.readd_nodes'):
                cmtynodes = self.readd_nodes(cmtynodes)
            cmtys = cmty.Communities(cmtynodes)
            return cmtys
        # Create initial classification vectors.  To detect q
        # communities, we need q-1 EVs.
        vecs = [ ]
        for index in range(q-1):
            # Some methods take EVs [0,q-1], some take [-1, -q-1], some take [-2, -q-1]
            ev_index = self._classification_index(index)
            vecs.append(self.cmty_classify(ev_index))

        # Create a list of all nodes
        nodes_set = set()
        for cl in vecs:
            nodes_set.update(cl.keys())
        nodes = sorted(nodes_set)
        #nodes = sorted(self._all_nodes())
        if len(nodes) != self.N:
            print "We seem to not spanned the graph.  Do all nodes have an incoming edge?"
        # Convert list of dictionaries into a list of tuples
        vecs = [
            tuple(cl[n] for n in nodes)
            for cl in vecs ]
        vecs = numpy.asarray(vecs).T

        # Do the clustering
        import scipy.cluster.vq
        def do_classify():
            """Function to do the clustering.  Within a function since
            we may apply it multiple times."""
            centroids, distortion = scipy.cluster.vq.kmeans(
                vecs,
                #scipy.cluster.vq.whiten(vecs),
                q,
                iter=10 # run k-means 10 times
                )
            #centroids, labels = scipy.cluster.vq.kmeans2(
            #    scipy.cluster.vq.whiten(vecs),
            #    q,
            #    )

            # Classify into the communities
            def dist(a, b):
                return numpy.linalg.norm(a-b)
            cmtynodes = dict((c, set()) for c in range(q))
            for i, node in enumerate(nodes):
                vec = vecs[i]
                dists = [ dist(vec, cent) for cent in centroids ]
                c = numpy.argmin(dists)
                cmtynodes[c].add(node)
            return  cmtynodes

        # Try the classification 10 times.  Stop when we have a
        # classification that has non-empty communities.  There is a
        # high chance to have empty communities.
        for i in range(10):
            cmtynodes = do_classify()
            if all(len(ns)>0 for ns in cmtynodes.itervalues()):
                break
        # Print a warning if we have empty communities
        if self.check_empty:
            def errormsg():
                if hasattr(self, 'q_circle'):
                    return "Detected some empty communities (q=%s,%s detected, q=%s sought) (%s)."%(
                        self.q_bulk(), self.q_circle(), q,
                        sorted([len(ns) for ns in cmtynodes.itervalues()]))
                return "Detected some empty communities: (q=%s sought, %s found) (%s)"%(
                    q, sum(1 for ns in cmtynodes.itervalues() if len(ns)>0),
                    sorted([len(ns) for ns in cmtynodes.itervalues()]))
            assert all(len(ns)>0 for ns in cmtynodes.itervalues()), errormsg()
        # And then remove the communities with no nodes.
        for cname in list(cmtynodes.keys()):
            if len(cmtynodes[cname]) == 0:
                del cmtynodes[cname]
        # re-add dangling tree nodes:
        if hasattr(self, 'self.readd_nodes'):
            cmtynodes = self.readd_nodes(cmtynodes)
        # Create and return final communities object.
        cmtys = cmty.Communities(cmtynodes=cmtynodes)
        #from fitz import interact ; interact.interact()
        return cmtys



class _SimpleSpectrum(Spectrum):
    #operator = networkx.adjacency_matrix
    def __init__(self, g, **kwargs):
        super(_SimpleSpectrum, self).__init__(**kwargs)

        self.N = len(g)  # number of nodes
        nodes = self.nodes = sorted(g.nodes())
        M = self.get_matrix(g, nodelist=nodes)

        self.ev, self.evec = self._calc_ev(M, herm=True)
        # Calculate
        self._calc_ev(M, herm=True)
        self._init()
        #E = g.number_of_edges()
    def _all_nodes(self):
        return self.nodes
    def cmty_classify(self, index, signs=False):
        """Classify nodes into communities using the specified eigenvector."""
        evec = self.get_evec_index(index)
        return dict(zip(self.nodes, evec))

class Adjacency(_SimpleSpectrum):
    def get_matrix(self, g, nodelist=None):
        return networkx.to_numpy_matrix(g, nodelist=nodelist)
    def relevant_ev(self, k):
        """We want the largest eigenvalues"""
        return [self.get_ev_index(i) for i in range(-1, -k-1, -1) ]
    def _classification_index(self, k):
        raise NotImplementedError("Check this first.")
        return -k-1  # 0: -1   1: -2, etc.
class Laplacian(_SimpleSpectrum):
    def get_matrix(self, g, nodelist=None):
        return networkx.laplacian_matrix(g, nodelist=nodelist)
    def relevant_ev(self, k):
        """We want the smallest eigenvalues"""
        return [self.get_ev_index(i) for i in range(0, k, 1) ]
    def _classification_index(self, k):
        return k+1  # 0: 1   1: 2, etc.
class Modularity(_SimpleSpectrum):
    def get_matrix(self, g, nodelist=None):
        adj = networkx.to_numpy_matrix(g, nodelist=nodelist)
        d = [g.degree(n) for n in nodelist]
        degree_products = numpy.outer(d, d)
        B = adj - degree_products/float(2.*g.number_of_edges())
        return B
    def _classification_index(self, k):
        raise NotImplementedError("Check this first.")
        return -k-1  # 0: -1   1: -2, etc.
    def relevant_ev(self, k):
        """We want the largest eigenvalues"""
        return [self.get_ev_index(i) for i in range(-1, -k-1, -1) ]



class _2Walk(Spectrum):
    check_num_ev = True
    def __init__(self, g, num_ev=None, **kwargs):
        for k, v in kwargs.iteritems():
            if not hasattr(self, k): raise ValueError("Unknown option '%s'"%k)
            setattr(self, k, v)

        # Return a copy of the graph suitable for 2-walk analysis
        # (without dangling trees).
        if self.verbosity >= 0:
            print "old graph size n, e:", len(g), g.number_of_edges()
        g = self.remove_dangling_trees(g)
        if self.verbosity >= 0:
            print "new graph size n, e:", len(g), g.number_of_edges()
        self.average_degree = 2. * g.number_of_edges() / g.number_of_nodes()
        self.average_dDd1 = numpy.mean([ g.degree(n)/float(g.degree(n)-1) for n in g.nodes_iter()])

        self.E = E = g.number_of_edges()
        self.N = len(g)  # number of nodes
        if num_ev is not None:
            self.sparse = True


        # Construct our ordering of edges to go into the matrix.
        B_edgelist = self.B_edgelist = [ ]
        #print n_neigh = [ ]
        for node, adj_dict in g.adjacency_iter():
            B_edgelist.extend((node, key) for key in adj_dict if key!=node)
            #for neigh in g.neighbors_iter(node):
            #    if (node, neigh) not in B_edgelist:
            #        raise
        assert len(B_edgelist) == 2*E - 2*len(g.nodes_with_selfloops())
        B_edgelist_lookup = dict((edge, i) for i,edge in enumerate(B_edgelist))
        #print len(B_edgelist_lookup), len(B_edgelist)
        assert len(B_edgelist_lookup) == len(B_edgelist)

        # Construct our actual matrix.  For each starting node, look
        # at neighbors, and then neighbors of neighbors.  This gives
        # us the 2-paths in the graph.
        # float32 required for ARPACK.
        if not self.sparse:
            B = numpy.zeros(shape=(2*E, 2*E), dtype=numpy.float32)
        else:
            #B = sparse.coo_matrix((2*E, 2*E), dtype=int)
            B = scipy.sparse.lil_matrix((2*E, 2*E), dtype=numpy.float32)
        for n1 in g.nodes_iter():
            for n2 in g.neighbors(n1):
                if n1 == n2: continue   # Exclude self-loops
                for n3 in g.neighbors(n2):
                    if n2 == n3: continue   # Exclude self-loops
                    if n1 == n3: continue   # Non backtracking
                    i1 = B_edgelist_lookup[(n1, n2)]
                    i2 = B_edgelist_lookup[(n2, n3)]
                    if not self.degree_corrected:
                        # self.degree_corrected for the Flow method
                        B[i1, i2] = 1
                    else:
                        B[i1, i2] = 1./(g.degree(n2)-1)
        if self.sparse:
            B = scipy.sparse.csr_matrix(B)
        # The matrix is not expected to be symmetric.
        #print numpy.max(B-B.transpose())
        #assert not numpy.any(numpy.abs(B - B.transpose()) > 1e-6)
        #from fitz import interact ; interact.interact()

        #self.M = B

        self.ev, self.evec = self._calc_ev(B, k=num_ev)
        self._init()
    def remove_dangling_trees(self, g):
        """Remove dangling trees and save state.

        To run walk spectra, we should remove all dangling trees from
        the graph.  This takes a graph, copies it, removes dangling
        trees, and saves the state on self.  Use the function
        `readd_nodes` to re-add nodes to a 'cmtynodes' object."""
        g = g.copy()
        self.allnodes = set(g.nodes())
        # change (the copy of) g in place.
        self.root_map = remove_trees(g)
        return g
    def readd_nodes(self, cmtynodes):
        """Re-add dangling tree nodes to communities detected.

        This modifies cmtynodes IN PLACE..."""
        new_cmtynodes = add_trees_to_cmtys(cmtynodes, self.root_map, self.allnodes)
        return new_cmtynodes
    def _all_nodes(self):
        """All nodes that should be detectable.

        This is used for q=1"""
        return set(to for from_, to in self.B_edgelist)


    def q(self, *args, **kwargs):
        return getattr(self, 'q_'+self.q_mode)(*args, **kwargs)

    def cmtys_q2(self, complete_graph=True):
        """Detect communities, assuming q=2.

        If complete_graph is True (default), then re-add extra nodes
        that were removed as part of dangling trees."""
        if self.q() != 2:
            print "Warning(\"We didn't detect two communities (%s actual)"%self.q()

        signs = self.cmty_classify(-2)

        cmtynodes = {0:set(), 1:set()}   # this case (q=2) case limited to two communities
        for node, sum_ in signs.iteritems():
            if sum_ < 0:
                cmtynodes[0].add(node)
            elif sum_ > 0:
                cmtynodes[1].add(node)
            else:
                # random placement in both of them... any better ideas?
                cmtynodes[random.choice((0,1))].add(node)

        #from fitz import interact ; interact.interact()
        if sum(len(ns) for ns in cmtynodes.itervalues()) != self.N:
            print "We seem to not spanned the graph.  Do all nodes have an incoming edge?"
        # re-add dangling tree nodes:
        if complete_graph:
            cmtynodes = self.readd_nodes(cmtynodes)
        # Make communities object.
        cmtys = cmty.Communities(cmtynodes=cmtynodes)
        return cmtys

    def _classification_index(self, k):
        """The kth relevant EV for clustering.

        For these 2-walk matrices, the first EV to use for clustering
        is the second highest, or -2 in Python syntax."""
        return -k-2  # 0: -2  1: -3, etc

    def q_circle(self):
        """Returns detected number of communities.

        Detected number of communities is the number of (real)
        eigenvalues greater than sqrt(self.c)."""
        limit = self.circle_radius
        q = len(tuple(ev for ev in self.get_real_ev() if ev > limit))
        num_ev = self.evec.shape[1]
        if self.check_num_ev:
            if q >= num_ev:
                raise ValueError("We did not compute enough eigenvalues so "
                                 "probably did not calculate q correctly. "
                                 "(q=%s, num_ev=%s)"%(q, num_ev))

        return q
    def q_bulk(self, tol=1e-6):
        """Detect number of communities based on bulk measures.

        Returns the number of real eigenvalues which have a magnitude
        greater than (the real part of the complex eigenvalue with the
        greatest real magnitude).

        tol: imaginary parts less than this tolerance (1e-6 default)
        are considered zero.
        """
        q = 0
        for i in range(len(self.ev_sorted)-1, -1, -1):
            ev = self.ev_sorted[i]
            if abs(ev.imag) < tol:
                q += 1
            else:
                break
        num_ev = self.evec.shape[1]
        if self.check_num_ev:
            if q >= num_ev:
                raise ValueError("We did not compute enough eigenvalues so "
                                 "probably did not calculate q correctly. "
                                 "(q=%s, num_ev=%s)"%(q, num_ev))
        #from fitz import interactnow
        return q


    def cmty_classify(self, index, signs=False):
        """Classify nodes into communities along specified eigenvalue.

        Returns a dict mapping node->sum of incoming character along
        the eigenvector 'index'.  Index should probably be negative,
        from -2 to -q inclusive."""

        evec = self.get_evec_index(index)
        node_sign_sum = collections.defaultdict(int)
        for i, edge in enumerate(self.B_edgelist):
            from_, to_ = edge
            if signs:
                if evec[i].real > 0:    sign = 1
                elif evec[i].real < 0:  sign = -1
                else:                   sign = 0
                node_sign_sum[to_] += sign
            else:
                node_sign_sum[to_] += evec[i].real
        return dict(node_sign_sum)

    def plot_spectrum(self):
        """Make a scatter plot (in the complex plane) of eigenvalues"""
        import pylab
        #from matplotlib import pyplot

        pylab.scatter(self.ev.real, self.ev.imag)

        circ = pylab.Circle((0,0), self.circle_radius, fill=False)
        pylab.gca().add_artist(circ)
    def plot_detected(self, axes=None, cmtys=None):
        if axes is None:
            ax1 = -2
            ax2 = -3
        else:
            ax1, ax2 = axes
        import matplotlib
        from matplotlib import pylab
        nodes, sums_ax1 = zip(*sorted(self.cmty_classify(ax1).items()))
        nodes, sums_ax2 = zip(*sorted(self.cmty_classify(ax2).items()))

        if cmtys:
            cm = matplotlib.cm.get_cmap('jet', lut=cmtys.q)
            nodecmtys = cmtys.nodecmtys_onetoone()
            colors = cm(tuple(nodecmtys[n] for n in nodes))
        else:
            colors=cmtys
        #pylab.scatter(sums, nodes, color=cm(nodes))
        pylab.scatter(sums_ax1, sums_ax2, color=colors)


class NonBacktrackingWalk(_2Walk):
    degree_corrected = False
    def c(self):
        """Return detected average degree = largest_ev.real"""
        assert self.ev_sorted[-1].imag < 1e-9
        return self.ev_sorted[-1].real
    @property
    def circle_radius(self):
        return math.sqrt(self.c())
    q_mode = 'circle'
    #q = _2Walk.q_circle
# compatability alias
KMMNEZZ = NonBacktrackingWalk


class Flow(_2Walk):
    #circle_radius = 1
    degree_corrected = True
    q_mode = 'bulk'
    #q = _2Walk.q_bulk
    @property
    def circle_radius(self):
        # this should be the below, but we need to calculate more averages.
        #math.sqrt( mean(d/(d-1)) / mean(d) )
        return math.sqrt(self.average_dDd1 / self.average_degree)


def remove_trees(g):
    # remove trees in-place.
    root_map = { }   # map of node -> community of it's root.
    nodes = g.nodes() # must be list
    while len(nodes) > 0:
        n = nodes.pop()
        if g.degree(n) > 1:
            continue
        # this will fail if degree = 0
        neighbor = g.neighbors(n)
        assert len(neighbor) == 1
        neighbor = neighbor[0]
        # update root map.
        root_map[n] = neighbor
        #if neighbor in root_map:
        # remove from graph.
        g.remove_node(n)
        if g.degree(neighbor) == 1:
            nodes.append(neighbor)
    return root_map
def find_root(root_map, n):
    # n not in root_map: return n
    # n in root map: recursive lookup, return the final value.
    while n in root_map:
        n = root_map[n]
    return n
def add_trees_to_cmtys(cmtynodes, root_map, allnodes):
    cmtys = cmty.Communities(cmtynodes)
    nodecmtys = cmtys.nodecmtys_onetoone()
    for n in allnodes:
        if n in root_map:
            assert n not in nodecmtys
            root_node = find_root(root_map, n)
            c = nodecmtys[root_node]
            cmtynodes[c].add(n)
    return cmtynodes
