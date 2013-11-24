# Richard Darst, December 2012

"""Class representing community object.

Communities represent groupings of nodes.  These groups of nodes can
optionally have some meaningful name.  This module contains classes
for represting community objects (stored in different forms, such as a
dict or file on disk).

Iterators
=========

The most general method of accessing communities is by iterating
through them.  This allows on to implement streaming schemes, for
example, from a file.  This allows partial loading of the community
information from files, allowing one to access larger scale networks.

The basic iterator methods are:
  cmtys.iteritems()
    iterate over (name, set(nodes)) pairs
    example: `for name, nodes in cmtys.iteritems(): print len(nodes)`
  cmtys.itervalues()
    like .iteritems() but only iterates over nodes
    example: `for nodes in cmtys.itervalues(): print len(nodes)`
  cmtys.iterkeys()
    like .iteritems() but only iterates over community  names
    example: `for name in cmtys.itervalues(): print name`
  len(cmtys)
    number of communities.  Note that this may require internally
    iterating through all communities to calculate this, thus it may
    not be fast process.  Some metheds will cache this value.

You will notice the similarity of this to dict methods.  This is by
design, so that simple functions can take either dictionaries or
Communities objects.

Normally, you would think the universe of nodes is the union of all
communities.  This would mean that the structure always covers the
entire graph.  Sometimes, the community sturcture doesn't span the
entire universe of nodes, and we need a way to represent that.  For
Communities objects, you can pass a nodes= parameter to the
constructor.  This distinction is relevant to the difference between
.nodes and .nodes_spanned().



Methods implemented using only iterators
========================================

Given the four iterators above, many things can be calculated
automatically.  These methods are implemented in the _CommunitiesBase
object, which groups all methods which take only the requirements
above (for the most part, this is not quite true yet).

cmtys.N: int
    Number of nodes
cmtys.q: int
    Number of communities
cmtys.to_dict(): dict
    All communities.  Like dict(self.iteritems())
cmtys.cmtysizes(): dict
    Mapping name->len(nodes)
cmtys.cmtysizes_sum() : int
    Sum of all community sizes (may be more than N if there are overlaps)
cmtys.nodecmtys(): dict
    Mapping node -> set(communities_of_node).  The reverse of
    iteritems() and to_dict()
cmtys.nodecmtys_onetoone(): dict
    Mapping node -> community.  Like above, but raises error if
    communities overlap.
cmtys.nodes: set
    Set of all nodes in this object.  If not explicitely given at
    creation time, this is union of all self.itervalues().  This
    differs from the below since Communities objects (the dict
    based one) can be given a nodes= option at creation, which
    represents the universe.  This is needed in order to support
    non-spanning systems.
cmtys.nodes_spanned(): set
    Union of all self.itervalues().  Compare to above.
cmtys.cmtyintmap(): dict
    Mapping community_names -> int.  This can be used as a mapping
    to ints.  May not be totally efficient in space or time (since
    you can't optimize both at the same time).
cmtys.nodeintmap(): dict
    Mapping community_names -> int.  This can be used as a mapping
    to ints.  May not be totally efficient in space or time.
fraction_covered(): float
    Fraction of network covered by at least one community.  Only makes
    sense for objects which have a node universe given at creation.
overlap(): float
    sum(all_community_sizes) / N
is_cover(): bool
is_non_overlapping(): bool
is_partition(): bool
    Exactly what they say.
Q(g): int
    Modularity calculation.  Requires argument of the networkx graph.


Input/output
~~~~~~~~~~~~

These constructors should work for any subclass that implements
from_iter():

from_iter(nodecmtys):  new instance
   Constructor from iterator over (nodes, communities).  This is the
   'base constructor' that all other subclasses should implement, if
   they want to be creatable from any of the below constructors.
from_nodecmtys(nodecmtys):  new instance
   Constructor from mapping node -> community
from_nodecmtys(nodecmtys):  new instance
   Constructor from mapping node -> set(communities)
from_pcd(G): new instance
   Constructor from pcd object G
from_networkx(g): new instance
   Constructor from networkx object with 'cmty' or 'cmtys' attributes
   on nodes.

to_dict(G): dict
   Convert to dictionary.
to_full(G): new Communities instance
   Return a Communities object that is based on a dictionary.  This
   is used to load any subclass into memory.
to_pcd(): pcd.Graph
to_pcd_cmtystate(): dict
load_pcd(pcd G):
load_networkx_custom(g):
    Add communities to networkx node data.
load_networkx_online(g):
    Add communities to networkx node data.
write_gml(g, fname):
to_networkx_label(g):
write_clusters(self, fname, ...):
    Write to a clusters file, consisting of one line per community.

Other
~~~~~

stats(): list
    Print out some stats on community sturcture
list_communities(): list
    List all communities
list_overlapping_nodes(): list

There is a lot more that isn't documented.


Other classes
=============

Communities:
    This is the main class which uses a dict as community storage.

CommunityFile:
    Iterative access of a one-line-per-community file.

CommunityFilter:
    Iterative filtering of communities.  Can apply transformations to
    communitiy names and nodes, such as largest component, splitting,
    ... .  No merging yet.

CommunityUnion:
    Union of several community objects, either allowing or not
    allowing duplicates.


Community interface
===================

Summary of the above copy and pasted from somewhere else.  The above
is much more complete.

comms.label

# What is the raw filename?
comms.fname

# iterate over comm_name, nodes:
comms.iteritems()

# iterate over only community names:
comms.iterkeys()

# iterate over only community names:
comms.itervalues()

# What if we want the full community object (not an iterator)?
# This loads everything into memory.  Once it is in memory, you
# can do comms[cname] to get the contents of community cname
# (random access, as opposed to iterative access).
comms.to_full()

# Convert this to a dictionary.  Generally, there is no need to do
# this instead of to_full
comms.to_dict()

# What if you want a mapping of node -> communities?
comms.nodecmtys()  # dict node -> set((c1, c2, ...))

"""

import collections
import numpy
import os
import re
import textwrap
import time

import networkx

import pcd.nxutil

class OverlapError(Exception):
    pass

# This is used to override the default open / os.path.exists to
# provide transparent compression and decompression.
import __builtin__
if 'open' not in globals():   open   = __builtin__.open
if 'exists' not in globals(): exists = os.path.exists

class _CommunitiesBase(object):
    """Base properties for community objects.

    This object contains methods which rely only on these assumptions:

    .iteritems()
    .iterkeys()
    .itervalues()
    __len__()
    __getitem__()
    copy()

    Methods you might want to implement, if you can do so efficiently:
    __repr__ (include additional cheap status information)
    .cmtynodes()
    .nodes()
    .cmtynames()

    """
    def __repr__(self):
        return '<%s object at %s>'%(self.__class__.__name__, hex(id(self)))
    # Some convenience functions about community structure.
    @property
    def N(self):
        """Number of nodes.

        Overlapping nodes are only counted once."""
        return len(self.nodes)
    @property
    def q(self):
        """Number of communities."""
        return len(self)
    def cmtynodes(self):
        return self.to_dict()
    def cmtynames(self):
        """Set of all community names."""
        return set(self.iterkeys())
    def cmtysizes(self):
        """Mapping of cmty -> len(cmty_nodes)"""
        return dict((c, len(ns)) for c,ns in self.iteritems())
    def cmtysizes_sum(self):
        """Total number of nodes in communities.  If there are
        overlaps, count the node multiple times."""
        return sum(len(v) for v in self.itervalues())
    def nodecmtys(self):
        """Return mapping node -> cmty set.

        Returns a dictionary {n0:set(c00,c01,...), n1:set(c10,c11), ...}
        """
        nodecmtys = { }
        for c, nodes in self.iteritems():
            for n in nodes:
                nodecmtys.setdefault(n, set())
                nodecmtys[n].add(c)
        return nodecmtys
    def nodecmtys_onetoone(self):
        """Return mapping node -> community if there are no overlaps.

        Unlike the `nodecmtys()` method, this returns a dict with
        values as strings or intergers (or whatever).  As such, this
        data structure can not handle overlap.  If there are overlaps,
        raise an exception.

        If there are nodes that are in no community, they will be
        missing from the return value.
        """
        nodecmtys = { }
        for c, nodes in self.iteritems():
            for n in nodes:
                if n in nodecmtys:
                    raise OverlapError("Overlapping: node %s in cmtys %s and %s"%
                                       (n, nodecmtys[n], c))
                nodecmtys[n] = c
        return nodecmtys
    def _has_nodes_universe(self):
        """Is there a explicit universe of nodes?

        If false, then .nodes is the same as .nodes_spanned(), which
        is the union of all communities.  If true, the user has
        specified some universe of nodes, so that the communities can
        possibly be non-covering."""
        return False
    @property
    def nodes(self):
        """Return all nodes spanned.

        This may be different from nodes_spanned, since this function
        may be overwritten to return a different 'universe' of nodes."""
        return self.nodes_spanned()
    def nodes_spanned(self):
        """Calculate the union of all communities.

        This should be the set of all nodes in the system, if the
        communities span the system.  If the system is not completly
        covered, then this is only the nodes within at least one
        community."""
        nodes = set()
        for ns in self.itervalues():
            nodes.update(ns)
        return nodes

    def cmtyintmap(self, consecutive=True, new=False):
        """Return a mapping from community names to integer indexes.

        Returns a dictionary mapping (cname) -> int in range(0,q-1).
        This might require iterating through all community names
        twice, unless new=True in which case it will only require
        iterating through communities once but always make a new
        dictionary.

        consecutive: bool, default true
           if false, do not require integers be consecutive.  This
           might make it easier to return the current  mapping
        new: bool, default false
            If true, always return a new mapping, which might not
            preserve existing names.
        """
        # Ensure everything is integers:
        if not new \
               and all(isinstance(c, int) for c in self.iterkeys()):
            q = len(self)
            # If we require consecutive.
            if not consecutive \
                   or all( 0 <= name < q  for name in self.iterkeys()):
                # TODO: return a shortcut to making a new dictionary
                return dict((c, c) for c in self.iterkeys())
        return dict((c, i) for i,c in enumerate(self.iterkeys()))
    def nodeintmap(self, consecutive=True, new=False):
        """Map from node names to integers.

        TODO: if nodes are already integers, do not construct a new mapping"""

        # Ensure everything is integers:
        if not new \
               and all(isinstance(n, int) for n in self.nodes):
            N = self.N
            # If we require consecutive.
            if not consecutive \
                   or all( 0 <= n < N  for n in self.nodes):
                # TODO: return a shortcut to making a new dictionary
                return dict((n, n) for n in self.nodes)

        return dict((n, i) for i,n in enumerate(self.nodes))

    def fraction_covered(self):
        """Fraction of networked covered by at least one community."""
        # if self._nodes is none, we have no total node universe, so
        # the universe is by definition the union of all communities,
        # so by definition everything is spanned.
        if not self._has_nodes_universe():
            return 1.0
        # Calculate actual span
        return len(self.nodes_spanned()) / float(self.N)
    def overlap(self):
        """Measure of degree of overlap.

        sum(cmty_sizes)/n_nodes.
        - A partition has overlap=1,
        - >1 means greater degree of overlapgreater degrees of overlap
        - <1 means that the system is not spanned (but the converse is not true)"""
        return sum(len(nodes) for nodes in self.itervalues())/float(self.N)

    def is_cover(self, nodes=None, strict=True):
        """Is every node in at least one community?

        nodes: if given, use this as the universe of nodes to be covered.
        strict: if true, raise an except if we span more than set of nodes.
        """
        if nodes is None:
            if not self._has_nodes_universe():
                # If we aren't given self._nodes, we always cover the
                # universe since the universe is defined as the union
                # of all communities.
                return 1.0
            nodes = self.nodes
        spannednodes = self.nodes_spanned()
        if strict:
            if len(spannednodes - nodes) > 0: # No spanned nodes left over.
                raise ValueError("We cover more than the universe of nodes, and strict mode is enabled.")
        return len(nodes - spannednodes) == 0 # Every node in spanned nodes.
    def is_non_overlapping(self):
        """Is no node is in more than one community?"""
        # Recalculate spannednodes here (instead of using
        # self.nodes_spanned()) since we need to add up the
        # sum of community sizes anyway.
        total_nodes = 0
        spannednodes = set()
        for ns in self.itervalues():
            spannednodes.update(ns)
            total_nodes += len(ns)
        return total_nodes == len(spannednodes)
    def is_partition(self, nodes=None):
        """Is every node in exactly one community?

        nodes: if given, use this as the universe of nodes to be
        covered.

        Should be equivalent to `is_non_overlapping() and is_cover()`."""
        # This should be identical to is_cover and non_overlapping

        # Recalculate spannednodes here (instead of using
        # self.nodes_spanned()) since we need to add up the
        # sum of community sizes anyway.
        spannednodes = set()
        total_nodes = 0
        for ns in self.itervalues():
            spannednodes.update(ns)
            total_nodes += len(ns)
        # is_non_overlapping:
        if not ( total_nodes == len(spannednodes) ):   # not (no node in multiple cmtys)
            return False
        # is_cover:
        # If neither nodes nor self._nodes is given, then we are not
        # given a universe, so we automatically span the system.
        if  nodes is not None  or  self._has_nodes_universe():
            if nodes is None:
                # This only happens if nodes is None and self._nodes is not.
                nodes = self.nodes
            if not (len(nodes - spannednodes) == 0): # not (every node covered)
                return False
        return True
    def _Q_pcd(self, g, gamma=1.0):
        """Network modularity, computed using pcd C functions"""
        if gamma != 1.0: raise NotImplementedError('gamma != 1.0')
        G = pcd.Graph.fromNetworkX(g)
        self.load_pcd(G)
        return G.modularity()
    def _Q_pcd_energy(self, g, gamma=1.0):
        """Network modularity, computed using pcd C functions"""
        if gamma != 1.0: raise NotImplementedError('gamma != 1.0')
        G = pcd.Graph.fromNetworkX(g)
        G.enableModularity(None)
        self.load_pcd(G)
        energy = G.energy(gamma)
        #print sum(G.rmatrix.diagonal())/g.number_of_edges()
        print "Q_with_diag=", ( -energy - sum(G.rmatrix.diagonal()/2.) ) / g.number_of_edges()
        print "Q=", -energy / float(g.number_of_edges())
        return -energy / float(g.number_of_edges())

    def _Q_cmty(self, g, gamma=1.0):
        """Modularity, computed in a community-centric method.

        This turns out to be much slowen than Q, due to the cost of
        repeatedly constructing the numpy matrices.  This will
        probably be more efficient when the number of communities is
        smaller, or when constructing the full modularity matrix is
        impractical."""
        def to_numpy_matrix(g, nodelist, dtype):
            """More efficient to_numpy_matrix.

            The networkx function of this uses the entire
            G.adjacency_iter() function, which is far too inefficient
            for use with every single community.  This is based off
            that function in networkx/convert.py, but far more
            efficient."""
            nlen = len(nodelist)
            assert nlen == len(set(nodelist)), "Duplicate nodes in community."
            assert not g.is_multigraph(), "This not valid for multigraphs (yet)"
            assert not g.is_directed(), "This not valid for multigraphs (yet)"
            index = dict(zip(nodelist, range(nlen)))
            M = numpy.zeros((nlen,nlen), dtype=dtype)
            for u, v, data in g.edges_iter(nodelist, data=True):
                #if v not in index: continue
                try:
                    # This is where we assume the graph is undirected.
                    # Duplicate this for loop for directed graphs, and
                    # place edges only once.
                    M[index[u],index[v]] = M[index[v],index[u]] = data.get('weight', 1)
                except KeyError:
                    pass
            return M

        mod = 0
        cmtynodes = self.cmtynodes()
        E = g.number_of_edges()
        # Old method creating numpy matrices for everything.
        #for c in self.cmtynames():
        #    nodeList = list(cmtynodes[c])
        #    #print c, len(nodeList)
        #    in_degrees = [ g.degree(n) for n in nodeList ]
        #    #out_degrees= [ g.degree(n) for n in nodeList ]
        #    out_degrees = in_degrees
        #
        #    expected_degree = numpy.multiply.outer(in_degrees, out_degrees) / float(2*E)
        #    #print expected_degree
        #    #adj_matrix = networkx.to_numpy_matrix(g, nodelist=nodeList, dtype=float)
        #    adj_matrix = to_numpy_matrix(g, nodelist=nodeList, dtype=float)
        #    #print adj_matrix
        #    mod_cmty = numpy.sum((adj_matrix - expected_degree)) / float(2*E)
        #    #print 'ma', numpy.sum(adj_matrix)/2., numpy.sum(expected_degree), mod_cmty
        #    mod += mod_cmty
        # New method computing only what is needed.
        for c in self.cmtynames():
            nodeList = list(cmtynodes[c])
            #print c, len(nodeList)
            in_degrees = [ g.degree(n) for n in nodeList ]
            #out_degrees= [ g.degree(n) for n in nodeList ]
            out_degrees = in_degrees

            tot_expected_degree = sum(in_degrees)*sum(out_degrees) / (float(2*E)**2)
            number_of_edges = g.subgraph(cmtynodes[c]).number_of_edges()
            mod_cmty = (number_of_edges / float(E)) - gamma*tot_expected_degree
            #print 'mb', number_of_edges, tot_expected_degree, mod_cmty
            mod += mod_cmty
        return mod


    # I/O
    #
    # These functions create a Communities object out of other objects, or
    # save the graph structure elsewhere.
    @classmethod
    def from_nodecmtys(cls, nodecmtys, **kwargs):
        """Create new Communities object from mapping nodes->cmty.

        This object takes a dictionary mapping from nodes to
        communities, makes the cmty->node dictionary, and returns an
        object."""
        cmtynodes = { }
        for n, c in nodecmtys.iteritems():
            cmtynodes.setdefault(c, set()).add(n)
        return cls.from_dict(cmtynodes=cmtynodes, **kwargs)
    @classmethod
    def from_nodecmtys_overlap(cls, nodecmtys, nodes=None, **kwargs):
        """Creat new Communities object from mapping nodes->set(cmtys).

        Differs from `from_nodecmtys` in that this assumes the mapping
        values are an iterable which contain all communities of the
        nodes.

        nodes: if given, this is the complete set of all nodes.  If
        not given, set complete set of all nodes from the nodecmtys
        keys."""
        cmtynodes = { }
        for n, cmtys in nodecmtys.iteritems():
            for c in cmtys:
                cmtynodes.setdefault(c, set()).add(n)
        if nodes is None:
            nodes = set(nodecmtys)
        return cls.from_dict(cmtynodes, nodes=nodes, **kwargs)
    @classmethod
    def from_clustersfile(cls, fname, converter=str, nodes=None):
        """Convert a clusters file into community structure.

        All of the hard work is done in the CommunitiesIterator class.
        That class is called, and the cmtynodes dictionary is
        constructed from it, and a new Communities object is created."""
        c = CommunityListIterator(fname, converter=converter)
        return cls.from_iter(c.iteritems(), nodes=nodes)
    @classmethod
    def from_pcd(cls, G):
        """Convert a pcd.Graph into Communities object.

        Communities are loaded from pcd.Graph, using G._nodeLabel."""
        d = G.cmtyDict()
        cmtynodes = { }
        if hasattr(G, '_nodeLabel'):
            for c, nodes in d.iteritems():
                cmtynodes[c] = set(G._nodeLabel[n] for n in nodes)
            nodes = set(G._nodeIndex.iterkeys())
            assert len(G._nodeIndex) == len(G._nodeLabel)
        else:
            for c, nodes in d.iteritems():
                cmtynodes[c] = set(nodes)
            nodes = set(range(G.N))
        return cls.from_dict(cmtynodes=cmtynodes, nodes=nodes)
    @classmethod
    def from_networkx(cls, g):
        """Create a new Communities object from a networkx.Graph object.

        The communities are stored in each node data dict (g.node[n])
        in one of these keys (the first takes priority):
          g.node[n]['cmty'] = c
          g.node[n]['cmtys'] = set((c1,c2,c3,...))
        """
        nodes = set(g.nodes_iter())
        cmtynodes = { }
        for node, d in g.nodes_iter(data=True):
            for c in pcd.nxutil._iterCmtys(d):
                if c not in cmtynodes: cmtynodes[c] = set()
                cmtynodes[c].add(node)
        return cls.from_dict(cmtynodes=cmtynodes, nodes=nodes)

    #
    # Converting to other formats
    #
    def to_dict(self):
        """Return dictionary of cmty->node_set"""
        return dict(self.iteritems())
    def to_full(self):
        return Communities(self.to_dict())
    def to_pcd(self, g=None, sparse=True):
        """Load this community structure into a pcd.Graph structure.

        Returns pcd.Graph with this community structure stored within
        it."""
        if g is None:
            G = pcd.Graph(N=self.N, sparse=sparse)
            G._makeNodeMap(self.nodes)
        else:
            G = pcd.Graph.fromNetworkX(g, sparse=sparse)
        G = self.load_pcd(G, clear=True)
        return G
    def to_pcd_cmtystate(self, G=None):
        """Return the community structure as pcd.Graph state.

        This is suitable for usage in pcd.Graph.setcmtystate().  Thin
        wrapper around self.to_pcd().getcmtystate() (but could someday
        be optimized beyond that)."""
        if G is not None:
            G = G.copy()
            self.load_pcd(G)
        else:
            G = self.to_pcd()
        return G.getcmtystate()

    def load_pcd(self, G, clear=True):
        """Loads communities from self into a pcd graph G.

        clear: if True, zero all present communities before loading.
        Otherwise, simply add all nodes to the new communities."""
        cmtynodes = self.cmtynodes()
        non_overlapping = self.is_non_overlapping()
        # Community map: mapping communities to integers, if not already done.
        assert len(cmtynodes) <= G.N
        if (    all(isinstance(c, int) for c in cmtynodes.iterkeys())
            and max(cmtynodes.iterkeys()) < G.N ):
            cmty_map = dict((c, c) for c in cmtynodes.iterkeys())
        else:
            # make a community map, community names -> integer indexes.
            cmty_map = self.cmtyintmap()

        if clear:
            G.cmtyListClear()


        G.oneToOne = 0
        if non_overlapping and self.is_cover():
            G.oneToOne = 1
        for c, nodes in cmtynodes.iteritems():
            cid = cmty_map[c]
            for n in nodes:
                if not clear and G.cmtyContains(cid, G._nodeIndex[n]):
                    # Ignore nodes which are already in the community
                    # (raises an error in pcd code)
                    continue
                G.cmtyListAddOverlap(cid, G._nodeIndex[n])
                if non_overlapping:
                    G.cmty[G._nodeIndex[n]] = cid
        return G
    def load_networkx_custom(self, g, attrnameset='cmtys', attrname=None, type_=None,
                             clear=True):
        """Load the communities from this object onto node attributes"""
        # Remove existing community labels
        if clear:
            for node, data in g.nodes_iter(data=True):
                data.discard(attrname)
                data.discard(attrnameset)
        nodecmtys = self.nodecmtys()
        for node, cmtys in nodecmtys.iteritems():
            if attrname is not None and len(cmtys) == 1:
                # non-overlapping, attrname is given
                if type_ is None:
                    g.node[node][attrname] = next(iter(cmtys))
                else:
                    g.node[node][attrname] = type_(next(iter(cmtys)))
            elif attrnameset is None and len(cmtys) > 1:
                # overlapping, attrnameset not given
                if type_ is None:
                    g.node[node][attrname] = cmtys
                elif type_ is str:
                    g.node[node][attrname] = ' '.join(str(c) for c in cmtys)
                else:
                    g.node[node][attrname] = type_(cmtys)
            else:
                # overlapping, attrnameset is given.
                # OR: non-overlapping and attrname is not given.
                if type_ is None:
                    g.node[node][attrnameset] = cmtys
                elif type_ is str:
                    g.node[node][attrnameset] = ' '.join(str(c) for c in cmtys)
                else:
                    g.node[node][attrnameset] = type_(cmtys)
    #def load_networkx(self, g):
    #    self.load_networkx_custom(g, attrname=None, type_=None)
    def load_networkx_online(self, g, attrnameset='cmtys', type_=None, clear=True):
        """Load these communities into a networkx graph as node attributes.

        This function returns None.  The graph g is modified so that
        each node has a set which lists all communities that node is
        within.

        The node community sets are written into the argument
        'attrnameset' option, by default 'cmtys'.  If clear is
        specified (default True), then all existing 'cmtys' attributes
        are removed first."""
        #from fitz import interact ; interact.interact('x')
        if clear:
            for node, data in g.nodes_iter(data=True):
                #data.pop(attrnameset, None)
                data[attrnameset] = set()
        nodes = g.node
        for cname, cnodes in self.iteritems():
            #print cname, cnodes
            for node in cnodes:
                #g.node[node].setdefault(attrnameset, set()).add(cname)
                g.node[node][attrnameset].add(cname)
    load_networkx = load_networkx_online

    def write_gml(self, g, fname):
        """Write a gml file, including labels."""
        #raise NotImplementedError
        g = g.copy()
        #nodecolors = self.nodecolors()
        #for n, color in nodecolors.iteritems():
        nodecmtys = self.nodecmtys_onetoone()
        for n, c in nodecmtys.iteritems():
            g.node[n]['label'] = str(c)
        networkx.write_gml(g, fname)
    def to_networkx_label(self, g, attrname='label'):
        """Modify a networkx graph, adding community IDs as node attributes.

        This method assumes and requires a one-to-one mapping node ->
        community.  If some nodes do not have communities, these are
        simply ignored.

        attrname: set this attribute name to be the community ID.
        Default: 'label'."""
        nodecmtys = self.nodecmtys_onetoone()
        for n, c in nodecmtys.iteritems():
            g.node[n]['label'] = c
    def write_clusters(self, fname, headers=[], mapping=None, raw=False,
                       write_names='separate'):
        """Write clusters to one-line-per-community file.

        This function does not write empty communities.

        fname: str, filename to write to (or file object).

        headers: list of headers to write at the top of the file.  A
        comment character and space '# ' is prepended.

        mapping: if given, node names are transformed using this
        mapping when writing.  Could be used to translate to integers,
        for example.  If given and 'int' or `int` (the type), then
        automatically greate a mapping to integers using
        self.nodeintmap().

        raw: if True, print no header lines or anything.

        write_names: 'separate'(default) or 'inline' or None.
           'separate': write another file with extension .names which has one line per community
        """
        # remove all empty communities:
        cmtynodes = self.cmtynodes()
        cnames = [cname for cname, cnodes in cmtynodes.iteritems() if len(cnodes) != 0 ]
        #cnames = self.cmtynames()

        f = fname
        if not hasattr(f, 'write'):
            f = open(f, 'w')
        # This writes out a separate file containing the true names
        # for the comunities.
        if write_names == 'separate':
            # Write the community names to disk:
            f_names = open(fname+'.names', 'w')
            print >> f_names, '# Community names for file %s'%fname
            print >> f, '#', repr(self)
            if isinstance(headers, str):
                print >> f_names, '#', headers
            else:
                for line in headers:
                    print >> f_names, '#', line
            # Write the label for self, if it exists:
            if hasattr(self, 'label'):
                print >> f_names, '# label:', self.label
            print >> f_names, '#', time.ctime()
            for cname in cnames:
                print >> f_names, cname
            f_names.close()
        if not raw:
            # Write representation of self: includes q and N.
            print >> f, '# Comunities, one community per line'
            print >> f, '#', repr(self)
            # Write any user-supplied headers.
            if isinstance(headers, str):
                print >> f, '#', headers
            else:
                for line in headers:
                    print >> f, '#', line
            # Write the label for self, if it exists:
            if hasattr(self, 'label'):
                print >> f, '# label:', self.label
            # Write the community names inline, if requested:
            if write_names == 'inline':
                print >> f, '# community names:', ' '.join(str(cname) for cname in cnames)
            # Write the current time.
            print >> f, '#', time.ctime()

        # if mapping='int', then force an integer mapping
        if (mapping == 'int' or mapping == int or raw) and mapping is None:
            mapping = self.nodeintmap()
            assert len(mapping) == self.N

        # Write all the communities.
        for c in cnames:
            if mapping:
                print >> f, ' '.join(str(x) for x in sorted(mapping[n] for n in cmtynodes[c]))
            else:
                print >> f, ' '.join(str(x) for x in sorted(cmtynodes[c]))



    #
    # Human-readable statistics on the structure.
    #
    def stats(self, width=120):
        """Human-readable statistics on the overall community sturucture."""
        cmtynodes = self.cmtynodes()
        cmtynodes = dict((k,v) for (k,v) in self.iteritems() if v)
        cmtynodes_nonsingle = dict((k,v) for (k,v) in self.iteritems()
                                   if len(v)>1)
        def iter_nonsingle():
            return ((k,v) for (k,v) in self.iteritems()
                    if len(v)>1)

        stats = [ ]
        stats.append("number of nodes: %d"%self.N)
        if hasattr(self, 'g'):
            stats.append("number of edges: %d"%self.g.number_of_edges())
        # Regular stats
        stats.append("number of communities: %d"%self.q)
        stats.append("mean community size: %f"%numpy.mean(
            [len(c) for c in self.itervalues()]))
        stats.append("std community size: %f"%numpy.std(
            [len(c) for c in self.itervalues()]))
        # Stats excluding singletons
        stats.append("number of non-singleton communities: %d"%(
            len(cmtynodes_nonsingle)))
        stats.append("mean community size (no singletons): %f"%numpy.mean(
            [len(c) for c in cmtynodes_nonsingle.values()]))
        stats.append("std community size (no singletons): %f"%numpy.std(
            [len(c) for c in cmtynodes_nonsingle.values()]))
        cmtys_by_rev_size = [c for (c,ns) in sorted(cmtynodes.iteritems(),
                                                    reverse=True,
                                                    key=lambda (c,ns):len(ns))]
        stats.append(textwrap.fill(("community sizes: %s"%' '.join(
                         str(len(cmtynodes[c])) for c in cmtys_by_rev_size)),
                     width=width,
                     initial_indent="", subsequent_indent="     ",
                     ))
        stats.append("fraction of nodes in largest community: %f"%(
            max(len(c) for c in cmtynodes.values())/float(self.N)))
        nodecmtys = self.nodecmtys()
        stats.append("fraction of network covered: %f"%self.fraction_covered())
        stats.append("fraction of network covered (non-singleton): %f"%(
            sum(1 for (n, cs) in nodecmtys.iteritems()
                if any((len(cmtynodes[c])>1) for c in cs))
            /float(self.N)))
        return stats
    def list_communities(self, width=120):
        """Human-readable list of community contents."""
        stats = [ ]
        stats.append("List of communities:")
        singleton_cmtys = [ ]
        cmtynodes = self.cmtynodes()
        cmtys_by_rev_size = [c for (c,ns) in sorted(cmtynodes.iteritems(),
                                                    reverse=True,
                                                    key=lambda (c,ns):len(ns))]
        for c in cmtys_by_rev_size:
            cmtysize = len(cmtynodes[c])
            if cmtysize == 1:
                singleton_cmtys.append(c)
                continue
            stats.append(textwrap.fill(
                (("%3s: "%c)+' '.join(str(n) for n in sorted(cmtynodes[c]))),
                width=width,
                initial_indent="", subsequent_indent="     ",
                ))
        stats.append(textwrap.fill(
            "Nodes in singleton communities: "+' '.join(
                str(next(iter(cmtynodes[c]))) for c in singleton_cmtys),
            width=width,
            initial_indent="", subsequent_indent="     ",
            ))
        return stats
    def list_overlapping_nodes(self, width=120):
        """Generate list of overlapping nodes."""
        stats = [ ]
        stats.append("How many communities does each node belong to?:")
        cmty_counts = collections.defaultdict(set)
        nodecmtys = self.nodecmtys()
        for node, cmtys in nodecmtys.iteritems():
            cmty_counts[len(cmtys)].add(node)
        for num_cmtys, nodes in cmty_counts.iteritems():
            stats.append(textwrap.fill(
                (("%4s: "%num_cmtys)+' '.join(str(n) for n in sorted(nodes))),
                width=width,
                initial_indent="", subsequent_indent="      ",
                ))
        return stats
    def nodeInfo(self, n, g):
        """Print detailed information about a node, including what
        communities it is connected to.

        Uses edges (not neighbors).
        """
        adjacency = [ ]
        adjacency.extend(neighbor for n,neighbor in g.edges(n))
        nodecmtys = self.nodecmtys()

        print "Node %s in communities %s"%(n, nodecmtys[n])
        #print "Edges connect to %s"%sorted(adjacency)
        cmtys = collections.defaultdict(list)
        #print "Connected communities for each adjacent node:"
        for n2 in sorted(adjacency):
            #print "  %10s: cmtys %s"%(n2, sorted(nodecmtys[n2]))
            for c in nodecmtys[n2]:
                cmtys[c].append(n2)
        print "Connections to each community:"
        for c, nodes in sorted(cmtys.iteritems()):
            print " %3s: %3d edges, %s"%(c, len(nodes), sorted(nodes))
    #
    # Operations on communities.
    #
    def nodecolors(self, colormap_name='gist_rainbow'):
        """Return a map of nodes->colors, based on communities.  Suitable for matplotlib."""
        #import matplotlib.cm as cm
        #import matplotlib.colors as mcolors
        #
        #colormap = cm.get_cmap(colormap_name)
        #normmap = mcolors.Normalize(vmin=0, vmax=self.q)
        #cmty_index = dict((c,i) for i,c in enumerate(self.cmtynames()))

        cmtycolors = self.cmtycolors(colormap_name=colormap_name)

        nodecolors = { }
        for node, c in self.nodecmtys_onetoone().iteritems():
            nodecolors[node] = numpy.asarray(cmtycolors[c])
        return nodecolors
    def cmtycolors(self, colormap_name='gist_rainbow'):
        """Return a map of communities -> colors.

        Create a mapping of communities to colors.  Every community is
        a different color, and there is no graph coloring employed.
        The result is suitable for matplotlib.

        Warning: this function is not finalized yet, the 'index'
        argument needs better handling."""
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors

        colormap = cm.get_cmap(colormap_name)
        normmap = mcolors.Normalize(vmin=0, vmax=self.q)
        cmty_index = dict((c,i) for i,c in enumerate(self.cmtynames()))

        cmtycolors = dict((c, colormap(normmap(cmty_index[c])))
                          for c in self.cmtynames())
        return cmtycolors

    #
    # Community-detection related statistics.  Comparing one cover with
    # another.
    #
    def cmty_mapping(self, otherCmtys, mode="overlap"):
        """For each community in self, find the community in g1 which most
        overlaps it.  Return a mapping g0 communities -> g1 communities.

        g0 is planted, g1 is detected.

        cmty_mapping[c0] = the c1 which contains the most intersected nodes

        Or: every community it g0 -> every best commmunity in g1.

        Comparison modes:
        overlap: which c1 has the most overlapping nodes with c0
        F1: which c0/c1 has the greatest F-score overlap.
        """
        cmtys0 = self.cmtynodes()
        cmtys0list = list(cmtys0.iterkeys())
        cmtys1 = otherCmtys.cmtynodes()
        cmtys1list = list(cmtys1.iterkeys())
        cmty0_to_cmty1 = dict()
        assigned_cmtys = set()
        plantedToDetectedMap = dict()

        overlaps = numpy.zeros(shape=(len(cmtys0), len(cmtys1)), dtype=int)
        #overlaps[detected, planted]

        for i0, c0 in enumerate(cmtys0list):
            bestScore = 0
            bestcmty = None
            c0nodes = cmtys0[c0]
            for i1, c1 in enumerate(cmtys1list):
                c1nodes = cmtys1[c1]
                overlap = len(c0nodes & c1nodes)
                overlaps[i0,i1] = overlap
                score = 0
                if mode == "overlap_forward":
                    # Using overlap / planted community size (note:
                    # constant divisor, so the division doesn't do
                    # anything really.)
                    if c1 in assigned_cmtys:
                        continue
                    if len(c0nodes) > 0:
                        score = float(overlap)/len(c0nodes)
                elif mode == "F1":
                    # Using f-score
                    precision = float(overlap)/len(c1nodes)
                    recall    = float(overlap)/len(c0nodes)
                    F1 = 2. * precision * recall / (precision + recall)
                    score = F1
                elif mode == "newman":
                    pass
                    #overlaps[c0,c1] = overlap
                else:
                    #raise ValueError("Unknown mode: %s"%mode)
                    pass #XXXXX

                if bestScore is None or score > bestScore:
                    bestScore = score
                    bestcmty = c1
            #cmty0_to_cmty1[c0] = bestcmty
            assigned_cmtys.add(c1)
        if mode == "overlap":
            cmty_map = { }
            sorted_overlaps = [ (overlaps[i,j], (i, j))
                                for i in xrange(overlaps.shape[0])
                                for j in xrange(overlaps.shape[1])
                                if overlaps[i,j] > 0]
            #if overlaps.shape[1] == 1: raise
            sorted_overlaps.sort(reverse=True, key=lambda x: x[0])
            planted_used = set()
            detected_used = set()
            for ov, (i,j) in sorted_overlaps:
                if i in planted_used or j in detected_used:
                    continue
                cmty_map[i] = j   # j is sometimes None XXXX
                if j is None: raise
                planted_used.add(i)
                detected_used.add(j)
            if None in cmtys0list or None in cmtys1list: raise
            for k, v in cmty_map.iteritems():
                cmty0_to_cmty1[cmtys0list[k]] = cmtys1list[v]
            if None in cmty0_to_cmty1.values(): raise
        elif mode == "newman":
            cmty0_to_cmty1 = { }
            # Newman algorithm in newman2004fast, footnote #19.
            plant_to_detect = collections.defaultdict(set)
            best_overlap_of_detected = overlaps.argmax(1)
            for i_planted, i_detected in enumerate(best_overlap_of_detected):
                c_planted = cmtys0list[i_planted]
                c_detected= cmtys1list[i_detected]
                plant_to_detect[c_detected].add(c_planted)
            for c_detect, c_planteds in plant_to_detect.iteritems():
                if len(c_planteds) == 1:
                    #print "mapping: c %s to c %s"%(c_planteds, c_detect)
                    cmty0_to_cmty1[c_planteds.pop()] = c_detect
                else:
                    #print "mapping: NOT done: %s all map to %s"%(
                    #    sorted(c_planteds), c_detect)
                    pass
            #from fitz import interactnow
        return cmty0_to_cmty1


    def frac_detected(self, detected, cmtyMapping, n_nodes=None, limit_nodes=None,
                      norm=False, full=False):
        """Compute fraction of nodes properly identified.

        self: planted configuration
        detected: detected communities, another Communities object.
        cmtyMapping: the mapping between detected -> planted communities.
            In order to calculate fraction detected, one must know some way
            to associate the detected communities with their corresponding planted
            communities.
        limit_nodes: if true, then should be a set which limits the number of nodes
            which are considered when calculating the fraction detected.  This can
            be used to implement a 'fraction of well-defined detected' measure or
            or some such.
        n_nodes: if specified, then assume that this is universe size.  Useful for
            partial mappings.  If false, require self.N == detected.N and use that
            value
        norm: if True, return ((frac-(1./q)) / (1-(1./q))).

        full: if True, then return tuple (fraction_detected, set(nodes_detected)).
            Otherwise just return fraction_detected.
        """
        if isinstance(cmtyMapping, str):
            cmtyMapping = self.cmty_mapping(cmtyMapping)
        if isinstance(cmtyMapping, tuple):
            cmtyMapping = self.cmty_mapping(cmtyMapping[0], **cmtyMapping[1])
        n_detected = 0
        n_total = 0
        cmtysPlanted = self.cmtynodes()
        cmtysDetected = detected.cmtynodes()
        allDetected = set()
        if n_nodes is None:
            assert self.N == detected.N
            n_nodes = self.N

        # for each planted communitiy select the best detected community
        for cPlanted, cDetected in cmtyMapping.iteritems():
            # How many nodes are planted?
            plantedNodes = cmtysPlanted[cPlanted]
            detectedNodes = cmtysDetected[cDetected]

            if limit_nodes is not None:
                plantedNodes = plantedNodes & limit_nodes  # explicit copy
            # how many nodes were correctly detected?  Add to correct number.
            #detected = cmtysPlanted[cPlanted] & cmtysDetected[cDetected]
            detected = plantedNodes & detectedNodes
            if limit_nodes is not None:
                detected.intersection_update(limit_nodes)
            # If we use the line below, we miscalculate the total
            # number of nodes we are trying to detect.
            n_total += len(plantedNodes)
            n_detected += len(detected)
            allDetected.update(detected)
            #print sorted(plantedNodes-detectedNodes)
            #print sorted(detectedNodes-plantedNodes)
            #print cPlanted, cDetected, len(plantedNodes), len(detected)

        # This is done because plantedNodes only includes planted
        # cmtys where there was any match.  We need to know the size of the entire universe.
        #if limit_nodes is not None:
        n_total = n_nodes
        if n_total == 0:
            #print limit_nodes
            if full:
                return float('nan'), allDetected
            else:
                return float('nan')
            #return 0
        frac_detected = float(n_detected) / n_total
        if norm:
            frac_detected = (frac_detected-(1./self.q)) / (1-(1./self.q))
        if full:
            return frac_detected, allDetected
        return frac_detected


    def illdefined(self, g,
                   mode='internaledgedensity',
                   insist_cover=True,
                   insist_non_overlapping=True):
        """Calculate the ill-defined nodes"""

        nodecmtys = self.nodecmtys()
        cmtynodes = self.cmtynodes()
        cmtysizes = self.cmtysizes()

        n_well_defined = 0
        n_ill_defined = 0
        illnodes = set()
        for n0, n0data in g.nodes_iter(data=True):
            #c0 = g.node[n0]['cmty']
            #print "zzz", n0, n0data
            #for c0 in _iterCmtys(n0data):
            if insist_cover and len(nodecmtys[n0]) < 1:
                raise ValueError("Node $r is in no communities"%n0)
            if insist_non_overlapping and len(nodecmtys[n0]) > 1:
                raise ValueError("Node $r is in multiple communities"%n0)
            if len(nodecmtys[n0]) > 1:
                raise NotImplementedError("We don't handle this case yet.")
            for c0 in nodecmtys[n0]:
                #print "zzz cmty", c0
                degree = 0
                degree_int = 0
                degree_ext = 0
                degree_ext_by_cmty = collections.defaultdict(int)
                #from fitz import interact ; interact.interact('b> ')
                # WARNING: fix for DiGraph or MultiGraph
                for n1 in g.neighbors_iter(n0):
                    n1cmtys = set(nodecmtys[n1])
                    #print "xx", n0, n1, n1cmtys
                    if insist_cover and len(n1cmtys) < 1:
                        raise ValueError("Node $r is in no communities"%n1)
                    if insist_non_overlapping and len(n1cmtys) > 1:
                        raise ValueError("Node $r is in multiple communities"
                                                                           %n1)
                    degree += 1
                    for c1 in n1cmtys:
                        if c0 == c1:
                            degree_int += 1
                        else:
                            degree_ext += 1
                            degree_ext_by_cmty[c1] += 1
                int_density = float(degree_int)/(cmtysizes[c0]-1)

                if len(degree_ext_by_cmty) == 0:
                    max_ext_degree = 0
                    max_ext_density = 0
                else:
                    max_ext_degree = max(degree_ext_by_cmty.itervalues())
                    density_ext_by_cmty = dict(
                        (_c, float(_deg)/cmtysizes[_c])
                        for _c,_deg in degree_ext_by_cmty.iteritems()
                        if _deg != 0
                        )
                    max_ext_density = max(density_ext_by_cmty.itervalues())


                #print 'y', int_density, max_ext_density

            # At what level of indention should this be?
            if mode == "degree":
                if degree_int > max_ext_degree:
                    n_well_defined += 1
                else:
                    n_ill_defined += 1
                    illnodes.add(n0)
            elif mode == "internaledgedensity":
                if int_density > max_ext_density:
                #if int_density >= max_ext_density:
                    #g.node[n0]['ill'] = False
                    n_well_defined += 1
                else:
                    #g.node[n0]['ill'] = True
                    n_ill_defined +=1
                    illnodes.add(n0)
            else:
                raise ValueError("Unrecognized mode: %s"%mode)
        assert n_well_defined + n_ill_defined == len(g)
        assert len(self.nodes) == len(g)
        assert n_ill_defined == len(illnodes)


        return dict(
            illnodes=illnodes,
            fracwell=float(n_well_defined)/len(g),
            )
    def illnodes(self, g, *args, **kwargs):
        """Return list of ill-defined nodes"""
        return self.illdefined(g, *args, **kwargs)['illnodes']


    #
    # These all implement partition-comparison measures.  They use
    # different techniques, so some are more efficient than others.
    #
    def VI(self, other):
        import pcd.util
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return pcd.util.VI(G1, G2)
    def ovIn(self, other):
        import pcd.util
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return pcd.util.N(G1, G2)
    def ovInw(self, other):
        import pcd.util
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return pcd.util.N(G1, G2, weighted=True)
    def ovIn_LF(self, other):
        """NMI using LF's code"""
        nmi = pcd.util.ovIn_LF(self, other)
        return nmi
    def In(self, other):
        import pcd.util
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return pcd.util.In(G1, G2)
    def F1(self, other, weighted=False):
        import pcd.F1
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return .5 * (  pcd.F1.F1(G1, G2)[0]
                     + pcd.F1.F1(G2, G1)[0])
    def F1w(self, other, weighted=False):
        """Like F1, but weighted by community size"""
        import pcd.F1
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return .5 * (  pcd.F1.F1(G1, G2, weighted=True)[0]
                     + pcd.F1.F1(G2, G1, weighted=True)[0])
    def Q(self, g, gamma=1.0):
        """Modularity computation using full numpy arrays."""
        if self.N > 1000:
            return self._Q_cmty(g, gamma=gamma)
        nodeList = g.nodes()
        nodecmty = self.nodecmtys_onetoone()
        cmtyList = tuple(nodecmty[n] for n in nodeList)

        cmtyintmap = self.cmtyintmap()
        cmtyList = tuple(cmtyintmap[c] for c in cmtyList)

        #nodeList = tuple(self._nodeLabel[i] for i in range(self.N))
        #cmtyList = tuple(self.cmty[i] for i in range(self.N))
        in_degrees = [ g.degree(n) for n in nodeList ]
        out_degrees= [ g.degree(n) for n in nodeList ]
        E = g.number_of_edges()

        expected_degree = numpy.multiply.outer(in_degrees, out_degrees) / float(2*E)
        #print expected_degree
        adj_matrix = networkx.to_numpy_matrix(g, nodelist=nodeList, dtype=float)
        #print adj_matrix
        same_cmty_mask = numpy.asarray(numpy.equal.outer(cmtyList, cmtyList), dtype=int)
        #print same_cmty_mask

        #print (adj_matrix - expected_degree) #* same_cmty_mask
        #print numpy.multiply((adj_matrix - expected_degree), same_cmty_mask)
        mod = numpy.sum(
            numpy.multiply((adj_matrix - gamma*expected_degree), same_cmty_mask)) \
            / float(2*E)

        #print mod, self._Q_cmty(g)
        return mod



class Communities(_CommunitiesBase):
    """Mapping-like communities storage.

    Store and operate on a community assignment.  The standard
    implementation of this object stores the entire community
    structure as a dictionary in the attribute self._cmtynodes.  When
    this dict is needed, call self.cmtynodes() to return a copy.  This
    entire object has a dictionary-like interface, with .iterkeys(),
    .itervalues(), and .iteritems() (and the default implementation is
    a thin mapping to self._cmtynodes).

    This particular object is not very efficient, and dynamical
    calculates most higher-order structures.  Consider a subclass if
    you need higher performance.

    This object does not store the full set of nodes spanned.  If the
    .nodes attribute is called, the union of all communities is
    calculated and returned.  In the case that this is NOT the correct
    set of all nodes on the graph (the community structure does not
    cover the graph), see the 'nodes' parameter to __init__.

    When a mapping nodes -> communities_containing is needed, call
    self.nodecmtys() and it will generate and return this.  This is
    not very efficient, so it is advisable to cache the result.

    This object, as implemented here, is not considered mutable.  You
    can, however, mutate the _cmtynodes object yourself and changes
    will be reflected since there is no other internal state.
    """
    _nodes = None
    def __init__(self, cmtynodes, nodes=None):
        """init

        cmtynodes: a dictionary mapping community_id -> set(cmty_nodes)

        nodes: if given, this is considered the full set of nodes in
        the graph.  If it is not given, then if .nodes is accessed,
        the object will calculate the union of all communities and return that.

        """
        if cmtynodes is None:
            raise ValueError('cmtynodes argument should not be None.')
        self._cmtynodes = cmtynodes

    def __repr__(self):
        return '<%s object with q=%d at %s>'%(self.__class__.__name__, self.q,
                                              hex(id(self)))

    # Mapping type emulation.
    def iterkeys(self):
        """Iterator over community names"""
        return self._cmtynodes.iterkeys()
    def itervalues(self):
        """Iterator over community contents (nodes)"""
        return self._cmtynodes.itervalues()
    def iteritems(self):
        """Iterator over community (names, nodes_within) pairs."""
        return self._cmtynodes.iteritems()
    def __getitem__(self, c):
        """Mapping emulation: return nodes within community c"""
        return self._cmtynodes[c]
    cmtycontents = __getitem__
    def __iter__(self):
        raise NotImplementedError("Use .iter{keys,values,items} instead")
    def __len__(self):
        """Number of communities"""
        return len(self._cmtynodes)
    # Other instance management
    def copy(self):
        """Copy of self, with a new cmtynodes dictionary.

        The community structure of the copy can be modified without
        the original's being affected, but all other mutable
        attributes are shared."""
        #return self.__class__(dict(self._cmtynodes))
        import copy
        new = copy.copy(self)
        new._cmtynodes = dict(self._cmtynodes)
        return new
    def to_full(self):
        """Return a full, dict-based copy of this object.

        This method is designed for subclasses to override.  When
        invoked, it should return a new community object with an
        explitit cmtynodes dictionary.  If called on a Communities
        object, return self.

        This does not return a dict, it returns a Communities object
        *based* on a dict."""
        return self

    @classmethod
    def from_dict(cls, cmtynodes, nodes=None):
        return cls(cmtynodes, nodes=nodes)
    @classmethod
    def from_iter(cls, nodes=None):
        return cls(dict(cmtynodes), nodes=nodes)

    # Return other data structures related to these communities.
    @property
    def nodes(self):
        """Set of all nodes spanned by this graph.

        If `nodes` was given when creating the object, use that set.
        If not, dynamically generate it from cmtynodes."""
        # The following 'if' line is for support of pickled objects
        # before this method existed.  They will have a 'nodes'
        # attribute saved on the self instance, so we should return
        # that directly.  Those don't have self._nodes (with
        # underscore).
        if 'nodes' in self.__dict__:
            return self.__dict__['nodes']
        elif getattr(self, '_nodes', None):
            return self._nodes
        else:
            return self.nodes_spanned()
    def _has_nodes_universe(self):
        if 'nodes' in self.__dict__: return True
        elif self._nodes: return True
        return False
    def cmtynames(self):
        """Set of all community names.

        Simply returns self._cmtynodes.keys()."""
        return set(self._cmtynodes.keys())
    def cmtynodes(self, copy=False):
        """Return mapping cmty -> node set.

        This is the core data structure, and how this object
        represents communities internally.  This method might be used
        when you want direct dict access, and not proxy everything
        through this method, for example within the body of an
        analysis script.  This method will simply return a reference
        to the self._cmtynodes object.  As such, do not modify the
        dictionary returned unless you pass copy=True which will make
        a copy.

        Returns a dictionary {c0:set(n00,n01,...), c1:set(n10,n11), ...}

        Warning: do not mutate the returned dictionary unless copy=True.

        For the Communities object, this is a thin mapping on top of
        self._cmtynodes.  Objects should only use this method if they
        need the complete dictionary.  Otherwise, they should use one
        of the self.iter* methods to produce an on-line algorithm, or
        self[c] to get the contents of one community.  These are more
        able to be implimented efficiently for other data structures.
        If this method is called on other objects, it will attempt to
        generate an entire dictionary, which will be large."""
        if not copy:
            return self._cmtynodes
        else:
            cmtynodes = { }
            for c, nodes in self._cmtynodes.iteritems():
                cmtynodes[c] = set(nodes)
            return cmtynodes
CmtyDict = Communities








class CommunityFile(_CommunitiesBase):
    """On-line community obejct, from a file.

    This object iterates through a file and returns communities.  The file is
    in the standard 'one line per community' format.

    This object implements a dictionary-like interface.  In
    particular, one can iterate through (cmty_name, cmty_nodes) pairs
    using the iteritems method, just as in dictionaries.  After the
    iteration is complete, self.q is set as the number of communities
    detected.  See `iteritems` method for more details.

    The number of communities attribute '.q' is None until one full
    iteration has completed.  Then, q is the true number of
    communities.

    """
    # These are all cached properties.
    _cmtynames = None
    _cmtynamesfile = None
    _label = None
    _nodes = None
    _N = None
    _q = None
    def __init__(self, fname, cmtynames=None, converter=str):
        """

        fname: input filename.

        converter: each node id in the file is passed through this
        function to convert it to a python object.  For example, to
        convert the nodes to integers, pass `int`.  Default: str."""
        self.fname = fname
        self.abspath_dir = os.getcwd()
        self.converter = converter
        if not exists(self.fname):
            raise ValueError("%s is not accessable"%self.fname)
        if cmtynames:
            if isinstance(cmtynames, str):
                if not exists(cmtynames):
                    raise ValueError("%s is not accessable"%cmtynames)
                self._cmtynamesfile = cmtynames
            else:
                self._cmtynames = cmtynames
        elif exists(self.fname+'.names'):
            self._cmtynamesfile = self.fname+'.names'
    def __repr__(self):
        """Repr of self: include q if it is known, otherwise don't."""
        if self._q is None:
            return '<%s(%s) object at %s>'%(self.__class__.__name__,
                                            self.fname, hex(id(self)))
        else:
            return '<%s(%s) object with q=%d at %s>'%(
                self.__class__.__name__, self.fname, self.q, hex(id(self)))
    def __len__(self):
        """Number of communities, or RuntimeError if not known yet."""
        if self._q is None:
            self._q = sum(1 for _ in self.iteritems())
        return self._q
    def __iter__(self):
        raise NotImplementedError("Use .iter{keys,values,items} instead")
    @property
    def N(self):
        """Calculate number of total nodes.  This runs the iterator fully."""
        if self._N is None:
            self._N = len(self.nodes)
        return self._N
    def _find_label(self):
        # At first, I was going to have the finding replace
        # self.__dict__['label'], but even when I do that it still
        # calls the descriptor from the class to find the attribute.
        # Since every lookup of .label calls this proxy, may as well
        # use the internal _label.
        if self._label is not None:
            return self._label
        # Search for label in the file
        data = open(self.fname).read(512)
        m = re.search(r'^# label: ([^\n]+)$', data, re.M|re.I)
        if m:
            label = self._label = m.group(1).strip()
            return label
        # Search for label in the names file.
        if self._cmtynamesfile is not None:
            data = open(self._cmtynamesfile).read(512)
            m = re.search(r'^# label: ([^\n]+)$', data, re.M|re.I)
            if m:
                label = self._label = m.group(1).strip()
                return label
        # If all else files, return the basename as the label.
        label = self._label = os.path.basename(self.fname)
        return label
    def _set_label(self, label):
        self._label = label
    label = property(fget=_find_label, fset=_set_label, doc="""Label of this community file""")
    def cmtynames(self):
        """Names of all communities.

        This is a function, because the Communities object has
        .cmtynames() as a function.

        Note: if no names can be found, return None, and user should
        use integer indexes."""
        # Cached copy, explicit list
        if self._cmtynames is not None:
            return self._cmtynames
        # Open from self._cmtynamesfile.  This is FILENAME.names if
        # not given explicitely and not explicitely given.
        if self._cmtynamesfile is not None:
            self._cmtynames = [
                x.strip() for x in
                   open(self._cmtynamesfile, 'rU').read().split('\n')
                if x.strip() and x[0]!='#' ]
            return self._cmtynames
        # No default, return None and user should just use integer
        # indexes.
        return None

    #def _set_cmtynames(self, cmtynames):
    #    print 'setting cmtynames'
    #    self.__dict__['_cmtynames'] = cmtynames
    #cmtynames = property(fget=_find_cmtynames, fset=_set_cmtynames,
    #                  doc="""Community names of this file.""")


    def cmtynodes(self):
        """Create the full dictionary of community structure."""
        return self.to_dict()
    def to_full(self, cls=Communities, nodes=None):
        """Convert to a full dict-based copy of this community structure.

        If nodes is given, use this as the node universe.  Be careful,
        if your file is large, this will use up much memory.
        """
        cmtys = cls(self.cmtynodes(), nodes=nodes)
        if self.label is not None:
            cmtys.label = self.label
        return cmtys
    def copy(self):
        """Copy of self.

        Returns shallow copy of self.  The actual community structure
        is immutable, thus there is no need to make a copy of it."""
        return copy.copy(self)

    def iterkeys(self):
        """Iterater over community names.

        This is simply a thin wrapper over self.iteritems()."""
        return (c for c,nodes in self.iteritems())
    def itervalues(self):
        """Iterate over community contents.

        This is simply a wrapper around self.iteritems()."""
        return (nodes for c,nodes in self.iteritems())

    def iteritems(self):
        """Iterate (cmty_name, cmty_nodes_list) pairs.

        This is the central iterator of the file.  Iterates through
        the file and returns (c, nodes) pairs, like
        dictionary.iteritems().  If we see a line that begins with
        exactly '# community names: ', this is taken as the true
        community namess, and these will be returned as community
        names `c`.

        After the iterator exits, this will save the total number of
        communities found as `self.q`"""
        number_of_cmty = 0
        converter = self.converter
        names = self.cmtynames()
        label = self.label

        file = open(os.path.join(self.abspath_dir, self.fname), 'rU')

        cmty_id = 0
        for lineno, line in enumerate(file):
            line = line.strip()
            if not line: continue
            if line[0] == '#':
                if not label and line.startswith('# label: '):
                    self.label = line[9:].strip() # maxsplit of 1
                if names is None and line.startswith('# community names: '):
                    names = self._cmtynames = line[19:].strip().split()
                # Comment line, do not parse it.
                continue

            # If we found labels, then use that as the name of the
            # community.  This is written out write_clusters.
            if names:   cname = names[cmty_id]
            else:       cname = cmty_id

            nodes = set(converter(x) for x in line.split())
            number_of_cmty += 1

            yield cname, nodes
            cmty_id += 1

        # Insert this into dict directly, since q is a property of a
        # superclass and doesn't support item assignment.
        self.__dict__['q'] = number_of_cmty
CommunityListIterator = CommunityFile

class CommunityFilter(_CommunitiesBase):
    def __init__(self, cmtys, filter):
        self._cmtys = cmtys
        self._filter = filter
    def __len__(self):
        return sum(1 for x in self.iteritems())
    def __repr__(self):
        return '<%s object with at %s>'%(self.__class__.__name__,
                                         hex(id(self)))
    def copy(self):
        return copy.copy(self)
    def cmtynodes(self):
        """Create the full dictionary of community structure."""
        return self.to_dict()
    def iterkeys(self):
        for cname, nodes in self.iteritems():
            yield cname
    def itervalues(self):
        for cname, nodes in self.iteritems():
            yield nodes
    def iteritems(self):
        filter = self._filter
        for cname, nodes in self._cmtys.iteritems():
            for new_name, new_nodes in filter(cname, nodes):
                yield new_name, new_nodes
#class CommunityProcess(Communities):
#    def __init__(self, cmtys, filter, remap=None):
#        self.cmtys = cmtys
#    def __repr__(self):
#        return '<%s object with q=%d at %s>'%(self.__class__.__name__, self.q,
#                                              hex(id(self)))
#    def iteritems(self):
#        filter = self.filter
#        for cname, nodes in self.cmtys.iteritems():
#            for new_name, new_nodes in filter(cname, nodes):
#                yield new_name, new_nodes


class CommunityUnion(_CommunitiesBase):
    """Combine multiple iterators into their union.

    This object takes multiple input Communities objects, and returns
    a new Communities object.  The new communities object iterates
    over each community object, and returns every community in every object.

    cmtys: the input Communities objects.  The interface should also
    work with dicts.

    dup_ok: if false, communities which exist in multiple levels are
    returned only one. (default false).  If true, the same community
    of nodes could be returned multiple times from this iterator."""
    def __init__(self, cmtys, dup_ok=False):
        self._cmtys = cmtys
        self._dup_ok = dup_ok
    def __len__(self):
        if self._dup_ok:
            return sum(len(c) for c in self._cmtys)
        else:
            return sum(1 for x in self.iteritems())
    def iterkeys(self):
        for cname, nodes in self.iteritems():
            yield cname
    def itervalues(self):
        for cname, nodes in self.iteritems():
            yield nodes
    def iteritems(self):
        dup_ok = self._dup_ok
        dups = set()
        for i, cmty in enumerate(self._cmtys):
            label = getattr(cmty, 'label', None)
            if label:
                label = "U-%02d-%s"%(i, label)
            else:
                label = "U-%02d"%i

            for cname, nodes in cmty.iteritems():
                cname = label+'-'+str(cname)
                # If duplicates are not allowed, do our duplicate checking.
                if not dup_ok:
                    nodes_hash = hash(frozenset(nodes))
                    if nodes_hash in dups:
                        continue
                    dups.add(nodes_hash)
                yield cname, nodes



def _test_interface(cmtys):
    repr(cmtys)
    #cmtys.__getitem__[]
    cmtys.nodes
    len(cmtys)
    cmtys.to_dict()
    cmtys.nodes
    cmtys.q
    cmtys.cmtynames()
    cmtys.nodes_spanned()
    cmtys.overlap()
    cmtys.is_cover()
    cmtys.is_non_overlapping()
    cmtys.is_partition()
    cmtys.fraction_covered()
    tuple(cmtys.iteritems())
    tuple(cmtys.iterkeys())
    tuple(cmtys.itervalues())

    assert set(cmtys.cmtyintmap().keys()) == set(cmtys.iterkeys())
    assert set(cmtys.cmtyintmap().values()) == set(range(len(cmtys)))
    assert set(cmtys.nodeintmap().keys()) == set(cmtys.nodes)
    assert set(cmtys.nodeintmap().values()) == set(range(cmtys.N))
    cmtys.cmtysizes()
    cmtys.cmtysizes_sum()




#def calc_sizes(cmtys):
#    cmtysizes = [ ]
#    for cname, nodes in cmtys.iteritems():
#        cmtysizes.append(len(nodes))
#
#
#cmtys1 = CommunityIterator('path/to/file.txt')
#cmtys2 = Communities(some_dict)
#
#calc_sizes(cmtys1)
#calc_sizes(cmtys2)
#
#
#
#class Graph(networkx.Graph):
#    def __init__(self, *args, **kwargs):
#        self.cmtys = Communities()

