Community representation
========================


Design decisions
----------------

There are many ways of representing communities.  PCD tries to be as
general as possible, and not many any particular assumptions about the
community labels, node labels, and so on.  Also, all tools must
support cases such as arbitrary node labels, non-integer nodes,
overlapping nodes, communities which do not span an entire graph.
Community labels should also be arbitrary Python objects.
Modifications to the structure should be efficient.  For this, a
simple array that is a community list *not* work.

Instead, communities are stored as a mapping from community names to
sets of community nodes.  For example::

  c = {
      0: set([3, 5, 9, "A"]),
      2: set([2, 3, 4]),
      "some_name": set(["b", "c", "d"]),
      }

One may think "this is inefficient", but pcd was originally made for
research into community detection.  Flexibility is the primary
importance here.  Some features, like the arbitrary node label names,
just make work more enjoyable.

Data is not really stored in straight-up dictionaries, but an object
with this type of interface.  All functions using community structure
should primarily operate on this type.  There are conversion functions
to and from many other community structures for use in other contexts.


Basic definitions
-----------------
A node can be any hashable Python object (any object that can serve as
a dictionary key).  This corresponds to the ``networkx`` convention.
in particular, nodes in a partition do not have to just be integers.
They do not have to be consecutive in any order.  You can remove the
"first" node without having to reindex everything.

A community label can be any hashable Python object.

The nodes within a community are stored in a Python ``set``.  This
allows fast lookup and set operations such as intersect/union.

The basic interface looks like a dictionary, and is a mapping from
community name to set of community nodes (the ``cmtynodes``).
Experience showed that this was the most useful storage method for
most cases.  There are methods to easily transform one of these
objects into a mapping from nodes to communities communities
containing that node.

There are two levels of interface.  A dictionary interface lets you
immediately get the nodes of any community as a very fast operation.
``pcd`` is also designed to allow processing of community structures
too big to fit in memory.  This is done via an iterative interface.
For communities loaded like this (using :py:class:`pcd.cmty.CommunityFile`)

Community objects
-----------------

The primary object is :py:class:`pcd.cmty.Communities`.  There are
different types of community objects, but this one uses a straight
dictionary for storage.  It has these basic operations:

The following interface is designed to look like a dictionary.

- ``cmtys[cmty_label]`` returns the ``set`` of nodes corresponding to
  a community label.

- ``cmtys.iteritems()`` - iterate through ``(community_label,
  node_set)`` tuples.

- ``cmtys.iterkeys()`` - iterate through ``community_label``\ s.

- ``cmtys.itervalues()`` - iterate through ``node_set``\ s.


API documentation
-----------------

For full reference (and an older tutorial), please see the API
documentation for :py:mod:`pcd.cmty`.  This was written a while
ago, so the format is not updated.

You will find these objects:

- :py:class:`pcd.cmty.Communities`: dictionary-based community object.
- :py:class:`pcd.cmty._CommunitiesBase`: The base communities object.
  This defines most basic operations, using only the iterator
  interface.  dictionary-based community object.
- :py:class:`pcd.cmty.CommunityFile`: a communities object that can
  iterate through a file.  Communities do not have to be fully loaded
  in memory.
- :py:class:`pcd.cmty.CommunityFilter`: transform one community object
  into another by applying a filter.
- :py:class:`pcd.cmty.CmtyLimitNodes`: Transform one community
  structure into another, limiting to only a subset of communities.
- :py:class:`pcd.cmty.CmtyLimitAnyNodes`:
- :py:class:`pcd.cmty.CommunityUnion`:

