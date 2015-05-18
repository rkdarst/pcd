Community partition comparisons
===============================
The :py:mod:`pcd.cmtycmp` implements measures such as Normalized
Mutual Information or Variation of Information for comparing the
similarity of two partitions.  It depends on :mod:`pcd.cmty` for
handling the communities themselves.

Multiple implementations
------------------------
This module is structured to contain multiple implementations of each
comparison measure.  In the unit tests (:mod:`test_cmtycmp.py`),
implementations are compared to each other for accuracy purposes.  Also,
one could test their own implementation against external
implementations.

The following standard implementation groups exist:

- ``_python``: Naive (generally :math:`O(n^2)`) implementations.  These
  are easy to write but slower.
- ``_python2``: The fastest pure-Python implementations.
- ``_pcd``: Old C implementations from the previous generation of
  :mod:`pcd`.  These still exist in the ``old/`` directory but require.
  compiling C code.
- ``_igraph``: These transform data and make calls to the ``igraph``
  library.
- and more.

At the bottom of the file, there is a ``measures`` dictionary that
defines all measures that should do the same thing.  Below that, the
"default" implementation for each measure is defined.

Unit tests
----------
Measures should be unit tested.  These can be put in the
``test_cmtycmp.py`` module.  Run with ``nosetests test_cmtycmp.py``.

There are some default tests defined which all methods should pass.
They are:

``test_same``:
    M(P, P) = 1.  A measure should return a similarity of 1 when
    applied to the same partition twice.  This is tested on a unitary
    partition (all nodes in same community) and random partitions.
    (Not all measures satisfy this, these are excluded)
``test_one``:
    Different implementations of the same measure should return the
    same result when applied to the same partitions.  This tests using
    the unitary partition.
``test_random``:
    Same, but comparing the same random partition
``test_random_one``:
    Same, but comparing a random partition with the unitary partition.
``test_random_random``:
    Same, but comparing two random partitions.

Above are the basic, generalize tests, but others can and should be
made with specific partitions and known results for specific methods.

Other notes: Some methods (especially externally written methods) will
fail in random ways with certain data, so these have to be excluded.
There is some basic machinery to do this.  Also, this uses the
``nosetests`` framework of generator tests.

Adding a new measure
--------------------
1. Decide a unique name for your measure.  This should not be so
   generic that there will be confusion later once more measures are
   added.  It should also not be too long.  Let's say it is ``mycmp``.
2. Define some unit tests in ``test_cmtycmp.py``.  This should be
   known results.
3. (Optional, do this if it helps.)  Write a naive implementation
   ``mycmp_python``.  This should be easy to understand and modify,
   and can be slow.  For example, this can directly apply the
   mathematical formula even if it is :math:`O(n^2)`.
4. Add your test to the ``measures`` dictionary.  This should be added
   as a new key ``"mycmp": [mycmp_python]``.  Run unit tests until
   everything passes.
5. Write an optimized ``mycmp_python2`` implementation and add it to
   the ``measures`` dict.  Test it using your custom tests and also
   that it matches the ``_python`` implementation.

