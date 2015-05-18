.. PCD documentation master file, created by
   sphinx-quickstart on Thu May 14 11:09:15 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PCD's documentation!
===============================

Contents:

.. toctree::
   :maxdepth: 2

   cmty
   cmtycmp
   api/modules


What is PCD?
============

PCD is a python package for networks and community detection.  It was
started by Richard Darst in 2010 as a community detection code, but
since then has evolved to be a general network, community detection,
and research support library.  It does not aim to replace things like
:py:mod:`networkx`, instead it adds on various features missing in
these.  While there are some overall themes, it can mostly be
considered "any research code that has become useful in more than one
place".

Tour of available modules
=========================

:mod:`pcd.bulk`
    Apply community detection to many files automatically.

:mod:`pcd.cache`
    (unfinished) Tools to implement data caches to have code only
    re-run what does not already exist.

:mod:`pcd.cli`
    Command line interface to exploring communities.  Used by
    ``python -m pcd.cli``.  Not very useful at the moment.

:mod:`pcd.cmtycmp`
    Community partition comparison functions

:mod:`pcd.cmty`
    Community structure storage class and tools.

:mod:`pcd.commstats`
    Tools to automatically make plots of distribution of properties of
    communities.

``pcd/data/``
    Directory containing standard and test data.

:mod:`pcd.draw`
    Obsolete drawing tools.

:mod:`pcd.graphs`
    Various standard graphs, with community structure.

:mod:`pcd.grow`
    Several growing network models with triadic closure.

:mod:`pcd.igutil`
    Igraph interface utilities.  Minimal currently.

:mod:`pcd.ioutil`
    Input/output utilities.  In particular, has the :func:`pcd.ioutil.open_any`
    function that will open "any" graph file by examining the content
    to figure out what format it is in (not just by extension!).  It
    also has utilities for dealing with automatically compressed and
    uncompressed files, and writing ``pajek`` graph with community
    structure colored nicely.

:mod:`pcd.layout`
    Obsolete layout module.

:mod:`pcd.lazy`
    Utilities for lazy initialization of objects, including
    ``__XXX__`` method proxying and pickle support.

:mod:`pcd.nxutil`
    Tools to interface with and extend the ``networkx`` module.
    Contains :func:`pcd.nxutil.graphcolor` which heuristically color a
    graph in a way suitable for visualization.

``pcd/old/``
    Mostly C code from an older generation of ``pcd``.  No longer used.

:mod:`pcd.sbmdegree`
    Make degree-controlled SBM graphs (Karer and Newman).

:mod:`pcd.sbmgraphs`
    Tools for making classing Stochastic Block Model graphs.  This
    module is very old and crufty and should be rewritten.  Has
    detectability limit equations.

:mod:`pcd.sqlite`
    Utilities related to sqlite.  Right now just a on-disk
    dictionary-like key-value store :class:`pcd.sqlite.SQLiteDict`

:mod:`pcd.stats`
    Somewhat obsolete or useless.  Computes graph statistics all in
    one place.

``pcd/support/``
    A bunch of support code.  The distinction is that things in
    ``support`` tend to interface to other utilities instead of
    standing on their own, but the distinction is not entirely clear.

:mod:`pcd.support.algorithms`
    Interface to Community Detction algorithm code.  Allows you to
    call a single function on any :class:`networkx.Graph` and get
    :class:`pcd.cmty.Communities` out for further analysis, with no
    extra setup overhead and work.  This only interfaces with other
    code, not implements algorithms itself.

:mod:`pcd.support.cavity` (cython)
    Implementation of the cavity method

:mod:`pcd.support.fdop`
    Do weird things with UNIX file descriptiors to support better
    logging in :mod:`pcd.support.algorithms`.

:mod:`pcd.support.fixq`
    Fix up community partitions by fixing a certain number of
    communities.

:mod:`pcd.support.lfm`
    Hackish implementation of the local fitness optimization method of
    community detection.  Written as an example, *not* for any
    practical use.  Not even complete.

:mod:`pcd.support.lfr`
    Python interface to creating LFR graphs from the upstream LFR
    code.  Handles temporary files, directories, reading, converting
    to :class:`networkx.Graph` with community attributes.

:mod:`pcd.support.matplotlibutil`
    :mod:`matplotlib` utilities.  Mainly shortcut for generating
    figure and axes objects without using the object-oriented
    interface everytime.

:mod:`pcd.support.powerlaw`
    Generate numbers according to a power law flexibly (specify any
    two of min, max, mean).

:mod:`pcd.support.spectral`
    Spectral methods of community detection.

:mod:`pcd.support.testutil`
    Some unit test assertions related to graphs.

:mod:`pcd.tcmty`
    Temporal community storage class.

``pcd.test_*``
    Unit tests for files.

:mod:`pcd.tnetsql`
    :mod:`sqlite3` implementation of a temporal network storage class.

:mod:`pcd.util`
    Random utilities.  :class:`pcd.util.Averager` to compute
    numerically stable averages and other statistics on-line.
    :func:`pcd.util.chdir_context` is a context manager for changing
    directories.  :func:`pcd.util.tmpdir_context` is a
    context manager for temporary directories.
    :func:`pcd.util.writefile_context` handles writing files
    safely.  :class:`pcd.util.WeightedChoice` makes
    weighted :func:`random.choice`.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

