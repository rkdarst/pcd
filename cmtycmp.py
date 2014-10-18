from collections import defaultdict
from functools import partial
from math import log
import subprocess

import pcd.cmty
# note: pcd.util imports names directly from this module, so be
# careful with circular imports.
import pcd.util

log2 = lambda x: log(x, 2)

def limit_to_overlap(cmtys1, cmtys2, nodes=None):
    if nodes is None:
        nodes = cmtys1.nodes & cmtys2.nodes
    cmtys1new = pcd.cmty.CmtyLimitNodes(cmtys1, nodes)
    cmtys2new = pcd.cmty.CmtyLimitNodes(cmtys2, nodes)
    return cmtys1new, cmtys2new


# Simple python-based implementations.  These are naive O(n^2)
# implementations, whose only purpose is for pedagogy and unit-testing
# comparison with the efficient python2 implementations.
#
def mutual_information_python(cmtys1, cmtys2):
    MI = 0.0
    N = len(cmtys1.nodes)
    assert N == len(cmtys2.nodes)
    for c1name, c1nodes in cmtys1.iteritems():
        for c2name, c2nodes in cmtys2.iteritems():
            n1 = len(c1nodes)
            n2 = len(c2nodes)
            n_shared = len(set(c1nodes) & set(c2nodes))
            if n_shared == 0:
                continue
            MI += (n_shared/float(N)) * log2(n_shared*N/float(n1*n2))
    return MI
def entropy_python(cmtys):
    H = 0.0
    N = float(len(cmtys.nodes))
    for cnodes in cmtys.itervalues():
        n = len(cnodes)
        if n == 0: continue
        H += n/N * log2(n/N)
    return -H
def vi_python(cmtys1, cmtys2):
    """Variation of Information"""
    I = mutual_information_python(cmtys1, cmtys2)
    VI = entropy_python(cmtys1) + entropy_python(cmtys2) - 2*I
    return VI
def nmi_python(cmtys1, cmtys2):
    I = mutual_information_python(cmtys1, cmtys2)
    Hs = entropy_python(cmtys1) + entropy_python(cmtys2)
    if Hs == 0 and I == 0:
        return 1.0
    if Hs == 0:
        return float('nan')
    In = 2.*I / Hs
    return In
In = nmi_python


# More efficient implementations using the confusion matrix.
def mutual_information_python2(cmtys1, cmtys2):
    N = len(cmtys1.nodes)
    assert N == len(cmtys2.nodes)
    # Compute cofusion matrix
    confusion = defaultdict(int)
    nodecmtys2 = cmtys2.nodecmtys_onetoone()
    for cname, cnodes in cmtys1.iteritems():
        for node in cnodes:
            confusion[(cname, nodecmtys2[node])] += 1
    # Compute actual MI
    MI = 0.0
    sizes1 = cmtys1.cmtysizes()
    sizes2 = cmtys2.cmtysizes()
    for (c1, c2), overlap in confusion.iteritems():
        MI += ((overlap / float(N))
               * log2(overlap*N / float(sizes1[c1]*sizes2[c2])))
    return MI
def vi_python2(cmtys1, cmtys2):
    """Variation of Information"""
    I = mutual_information_python2(cmtys1, cmtys2)
    VI = entropy_python(cmtys1) + entropy_python(cmtys2) - 2*I
    return VI
def nmi_python2(cmtys1, cmtys2):
    I = mutual_information_python2(cmtys1, cmtys2)
    Hs = entropy_python(cmtys1) + entropy_python(cmtys2)
    if Hs == 0 and I == 0:
        return 1.0
    if Hs == 0:
        return float('nan')
    In = 2.*I / Hs
    return In
def _nmi_LFK_python2_HXY_norm(cmtys1, cmtys2):
    """One half of the LFK NMI."""
    # Compute cofusion matrix
    confusion = defaultdict(lambda: defaultdict(int))
    nodecmtys2 = cmtys2.nodecmtys()
    for cname, cnodes in cmtys1.iteritems():
        for node in cnodes:
            for c2 in nodecmtys2[node]:
                confusion[cname][c2] += 1
    c1sizes = cmtys1.cmtysizes()
    c2sizes = cmtys2.cmtysizes()
    # Total community sizes
    N1 = float(len(cmtys1.nodes))
    N2 = float(len(cmtys2.nodes))
    N  = float(len(cmtys1.nodes | cmtys2.nodes))

    # These definitions correspond to the definitions in the LFK
    # paper, pretty much.
    def h(p):
        if p==0 or p==1: return 0
        return -p * log(p,2)
    def H(size, N):
        return h(size/N) + h((N-size)/N)
    def H2(c1size, c2size, overlap, N):
        hP11 = h((         overlap  )/N)
        hP10 = h((c1size - overlap  )/N)
        hP01 = h((c2size - overlap  )/N)
        hP00 = h((N -  c1size - c2size + overlap)/N)
        if hP11 + hP00 <= hP01 + hP10:
            # We want to exclude ones matching this condition.  Return
            # 'inf' and it will be excluded from the min() function
            # wrapping this one.
            return float('inf')
        hPY1 = h(  (  c2size) / N)
        hPY0 = h(  (N-c2size) / N)
        return hP11+hP00+hP01+hP10 - hPY1 - hPY0

    # The main NMI loop.  For every community in cmtys1, find the best
    # matching community (as seen by H2) and take the arithmetic mean.
    # This will later be upgraded to consider weighting by community
    # size.
    from pcd.util import Averager
    HX_Y = Averager()
    for c1name, c2name_overlaps in confusion.iteritems():
        c1size = c1sizes[c1name]
        # Best matching other community.
        HX_Yhere = min(H2(c1size, c2sizes[c2name], overlap, N)
                       for c2name, overlap in c2name_overlaps.iteritems())
        # Handle the case of no matches.
        if HX_Yhere == float('inf'):
            HX_Yhere = H(c1size, N1)
        # Norm and return.
        _ = H(c1size, N1)
        if _ == 0:
            HX_Y += 0
        else:
            HX_Y += HX_Yhere / _
    return HX_Y.mean
def nmi_LFK_python2(cmtys1, cmtys2):
    """NMI of LFK (overlapping)

    This implementation matches the results of the LF-distributed one.
    It is not yet maximally efficient, and also does not yet have the
    weighted options.
    """
    #print _nmi_LFK_python2_HXY_norm(cmtys1, cmtys2)
    #print _nmi_LFK_python2_HXY_norm(cmtys2, cmtys1)
    N = 1 - .5 * (_nmi_LFK_python2_HXY_norm(cmtys1, cmtys2)
                  + _nmi_LFK_python2_HXY_norm(cmtys2, cmtys1) )
    return N


def recl_python2(cmtys1, cmtys2, weighted=True):
    """Recall: how well cmtys2 is matched by cmtys1.

    Consider cmtys1 to be 'true' communities, and cmtys2 to be
    detected communities.  Recall measures how well the known
    communities are detected.  Extra detected communities do not
    affect this result."""
    confusion = defaultdict(lambda: defaultdict(int))
    nodecmtys2 = cmtys2.nodecmtys()
    for cname, cnodes in cmtys1.iteritems():
        for node in cnodes:
            for c2 in nodecmtys2.get(node, []):
                confusion[cname][c2] += 1
    sizes1 = cmtys1.cmtysizes()
    sizes2 = cmtys2.cmtysizes()

    recls = [ ]
    #sizes = [ ]
    for c1name, others in confusion.iteritems():
        size = sizes1[c1name]
        #if weighted:
        #    sizes.append(size)
        others = [ overlap/float(size + sizes2[c2name] - overlap)
                   for c2name, overlap in others.iteritems() ]
        if not others:
            recls.append(0.0)
        else:
            recls.append(max(others))
    if weighted == 2:
        return (sum(x*sizes1[cname]
                    for x, cname in zip(recls, confusion.iterkeys()))
               / sum(sizes1.itervalues()))
    elif weighted:
        if sum(recls) == 0:
            return 0.0
        sizes1weighted = { }
        nodecmtys1 = cmtys1.nodecmtys()
        for cname, cnodes in cmtys1.iteritems():
            sizes1weighted[cname] = sum(1./len(nodecmtys1[n]) for n in cnodes)
        return (sum(x*sizes1weighted[cname]
                    for x, cname in zip(recls, confusion.iterkeys()))
               / float(len(nodecmtys1)))
    else:
        return sum(recls) / float(len(recls))
def prec_python2(cmtys1, cmtys2, weighted=True):
    """Precision: hew well cmtys1 matches cmtys2.

    Consider cmtys1 to be 'true' communities, and cmtys2 to be
    detected communities.  Precision measures how well the detected
    communities all match the known.  Extra known communities do not
    affect this result."""
    return recl_python2(cmtys2, cmtys1, weighted=weighted)
def F1_python2(cmtys1, cmtys2, weighted=True):
    """The harmonic mean of recall and precision.

    This is a symmetric measure."""
    recl = recl_python2(cmtys1, cmtys2, weighted=weighted)
    prec = prec_python2(cmtys1, cmtys2, weighted=weighted)
    return (2.0 * prec * recl) / (prec + recl)
def recl_uw_python2(cmtys1, cmtys2):
    """Unweighted recall"""
    return recl_python2(cmtys1, cmtys2, weighted=False)
def prec_uw_python2(cmtys1, cmtys2):
    """Unweighted precision"""
    return prec_python2(cmtys1, cmtys2, weighted=False)
def F1_uw_python2(cmtys1, cmtys2):
    """Unweighted F1"""
    return F1_python2(cmtys1, cmtys2, weighted=False)


#
# Old implementation from my pcd C code
#
def _cmtys_to_pcdgraph(cmtys):
    return cmtys.to_pcd()
def mutual_information_pcd(cmtys1, cmtys2):
    from .old.util import mutual_information_c
    return mutual_information_c(
        _cmtys_to_pcdgraph(cmtys1), _cmtys_to_pcdgraph(cmtys2))
def nmi_pcd(cmtys1, cmtys2):
    from .old.util import In
    return In(
        _cmtys_to_pcdgraph(cmtys1), _cmtys_to_pcdgraph(cmtys2))
def vi_pcd(cmtys1, cmtys2):
    from .old.util import VI
    return VI(
        _cmtys_to_pcdgraph(cmtys1), _cmtys_to_pcdgraph(cmtys2))
def nmi_LFK_pcd(cmtys1, cmtys2):
    from .old.util import mutual_information_overlap
    return mutual_information_overlap(
        _cmtys_to_pcdgraph(cmtys1), _cmtys_to_pcdgraph(cmtys2))
def nmi_LFK_pcdpy(cmtys1, cmtys2):
    from .old.util import mutual_information_overlap_python
    return mutual_information_overlap_python(
        _cmtys_to_pcdgraph(cmtys1), _cmtys_to_pcdgraph(cmtys2))


#
# External code: Overlap NMI by L and F
#
def nmi_LFK_LF(cmtys1, cmtys2, check=True, use_existing=False):
    """Compute NMI using the overlap-including definition.

    This uses the external code 'mutual3/mutual' to calculate the
    overlap-using In.

    If the communities object is a pcd.cmty.CommunityFile, has a fname
    attribute, and it exists and doesn't end in .gz or .bz2, and the
    argument use_existing is true, then this file will be used for the
    input to the NMI code.  The NMI code is very dumb, and does not
    remove comments, and uses only floating-point or integer node IDs,
    and does not give any type of error if these conditions are not
    met.  Thus, only enable use_existing if you can vouch for the.

    check: bool, default True
        If true, check that all node names are either representable as
        integers or floating point numbers.  This can either be python
        integers or floats, or strings of those.
    use_existing: bool, default False
        If true, and community object has a '.fname' attribute, and
        that file exists, use that as the input filename instead of
        re-writing the file.  WARNING: if you use this, you must
        ensure that there are no comments within the file, or anything
        else which make cause the program to break.
        """
    from pcd.support.algorithms import _get_file
    from pcd.cmty import CommunityFile
    binary = _get_file('mutual3/mutual')

    def _is_float_str(x):
        """Is this a string representation of a float?"""
        try:
            float(x)
            return True
        except ValueError:
            return False
    # check that community files are valid for the program:
    if check:
        for nodes in cmtys1.itervalues():
            assert all(isinstance(x, (int, float)) or _is_float_str(x) for x in nodes )
        for nodes in cmtys2.itervalues():
            assert all(isinstance(x, (int, float)) or _is_float_str(x) for x in nodes )

    # We must use os.path.abspath *outside* of the tmpdir_context, or
    # else the absolute path will be wrong.
    args = [ binary ]
    if (use_existing
        and isinstance(cmtys1, CommunityFile)
        and hasattr(cmtys1, 'fname')
        and os.path.exists(cmtys1.fname)
        and not (cmtys1.fname.endswith('.bz2')
                 or cmtys1.fname.endswith('.gz'))):
        args.append(os.path.abspath(cmtys1.fname))
    else:
        args.append(None)

    if (use_existing
        and isinstance(cmtys1, CommunityFile)
        and hasattr(cmtys2, 'fname')
        and os.path.exists(cmtys2.fname)
        and not (cmtys2.fname.endswith('.bz2')
                 or cmtys2.fname.endswith('.gz'))):
        args.append(os.path.abspath(cmtys2.fname))
    else:
        args.append(None)

    with pcd.util.tmpdir_context(chdir=True, dir='.', prefix='tmp-nmi-'):
        # Write community files, if they do not already exist.  These
        # must be written inside of the tmpdir_context context because
        # only in here does it know the right
        if args[1] is None:
            cmtys1.write_clusters('cmtys1.txt', raw=True)
            args[1] = 'cmtys1.txt'
        #else:
        #    # Symlink it instead of run in-place (program can fail if
        #    # filenames are weird!)
        #    os.symlink(args[1], 'cmtys1.txt')
        #    args[1] = 'cmtys1.txt'
        if args[2] is None:
            cmtys2.write_clusters('cmtys2.txt', raw=True)
            args[2] = 'cmtys2.txt'
        #else:
        #    os.symlink(args[2], 'cmtys2.txt')
        #    args[1] = 'cmtys2.txt'
        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        ret = p.wait()
        stdout = p.stdout.read()
        if ret != 0:
            print stdout
            raise RuntimeError("The program '%s' returned non-zero: %s"%(
                args[0], ret))
        nmi = float(stdout.split(':', 1)[1])
    return nmi
# backwards compatability aliases
NMI_LF = nmi_LFK_LF
ovIn_LF = nmi_LFK_LF



#
# Igraph-based implementations
#
def to_membership_list(*cmtys_list):
    nodes = set.union(*(c.nodes for c in cmtys_list))
    results = [ ]
    for cmtys in cmtys_list:
        nodecmtys = cmtys.nodecmtys_onetoone()
        clist = [ nodecmtys[n] for n in nodes ]
        results.append(clist)

    return results, nodes
def _similarity_igraph(cmtys1, cmtys2, method):
    """Calculate community similarity using igraph

    Available methods:
    'vi' or 'meila'
    'nmi' or 'danon'
    'split-join'
    'rand'
    'adjusted_rand'

    Quirk/bug: Only the first character of these names is used!
    """
    import igraph
    (mlist1, mlist2), nodes = to_membership_list(cmtys1, cmtys2)
    val = igraph.compare_communities(mlist1, mlist2, method=method)
    if method == 'meila': val /= log(2)
    return val


def vi_igraph(cmtys1, cmtys2, **kwargs):
    return _similarity_igraph(cmtys1, cmtys2, method='meila', **kwargs)
def nmi_igraph(cmtys1, cmtys2, **kwargs):
    return _similarity_igraph(cmtys1, cmtys2, method='danon', **kwargs)
def rand_igraph(cmtys1, cmtys2, **kwargs):
    return _similarity_igraph(cmtys1, cmtys2, method='rand', **kwargs)
def adjusted_rand_igraph(cmtys1, cmtys2, **kwargs):
    return _similarity_igraph(cmtys1, cmtys2, method='adjusted_rand', **kwargs)


measures = {
    'vi': ['vi_python', 'vi_igraph', 'vi_pcd', 'vi_python2'],
    'mutual_information':
        ['mutual_information_python', 'mutual_information_pcd',
         'mutual_information_python2'],
    'nmi': ['nmi_python', 'nmi_igraph', 'nmi_pcd', 'nmi_python2'],
    'nmi_LFK': ['nmi_LFK_LF', #'nmi_LFK_pcd',
                'nmi_LFK_pcdpy', 'nmi_LFK_python2'],
    'rand': ['rand_igraph'],
    'adjusted_rand': ['adjusted_rand_igraph'],
    'F1':   ['F1_python2'],
    'recl': ['recl_python2'],
    'prec': ['prec_python2'],
    'F1_uw':   ['F1_uw_python2'],
    'recl_uw': ['recl_uw_python2'],
    'prec_uw': ['prec_uw_python2'],
    }

# Standard implementations

nmi = nmi_python
vi = vi_python
mutual_information = mutual_information_python
nmi_LFK = nmi_LFK_LF
rand = rand_igraph
adjusted_rand = adjusted_rand_igraph
F1 = F1_python2
recl = recl_python2
prec = prec_python2
