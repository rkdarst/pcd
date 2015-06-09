from collections import defaultdict
from functools import partial
from copy import deepcopy
import itertools
from math import log, sqrt
import subprocess
import numpy

import pcd.cmty
# note: pcd.util imports names directly from this module, so be
# careful with circular imports.
import pcd.util


# log2: from lambda to function (pep8)
def log2(x): return log(x, 2)

def comb2(n): return n*(n-1.0)/2.0

# equal in the sense of what (most) measures consider equal
def is_equal(cmtys1, cmtys2):
    if cmtys1.q == cmtys2.q and cmtys1.N == cmtys2.N:
        node_n1 = sorted(cmtys1.cmtysizes().values())
        node_n2 = sorted(cmtys2.cmtysizes().values())
        if node_n2 == node_n1:
            return True    
    return False
    
# equality check for measures that care about specific nodes 
# for example pair counting measures
def is_nodewise_equal(cmtys1, cmtys2):
    if cmtys1.q == cmtys2.q and cmtys1.N == cmtys2.N:
        nodesets = set(frozenset(i) for i in cmtys1.itervalues())
        for v in cmtys2.itervalues():
            nodesets.discard(set(v))
        if not nodesets:
            return True
    return False

def limit_to_overlap(cmtys1, cmtys2, nodes=None):
    if nodes is None:
        nodes = cmtys1.nodes & cmtys2.nodes
    cmtys1new = pcd.cmty.CmtyLimitNodes(cmtys1, nodes)
    cmtys2new = pcd.cmty.CmtyLimitNodes(cmtys2, nodes)
    return cmtys1new, cmtys2new

def _get_data(cmtys1, cmtys2):
    """Return important data between two communities. This is
    designed to be a general function for use within other comparison
    functions.

    Returns a tuple: (confusion, sizes1, sizes2).

    confusion: the confusion matrix
    sizes1:    dictionary cname->len(cmty) for cmtys1
    sizes1:    dictionary cname->len(cmty) for cmtys2.
    """
    confusion = { }
    sizes1 = { }
    sizes2 = { }
    # calculate nodecmtys2 here for greater efficiency, so we can
    # calculate sizes2 at the same time.
    #nodecmtys2 = cmtys2.nodecmtys()
    nodecmtys2 = defaultdict(set)
    for cname2, nodes2 in cmtys2.iteritems():
        sizes2[cname2] = len(nodes2)
        for n2 in nodes2:
            nodecmtys2[n2].add(cname2)
    # Main loop
    for cname1, nodes1 in cmtys1.iteritems():
        sizes1[cname1] = len(nodes1)
        for n1 in nodes1:
            for cname2 in nodecmtys2[n1]:
                if cname1 not in confusion:
                    confusion[cname1] = { }
                if cname2 not in confusion[cname1]:
                    confusion[cname1][cname2] = 0
                confusion[cname1][cname2] += 1
    #
    return confusion, sizes1, sizes2

def confusion(cmtys1, cmtys2):
    """Return the confusion matrix between two communities.

    Communities can include overlapping nodes.
    This function is fully general.

    Returns: dict of cmty_id_1->(dicts of cmty_id_2->overlap count)
    """
    d = _get_data(cmtys1, cmtys2)
    return d[0]

# Simple python-based implementations. These are naive O(n^2)
# implementations, whose only purpose is for pedagogy and unit-testing
# comparison with the efficient python2 implementations.
#
def mutual_information_python(cmtys1, cmtys2):
    MI = 0.0
    N = cmtys1.N
    assert N == cmtys2.N
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
    N = float(cmtys.N)
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
    
def vi_norm_python(cmtys1, cmtys2):
    """Normalized Variation of Information.

    This is variation of information, divided by log(N), which is the
    upper bound.
    """
    VI = vi_python(cmtys1, cmtys2)
    # mutual_information_python2 ensures that cmtys1.nodes ==
    # cmtys2.nodes.
    N = cmtys1.N
    NVI = VI / log2(N)
    # similarity would be 1 - NVI here
    return NVI
    
def nmi_python(cmtys1, cmtys2):
    I = mutual_information_python(cmtys1, cmtys2)
    Hs = entropy_python(cmtys1) + entropy_python(cmtys2)
    if Hs == 0 and I == 0:
        return 1.0
    if Hs == 0:
        return float('nan')
    In = 2.*I / Hs
    return In

def nmi_max_python(cmtys1, cmtys2):
    assert cmtys1.N == cmtys2.N
    if is_equal(cmtys1, cmtys2) or (cmtys1.q == 1 and cmtys2.q == 1):
        return 1.0
    I = mutual_information_python(cmtys1, cmtys2)
    H = max(entropy_python(cmtys1), entropy_python(cmtys2))
    return I/H
    
In = nmi_python # what is this for?
def _frac_detected(cmtys1, cmtys2, ovmode='overlap'):
    # overlap modes: overlap, newman
    return cmtys1.frac_detected(cmtys2, ovmode)
    
def jaccard_python(cmtys1, cmtys2):
    set1 = set()
    set2 = set()

    for cname, cnodes in cmtys1.iteritems():
        set1.update(frozenset((a,b)) 
                        for a,b in itertools.combinations(cnodes, 2))

    for cname, cnodes in cmtys2.iteritems():
        set2.update(frozenset((a,b)) 
                        for a,b in itertools.combinations(cnodes, 2))

    return len(set1 & set2) / float(len(set1 | set2))

# More efficient implementations using the confusion matrix.
def mutual_information_python2(cmtys1, cmtys2):
    N = cmtys1.N
    assert N == cmtys2.N
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

def vi_norm_python2(cmtys1, cmtys2):
    """Normalized Variation of Information.

    This is variation of information, divided by log(N), which is the
    upper bound.
    """
    VI = vi_python2(cmtys1, cmtys2)
    # mutual_information_python2 ensures that cmtys1.nodes ==
    # cmtys2.nodes.
    N = cmtys1.N
    NVI = VI / log2(N)
    return NVI

def vi_norm_python2_rev(cmtys1, cmtys2):
    return 1 - vi_norm_python2(cmtys1, cmtys2)
   
def nmi_python2(cmtys1, cmtys2):
    I = mutual_information_python2(cmtys1, cmtys2)
    Hs = entropy_python(cmtys1) + entropy_python(cmtys2)
    if Hs == 0 and I == 0:
        return 1.0
    if Hs == 0:
        return float('nan')
    In = 2.*I / Hs
    return In

def nmiG_python2(cmtys1, cmtys2):
    """NMI with geometric mean in denominator.

    Reference:

    Cluster Ensembles : A Knowledge Reuse Framework for Combining
    Multiple Partitions, A. Strehl and J. Ghosh, Journal of Machine
    Learning Research 3 (2002) 583-617
    """
    I = mutual_information_python2(cmtys1, cmtys2)
    Hs = entropy_python(cmtys1) * entropy_python(cmtys2)
    if Hs == 0 and I == 0:
        return 1.0
    if Hs == 0:
        return float('nan')
    In = I / sqrt(Hs)
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
        return -p * log2(p)
    def H(size, N):
        return h(size/N) + h((N-size)/N)
    def H2(c1size, c2size, overlap, N):
        hP11 = h((         overlap  )/N)
        hP10 = h((c1size - overlap  )/N)
        hP01 = h((c2size - overlap  )/N)
        hP00 = h((N -  c1size - c2size + overlap)/N)
        if hP11 + hP00 <= hP01 + hP10:
            # We want to exclude ones matching this condition. Return
            # 'inf' and it will be excluded from the min() function
            # wrapping this one.
            return float('inf')
        hPY1 = h(  (  c2size) / N)
        hPY0 = h(  (N-c2size) / N)
        # -(hPY1 + hPY0) = -H(c2size, N)
        return hP11+hP00+hP01+hP10 - hPY1 - hPY0

    # The main NMI loop. For every community in cmtys1, find the best
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
    detected communities. Recall measures how well the known
    communities are detected. Extra detected communities do not
    affect this result.
    """
    confusion = defaultdict(lambda: defaultdict(int))
    nodecmtys2 = cmtys2.nodecmtys()
    for cname, cnodes in cmtys1.iteritems():
        for node in cnodes:
            for c2 in nodecmtys2.get(node, []):
                confusion[cname][c2] += 1
    sizes1 = cmtys1.cmtysizes()
    sizes2 = cmtys2.cmtysizes()

    recls = []
    #sizes = []
    for c1name, others in confusion.iteritems():
        size = sizes1[c1name]
        #if weighted:
        #    sizes.append(size)
        others = [overlap/float(size + sizes2[c2name] - overlap)
                   for c2name, overlap in others.iteritems()]
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
    """Precision: how well cmtys1 matches cmtys2.

    Consider cmtys1 to be 'true' communities, and cmtys2 to be
    detected communities. Precision measures how well the detected
    communities all match the known. Extra known communities do not
    affect this result.
    """
    return recl_python2(cmtys2, cmtys1, weighted=weighted)
    
def F1_python2(cmtys1, cmtys2, weighted=True):
    """The harmonic mean of recall and precision.

    This is a symmetric measure.
    """
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
    
def count_pairs(cmtys1, cmtys2):
    """Counts pairs that are used in pair counting methods.
    a = pairs of nodes that are in the same community in both partitions
    b = pairs of nodes that are in the same community in cmtys1 and in
        different communities in cmtys2
    c = pairs of nodes that are in different communities in cmtys1 and
        in the same community in cmtys2
    d = pairs of nodes that are in different communities in both 
        partitions    
    all = all possible pairs
    """    
    # pair counting methods require nodes to be the same
    assert cmtys1.N == cmtys2.N
    
    pairs = {'a':0, 'b':0, 'c':0, 'd':0, 'all':0}
    ab = 0
    ac = 0

    confusion, sizes1, sizes2 = _get_data(cmtys1, cmtys2)
    a = sum(comb2(overlap) for c1name, others in confusion.iteritems()
                                     for overlap in others.itervalues())
    ab = sum(comb2(n) for n in sizes1.itervalues())
    ac = sum(comb2(n) for n in sizes2.itervalues())
    
    
    pairs['all'] = comb2(cmtys1.N)   
    pairs['a'] = a
    pairs['b'] = ab-a
    pairs['c'] = ac-a
    pairs['d'] = pairs['all'] - ab - ac + a
    
    return pairs
    
def jaccard_python2(cmtys1, cmtys2):
    """Jaccard index of two partitions.

    I do not have a reference for this measure, and the name is not
    very descriptive of what it does. It is said to be that
    implemented in radatools, but I have not checked this myself:
    http://deim.urv.cat/~sergio.gomez/radatools.php.

    Mathematical explanation:

    For cmtys1, take the set S1 of all pairs of nodes a,b where a and
    b are in the same community. Do the same for cmtys2 to produce
    S2. The pairs (a,b) are unordered.

    The jaccard score is  |S1 intersect S2| / |S1 union S2|

    Returns:
        the score (float)
    """

    # reusable functionality moved to count_pairs()    
      
#    same = 0
#    c1pairs = 0
#    c2pairs = 0
#
#    def comb2(n):
#        return n*(n-1)/2.0
#
#    confusion, sizes1, sizes2 = _get_data(cmtys1, cmtys2)
#    #for c1name, others in confusion.iteritems():
#    #    for c2name, overlap in others.iteritems():
#    #        same += comb2(overlap)
#    ovPairs = sum(comb2(overlap) for c1name, others in confusion.iteritems()
#                                     for overlap in others.itervalues())
#    c1pairs = sum(comb2(n) for n in sizes1.itervalues())
#    c2pairs = sum(comb2(n) for n in sizes2.itervalues())
#
#    return ovPairs / float(c1pairs + c2pairs - ovPairs)
    
    pairs = count_pairs(cmtys1, cmtys2)
    return pairs['a'] / (pairs['all'] - pairs['d'])

# for igraph implementation comparisons
def rand_python(cmtys1, cmtys2):
    pairs = count_pairs(cmtys1, cmtys2)
    return (pairs['a'] + pairs['d']) / pairs['all']
    
def adjusted_rand_python(cmtys1, cmtys2):
    # method should return 1.0 for the same partitions, but there are
    # some exceptions where the method isn't defined, so let's just
    # return 1.0 for the same partitions always    
    assert cmtys1.N == cmtys2.N

    if is_nodewise_equal(cmtys1, cmtys2) or (cmtys1.q == 1 and cmtys2.q == 1):
        return 1.0
        
    pairs = count_pairs(cmtys1, cmtys2)
    ab = pairs['a'] + pairs['b']
    ac = pairs['a'] + pairs['c']
    M = float(ab*ac)/pairs['all']
    
    return (pairs['a'] - M) / (0.5*(ab+ac) - M)

def omega_index_python(cmtys1, cmtys2):
    assert cmtys1.N == cmtys2.N
    
    if (cmtys1.q == 1 and cmtys2.q == 1) or is_nodewise_equal(cmtys1, cmtys2):
        # need to check if this measure would scale so that same
        # partitions would mean value 1.0 rather than returning 'nan'
        return float('nan')
    
    pairs = {frozenset((a,b)): 0
                        for a,b in itertools.combinations(cmtys1.nodes, 2)}

    pair_freq1 = deepcopy(pairs)
    pair_freq2 = deepcopy(pairs)
    
    for cname, cnodes in cmtys1.iteritems():
        for i in (frozenset((a,b)) 
                for a,b in itertools.combinations(cnodes, 2)):
            pair_freq1[i] += 1

    for cname, cnodes in cmtys2.iteritems():
        for i in (frozenset((a,b)) 
                for a,b in itertools.combinations(cnodes, 2)):
            pair_freq2[i] += 1 
    
    # compute how many times each node pair with nodes in the same
    # community occur
    
    # reverse mapping of the pair frequency (t_j(C) as a dict with j as
    # the key)
    t1 = {}
    for k,v in pair_freq1.iteritems():
        t1.setdefault(v, set()).add(k)
    
    t2 = {}
    for k,v in pair_freq2.iteritems():
        t2.setdefault(v, set()).add(k)
    
    omega_u = 0.0
    omega_e = 0.0
    for i in xrange(max(cmtys1.q, cmtys2.q)+1):
        omega_u += len(t1.get(i, set()) & t2.get(i, set()))
        omega_e += len(t1.get(i, set())) * len(t2.get(i, set()))
    omega_u = omega_u / comb2(cmtys1.N)
    omega_e = omega_e / comb2(cmtys1.N)**2
    
    return (omega_u - omega_e) / (1 - omega_e)
    
def fowlkes_mallows_python(cmtys1, cmtys2):
    pairs = count_pairs(cmtys1, cmtys2)
    ab = pairs['a'] + pairs['b']
    ac = pairs['a'] + pairs['c']
    if not (ab and ac): return float('nan') 
    return pairs['a'] / sqrt(ab*ac)

# this measure is not symmetric
def minkowski_coeff_python(cmtys1, cmtys2):
    # need to still check what should be returned if cmtys1==cmtys2 and 
    # q==N
    if cmtys1.q == cmtys1.N:
        return float('nan')
    pairs = count_pairs(cmtys1, cmtys2)
    
    if not (pairs['a'] or pairs['b']):
        return float('nan')
    
    return sqrt((pairs['b'] + pairs['c']) / (pairs['a'] + pairs['b']))


def gamma_coeff_python(cmtys1, cmtys2):
    pairs = count_pairs(cmtys1, cmtys2)
    ab = pairs['a'] + pairs['b']
    ac = pairs['a'] + pairs['c']
    
    if cmtys1.q == 1 or cmtys2.q == 1 or \
        cmtys1.q == cmtys1.N or cmtys2.q == cmtys2.N:
        return float('nan')
    elif is_nodewise_equal(cmtys1, cmtys2):
        return (pairs['all']-pairs['a'])/pairs['d']

    return ((pairs['all'] * pairs['a'] - ab*ac) / 
                sqrt(ab*ac*(pairs['all']-ab)*(pairs['all']-ac)))

# there is something similar in cmty.cmty_mapping but it doesn't
# seem to maximize the overlap sum
# there is probably some nicer solution for this (hungarian algorithm)
def classification_accuracy_python(cmtys1, cmtys2):
    assert cmtys1.N == cmtys2.N

    overlaps = numpy.zeros(shape=(cmtys1.q, cmtys2.q), dtype=int)
    
    for i1, (c1, n1) in enumerate(cmtys1.iteritems()):
        for i2, (c2, n2) in enumerate(cmtys2.iteritems()):
            overlaps[i1,i2] = len(set(n1) & set(n2))

    cmtys2_indx = range(cmtys2.q)
    def recr(i, j, prev):
        left = [n for n in cmtys2_indx if n not in prev]
        if left and i+1 < cmtys1.q:
            for k in left:
               for l in recr(i+1, k, prev+[k]):
                   yield l
        else:
            sum_ = 0
            for i_,j_ in enumerate(prev):
                sum_ += overlaps[i_,j_]     
            yield sum_
    
    max_ = 0
    for j in xrange(cmtys2.q):
        maxpath = max(recr(0, j, [j]))
        if maxpath > max_:
            max_ = maxpath
    
    return float(max_)/cmtys1.N

def classification_error_python(cmtys1, cmtys2):
    return 1 - classification_accuracy_python(cmtys1, cmtys2)

def classification_accuracy_python2(cmtys1, cmtys2):
    assert cmtys1.N == cmtys2.N

    overlaps = numpy.zeros(shape=(cmtys1.q, cmtys2.q), dtype=int)    

    for i1, (c1, n1) in enumerate(cmtys1.iteritems()):
        for i2, (c2, n2) in enumerate(cmtys2.iteritems()):
            overlaps[i1,i2] = len(set(n1) & set(n2))

    maxq = max(cmtys1.q, cmtys2.q) 
    minq = min(cmtys1.q, cmtys2.q)
    dq = maxq-minq
    
    # - overlaps, because we're trying to find maximum path
    M = -overlaps
    if dq:
        if cmtys2.q > cmtys1.q:
            shape = (dq,maxq)
            ax = 0
        else:    
            shape = (maxq,dq)
            ax = 1
        # add -maxval rows/columns to fill the matrix to square matrix
        z = numpy.full(shape, M.min())
        M = numpy.concatenate((M,z), axis=ax)
    

    def step3_1(mat):
        assigned = (mat == 0)
        assigned_rows = set()
        assigned_cols = set()
        
        row_q = range(maxq)
        while row_q:
            row_q.sort(key=lambda x: len(numpy.where(assigned[x,:])[0]))
            i = row_q.pop(0)
            
            assigned_rows.add(i)
            row = assigned[i,:]
            
            cols = set(numpy.where(row)[0])
            if not (cols - assigned_cols): continue
           
            j = list(cols - assigned_cols)[0]
            # choose column with least amount of nonassigned zeros
            t = maxq
            for c in (cols - assigned_cols):
                col_zeros = len(numpy.where(assigned[:,c])[0])
                if col_zeros < t:
                    j = c
                    t = col_zeros
            assigned_cols.add(j)

            col = assigned[:, j]
            row[:] = False
            col[:] = False
            row[j] = True    
        
        return assigned
    
    def step3_2(mat, assigned):
        zeros = (mat == 0)
        marked_cols = set()
        unmarked_rows = set(range(maxq))
        # 'drawing lines'
        
        newly_marked_rows = [i for i,row in enumerate(assigned) \
                                        if numpy.count_nonzero(row)==0]
        while newly_marked_rows:
            i = newly_marked_rows.pop()
            unmarked_rows.discard(i)
            for c in numpy.where(zeros[i,:])[0]:
                if c in marked_cols: continue
                marked_cols.add(c)
                for r in numpy.where(assigned[:, c])[0]:
                    if r in unmarked_rows:
                        newly_marked_rows.append(r)
                        
        return unmarked_rows, marked_cols
        
    def step4(rows, cols):
        if len(rows) + len(cols) == maxq:
            return True
        return False
    
    def step5(mat, rows, cols):
        # find minimum value from elements left over
        leftover = numpy.ones((maxq,maxq))
        # addition matrix
        addition = -numpy.ones((maxq,maxq))
        # check all marked cols/unmarked rows, +1 to those
        # so double marked get 1 as the multiplier while marked get 0
        for r in rows:
            leftover[r, :] = 0
            addition[r, :] += 1
        for c in cols:
            leftover[:, c] = 0
            addition[:, c] += 1
        minval = numpy.min(mat[numpy.nonzero(leftover)])
        mat += minval*addition
        return mat
    
    # step1
    # find row minimum for each row and subtract it from the row       
    M = (M.T - M.min(axis=1)).T
    # check if we can assign
    a = step3_1(M)
    rs,cs = step3_2(M, a)
    # 'if not numpy.count_nonzero(a) == maxq' should probably be enough
    if not step4(rs,cs):
        # step2
        # find column minimum for each column and subtract it from the
        # column
        M = M - M.min(axis=0)
        
        a = step3_1(M)
        rs,cs = step3_2(M, a)
        while not step4(rs,cs):
            # step3
            # mark columns and rows in matrix, so that all the
            # zeros are covered with minimum amount of covered
            # columns and rows
            # step4
            # check number of covered columns and rows,
            # if maxq then assign, else continue to step5
            # step5
            # find smallest unmarked value and subtract it from
            # all the unmarked values and add it to all double-
            # marked values, goto step3
            M = step5(M,rs,cs)
            
            a = step3_1(M)
            rs,cs = step3_2(M, a)

    if numpy.count_nonzero(a) == maxq:
        sum_ = (a[:cmtys1.q, :cmtys2.q] * overlaps).sum()
    else:
        print "something went wrong"
        print "assigned:", numpy.count_nonzero(a)
        return 0
    
    return float(sum_)/cmtys1.N

def classification_error_python2(cmtys1, cmtys2):
    return 1 - classification_accuracy_python2(cmtys1, cmtys2)

def norm_van_dongen_python(cmtys1, cmtys2):
    assert cmtys1.N == cmtys2.N    
    overlaps = numpy.zeros(shape=(cmtys1.q, cmtys2.q), dtype=int)
    
    for i1, (c1, n1) in enumerate(cmtys1.iteritems()):
        for i2, (c2, n2) in enumerate(cmtys2.iteritems()):
            overlaps[i1,i2] = len(set(n1) & set(n2))

    return 1- (1.0/(2*cmtys1.N))*(sum([max(row) for row in overlaps]) +
                                 sum([max(col) for col in overlaps.T]))

def norm_van_dongen_python_rev(cmtys1, cmtys2):
    return 1 - norm_van_dongen_python(cmtys1, cmtys2)

# distance metrics do not handle overlapping currently
# need to check if they should/can
def calculate_meet(cmtys1, cmtys2):
    meet = set()
    for c1, n1 in cmtys1.iteritems():
        for c2, n2 in cmtys2.iteritems():
            if set(n1) & set(n2):
                meet.add(frozenset(set(n1) & set(n2)))
            
    return meet
    
#def distance_moved_python(cmtys1, cmtys2):
#    assert cmtys1.N == cmtys2.N
#    meet = calculate_meet(cmtys1, cmtys2)
#    # NO IDEA IF THIS WORKS GENERALLY
#    m = 2.0*len(meet) - (cmtys1.q + cmtys2.q) + 2*abs(cmtys1.q - cmtys2.q)
#    
#    return 1 - m/cmtys1.N

def distance_moved_python(cmtys1, cmtys2):
    assert cmtys1.N == cmtys2.N
    meet = calculate_meet(cmtys1, cmtys2)
    mov = 0.0
    
    for c1,n1 in cmtys1.iteritems():
        n = set(n1)
        subsets = []
        for m in meet:
            if not n: break                
            if m.issubset(n) and m != n:
                subsets.append(m)
        subsets.sort(key=lambda x: len(x), reverse=True)
        if subsets:
            for set_ in subsets[1:]:
                mov += len(set_)
        
    for c2,n2 in cmtys2.iteritems():
        n = set(n2)
        subsets = []
        for m in meet:
            if not n: break                
            if m.issubset(n) and m != n:
                subsets.append(m)
        subsets.sort(key=lambda x: len(x), reverse=True)
        if subsets:
            for set_ in subsets[1:]:
                mov += len(set_)
    
    return 1 - mov/cmtys1.N
    
def distance_division_python(cmtys1, cmtys2):
    assert cmtys1.N == cmtys2.N
    meet = calculate_meet(cmtys1, cmtys2)
    d = 0.0
    
    for c1,n1 in cmtys1.iteritems():
        n = set(n1)
        for m in meet:
            if not n: break                
            if m.issubset(n) and m != n:
                d += 1
                n = n-m
    for c2,n2 in cmtys2.iteritems():
        n = set(n2)
        for m in meet:
            if not n: break
            if m.issubset(n) and m != n:
                d += 1
                n = n-m
                
    return 1 - d/cmtys1.N


# binary in- and exclusivity and scaled inclusivity
# only give scores for node classification similarity
def compute_cluster_comparison_matrix(cmtys1, cmtys2):
    assert cmtys1.N == cmtys2.N
    # cmtys1 is the reference partition
    X = numpy.zeros(shape=(cmtys2.q, cmtys1.q), dtype=float)
    cmty_pairs = []
    
    for i1, (c1, n1) in enumerate(cmtys2.iteritems()):
        row = []
        for i2, (c2, n2) in enumerate(cmtys1.iteritems()):
            set1 = set(n1)
            set2 = set(n2)
            X[i1,i2] = len(set1 & set2)**2 / (len(set1)*len(set2))
            row.append((n1,n2))
        cmty_pairs.append(row)
    return X,cmty_pairs

# what to do in the case of two equally good clusters?, which nodes get points
def binary_exclusive_python(cmtys1, cmtys2):
    X, cmty_pairs = compute_cluster_comparison_matrix(cmtys1, cmtys2)
    nodescores = {n:0 for n in cmtys2.nodes}
    col_maxs = numpy.argmax(X, axis=0)
    for j,i in enumerate(col_maxs):
        for n in cmty_pairs[i][j][0]:
            nodescores[n] += 1
            
    return nodescores

def binary_inclusive_python(cmtys1, cmtys2, threshold):
    X, cmty_pairs = compute_cluster_comparison_matrix(cmtys1, cmtys2)
    nodescores = {n:0 for n in cmtys2.nodes}
    rows, cols = numpy.where(X >= threshold)
    for i in rows:
        for j in cols:
            for n in cmty_pairs[i][j][0]:
                nodescores[n] += 1
                
    return nodescores

# can't do bias correction without comparing multiple partitions
# at the same time
def scaled_inclusivity_python(cmtys1, cmtys2):
    X, cmty_pairs = compute_cluster_comparison_matrix(cmtys1, cmtys2)
    nodescores = {n:0 for n in cmtys2.nodes}
    col_maxs = numpy.argmax(X, axis=0)
    for j,i in enumerate(col_maxs):
        for pair in cmty_pairs[i][j]:
            for n in set(pair[0]) & set(pair[1]):
                nodescores[n] += X[i,j]
            
    return nodescores

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
    input to the NMI code. The NMI code is very dumb, and does not
    remove comments, and uses only floating-point or integer node IDs,
    and does not give any type of error if these conditions are not
    met. Thus, only enable use_existing if you can vouch for the.

    check: bool, default True
        If true, check that all node names are either representable as
        integers or floating point numbers. This can either be python
        integers or floats, or strings of those.
    use_existing: bool, default False
        If true, and community object has a '.fname' attribute, and
        that file exists, use that as the input filename instead of
        re-writing the file. WARNING: if you use this, you must
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
            assert all(isinstance(x, (int, float))
                           or _is_float_str(x) for x in nodes )
        for nodes in cmtys2.itervalues():
            assert all(isinstance(x, (int, float)) 
                           or _is_float_str(x) for x in nodes )

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
        # Write community files, if they do not already exist. These
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
    'vi_norm': ['vi_norm_python', 'vi_norm_python2'],
    'mutual_information':
        ['mutual_information_python', 'mutual_information_pcd',
         'mutual_information_python2'],
    'nmi': ['nmi_python', 'nmi_igraph', 'nmi_pcd', 'nmi_python2'],
    'nmiG': ['nmiG_python2', ],
    'nmi_LFK': ['nmi_LFK_python2'],# 'nmi_LFK_LF', 'nmi_LFK_pcd',
                #'nmi_LFK_pcdpy'],
    'rand': ['rand_igraph', 'rand_python'],
    'adjusted_rand': ['adjusted_rand_igraph', 'adjusted_rand_python'],
    'F1':   ['F1_python2'],
    'recl': ['recl_python2'],
    'prec': ['prec_python2'],
    'F1_uw':   ['F1_uw_python2'],
    'recl_uw': ['recl_uw_python2'],
    'prec_uw': ['prec_uw_python2'],
    'jaccard': ['jaccard_python2', 'jaccard_python'],
    'omega': ['omega_index_python'],
    'nmi_max': ['nmi_max_python'],
    'fowlkes_mallows': ['fowlkes_mallows_python'],
    'minkowski': ['minkowski_coeff_python'],
    'gamma_coeff': ['gamma_coeff_python'],
    'classification_error': ['classification_error_python2'],
    'nvd': ['norm_van_dongen_python'],
    'distance_m': ['distance_moved_python'],
    'distance_d': ['distance_division_python']
    }

# Standard implementations

mutual_information = mutual_information_python2
vi = vi_python2
vi_norm = vi_norm_python2
nmi = nmi_python2
nmiG = nmiG_python2
nmi_LFK = nmi_LFK_python2
nmi_max = nmi_max_python

rand = rand_igraph
adjusted_rand = adjusted_rand_python
F1 = F1_python2
recl = recl_python2
prec = prec_python2
F1_uw = F1_uw_python2
recl_uw = recl_uw_python2
prec_uw = prec_uw_python2
jaccard = jaccard_python2
omega = omega_index_python
fowlkes_mallows = fowlkes_mallows_python
minkowski = minkowski_coeff_python
gamma_coeff = gamma_coeff_python
classification_error = classification_error_python2
nvd = norm_van_dongen_python
distance_m = distance_moved_python
distance_d = distance_division_python