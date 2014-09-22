from math import log

from pcd.util import Averager
import cmodels

log2 = lambda x: log(x, 2)


def cmtyIntersect(G1, c1, G2, c2):
    """Number of nodes in intersect of c1 and c2."""
    return cmodels.cmtyIntersect(G1._struct_p, c1, G2._struct_p, c2)
def cmtyUnion(G1, c1, G2, c2):
    """Number of nodes in union of c1 and c2."""
    return cmodels.cmtyUnion(G1._struct_p, c1, G2._struct_p, c2)


def mutual_information_python(G0, G1):
    """Calculate mutual information between two graphs.

    This is I, direct mutual information."""
    assert G0.N == G1.N
    N = G0.N

    MI = 0.0
    for c0 in [    c for c in range(G0.Ncmty) if G0.cmtyN[c] != 0]:
        for c1 in [c for c in range(G1.Ncmty) if G1.cmtyN[c] != 0]:
            # We correlate c0 in G0 and c1 in G1 according to the formula.
            n0 = G0.cmtyN[c0]
            n1 = G1.cmtyN[c1]

            # number of shared particles?
            n_shared = 0
            #for n in    G0.cmtyll[c0, :n0]:
            #    if n in G1.cmtyll[c1, :n1]:
            #        n_shared += 1
            s0 = set(G0.cmtyContents(c0))
            s1 = set(G1.cmtyContents(c1))
            n_shared = len(s0 & s1)
            #assert n_shared == len(s0 & s1)

            if n_shared == 0:
                continue

            MI += (n_shared/float(N)) * log2(n_shared*N/float(n0*n1))

    return MI
def mutual_information_c(G0, G1):
    return cmodels.mutual_information(G0._struct_p, G1._struct_p)
mutual_information = mutual_information_c
I = mutual_information

def _entropy(G):
    if G.oneToOne: return G.entropy
    return 1
def VI(G0, G1):
    """Variation of Information"""
    I = mutual_information(G0, G1)
    VI = G0.entropy_c + G1.entropy_c - 2*I
    return VI
def In(G0, G1):
    """Normalized Mutual Information"""
    I = mutual_information(G0, G1)
    Hs = G0.entropy_python + G1.entropy_python
    if Hs == 0 and I == 0:
        return 1.0
    In = 2.*I / Hs
    return In


#
# Overlap-using mutual information
#
def h(p):
    if p==0 or p==1: return 0
    return -p * log(p,2)
def H(G, c):
    N = float(G.N)
    return h(G.cmtyN[c]/N) + h((N-G.cmtyN[c])/N)
def p11(G1,G2,c1,c2): return (           len(cXm & cYm)  )/float(N)
def p10(G1,G2,c1,c2): return (len(cXm) - len(cXm & cYm)  )/float(N)
def p01(G1,G2,c1,c2): return (len(cYm) - len(cXm & cYm)  )/float(N)
def p00(G1,G2,c1,c2): return (  N      - len(cXm | cYm)  )/float(N)
def H2_python(GX, GY, cX, cY):
    """H(X_k | Y_l)."""
    cXm = set(GX.cmtyContents(cX)) # cmty X members
    cYm = set(GY.cmtyContents(cY)) # cmty Y members
    assert GX.N == GY.N
    N = float(GX.N)  # make it float so that division below works.
    #print "  p %d %d"%(len(cXm & cYm), len(cXm | cYm))
    hP11 = h((           len(cXm & cYm)  )/N)
    hP10 = h((len(cXm) - len(cXm & cYm)  )/N)
    hP01 = h((len(cYm) - len(cXm & cYm)  )/N)
    hP00 = h((  N      - len(cXm | cYm)  )/N)
    if hP11 + hP00 <= hP01 + hP10:
        # We want to exclude ones matching this condition.  Return
        # 'inf' and it will be excluded from the min() function
        # wrapping this one.
        return float('inf')
    hPY1 = h(  (  len(cYm)) / N)
    hPY0 = h(  (N-len(cYm)) / N)
    return hP11+hP00+hP01+hP10 - hPY1 - hPY0
def H2_c(GX, GY, cX, cY):
    return cmodels.H2(GX._struct_p, GY._struct_p, cX, cY)
H2 = H2_c

def HX_Ynorm_python(GX, GY):
    """ """
    HX_Y = Averager()
    for cX in GX.cmtys():
        HX_Yhere = min(H2(GX, GY, cX, cY) for cY in GY.cmtys())
        if HX_Yhere == float('inf'):
            HX_Yhere = H(GX, cX)
        _ = H(GX, cX)
        if _ == 0:
            HX_Y += 0
        else:
            HX_Y += HX_Yhere / _
    return HX_Y.mean
def HX_Ynorm_c(GX, GY, weighted=False):
    return cmodels.HX_Ynorm(GX._struct_p, GY._struct_p, weighted)
HX_Ynorm = HX_Ynorm_c



def mutual_information_overlap_python(G0, G1):
    """LF overlap-including normalized mutual information."""
    N = 1 - .5 * (HX_Ynorm_python(G0, G1) + HX_Ynorm_python(G1, G0) )
    return N
def mutual_information_overlap(G0, G1, weighted=False):
    """LF overlap-including normalized mutual information."""
    N = 1 - .5 * (HX_Ynorm(G0, G1, weighted) + HX_Ynorm(G1, G0, weighted) )
    return N
N = mutual_information_overlap

