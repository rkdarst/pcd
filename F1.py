# Richard Darst, December 2011.

import ctypes
import numpy

import pcd
import pcd.cmodels as cmodels
#import pcd.multiresolution

def PO_one(G0, c0, G):
    assert G0.N == G.N
    c0contents = numpy.asarray([G0.cmtyContains(c0, n) for n in range(G0.N)],
                               dtype=bool)
    P = 0
    for c in G.cmtys():
        cContents = numpy.asarray([G.cmtyContains(c, n) for n in range(G.N)],
                                  dtype=bool)
        overlap = c0contents == cContents
        Noverlap = numpy.sum(overlap)
        #print c, Noverlap/float(G0.N)
        P = max(P, Noverlap/float(G0.N))
    return P
def PO(G0, G):
    return numpy.mean([PO_one(G0, c, G) for c in G0.cmtys()])

def F1_one(G0, c0, G):
    assert G0.N == G.N
    c0contents = numpy.asarray([G0.cmtyContains(c0, n) for n in range(G0.N)],
                               dtype=bool)
    P = 0
    for c in G.cmtys():
        cContents = numpy.asarray([G.cmtyContains(c, n) for n in range(G.N)],
                                   dtype=bool)
        overlap = c0contents == cContents
        # How many correct ones in c, compared to how many in c
        precision = numpy.sum(overlap[cContents])/float(numpy.sum(cContents))
        # how many correct ones recalled, per number total.
        recall    = numpy.sum(overlap[cContents])/float(numpy.sum(c0contents))
        #Noverlap = sum(overlap)
        #print c, Noverlap/float(G0.N)
        #if precision + recall == 0:
        #    continue
        F1 = 2 * precision * recall / (precision + recall)
        P = max(P, F1)
    return P
def F1c2(G0, c0, G):
    #return numpy.mean([F1_one(G0, c, G) for c in G0.cmtys()])
    precision = cmodels.c_double(0)
    recall = cmodels.c_double(0)
    F1 = cmodels.F1_one(G0._struct_p, c0, G._struct_p,
                        ctypes.pointer(precision), ctypes.pointer(recall))
    return F1, precision.value, recall.value

def F1python(G0, G):
    ret = [ F1c2(G0, c0, G) for c0 in G0.cmtys() ]
    vals = zip(*ret)
    return (numpy.mean(vals[0]),  # F1
            numpy.mean(vals[1]),  # precision
            numpy.mean(vals[2]),) # recall
def F1c(G0, G):
    F1        = numpy.zeros(dtype=ctypes.c_double, shape=G0.q)
    precision = numpy.zeros(dtype=ctypes.c_double, shape=G0.q)
    recall    = numpy.zeros(dtype=ctypes.c_double, shape=G0.q)
    cmodels.F1(G0._struct_p, G._struct_p,
               F1.ctypes.data,
               precision.ctypes.data,
               recall.ctypes.data,
               G0.q
               )
    #assert abs(F1.mean() - F1python(G0, G)[0]) < 1e-5
    return F1.mean(), precision.mean(), recall.mean()
F1 = F1c

def calc_P0(self, data, settings):
    Gs = data['Gs']
    overlapGs = data['ovGs']

    POs = [PO(G0, G_) for G_ in Gs ]
    PO_ovs = [PO(G0, G_) for G_ in overlapGs ]
    PO_mean = numpy.mean(POs)
    PO_std  = numpy.std(POs)
    PO_ov_mean = numpy.mean(PO_ovs)
    PO_ov_std  = numpy.std(PO_ovs)
    return [('PO',    PO_mean),
            ('PO_std',    PO_std),
            ('PO_ov', PO_ov_mean),
            ('PO_ov_std', PO_ov_std),
            ]
#pcd.MultiResolutionCorrelation.calcMethods.append(calc_P0)
def calc_F1(self, data, settings):
    G0 = getattr(self, 'G0', None)
    Gs = data['Gs']
    overlapGs = data.get('ovGs', None)
    returns = [ ]


    if G0:
        F1_ret = [ F1(G0, G_) for G_ in Gs ]
        F1s, precs, recls = zip(*F1_ret)
        returns.extend([
            ('F1',        numpy.mean(F1s)),
            ('F1_std',    numpy.std (F1s)),
            ('F1_prec',   numpy.mean(precs)),
            ('F1_recl',   numpy.mean(recls)),
            ])
        if overlapGs:
            F1_ov_ret = [ F1(G0, G_) for G_ in overlapGs ]
            F1_ovs, prec_ovs, recl_ovs = zip(*F1_ov_ret)
            returns.extend([
                ('F1_ov',        numpy.mean(F1_ovs)),
                ('F1_ov_std',    numpy.std (F1_ovs)),
                ('F1_ov_prec',   numpy.mean(prec_ovs)),
                ('F1_ov_recl',   numpy.mean(recl_ovs)),
                ])

    s_F1_ret =      [F1(Gs[x], Gs[y])      for x,y in data['pairIndexes'] ]
    s_F1_ret.extend([F1(Gs[y], Gs[x])      for x,y in data['pairIndexes'] ])
    s_F1s, s_precs, s_recls = zip(*s_F1_ret)
    returns.extend([
        ('s_F1',        numpy.mean(s_F1s)),
        ('s_F1_std',    numpy.std (s_F1s)),
        ('s_F1_prec',   numpy.mean(s_precs)),
        ('s_F1_recl',   numpy.mean(s_recls)),
        ])

    if overlapGs:
        s_F1_ov_ret =      [ F1(overlapGs[x], overlapGs[y])
                             for x,y in data['pairIndexes'] ]
        s_F1_ov_ret.extend([ F1(overlapGs[y], overlapGs[x])
                             for x,y in data['pairIndexes'] ])
        s_F1_ovs, s_prec_ovs, s_recl_ovs = zip(*s_F1_ov_ret)
        returns.extend([
            ('s_F1_ov',        numpy.mean(s_F1_ovs)),
            ('s_F1_ov_std',    numpy.std (s_F1_ovs)),
            ('s_F1_ov_prec',   numpy.mean(s_prec_ovs)),
            ('s_F1_ov_recl',   numpy.mean(s_recl_ovs)),
            ])
    return returns
#pcd.MultiResolution.calcMethods.append(calc_F1)
