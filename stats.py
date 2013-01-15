# Richard Darst, May 2012

import numpy

import collections

def G_info(G):
    print "N:", G.N
    print "Ncmty:", G.Ncmty
    print "q:", G.q
    print "all communities:", sorted(G.cmtys())
    print "community sizes:", [len(G.cmtyContents(c)) for c in G.cmtys()]
    print "community size mean:", numpy.mean([len(G.cmtyContents(c)) for c in G.cmtys()])
    print "community size std:", numpy.std([len(G.cmtyContents(c)) for c in G.cmtys()])
    print "oneToOne:", G.oneToOne

def g_stats(g):
    N = len(g)
    N_well_defined = 0
    N_int_greater_than_half = 0
    for n0 in g.nodes_iter():
        c0 = g.node[n0]['cmty']
        degree = 0
        degree_int = 0
        degree_ext_by_cmty = collections.defaultdict(int)
        for n1 in g.neighbors_iter(n0):
            c1 = g.node[n1]['cmty']
            degree += 1
            if c0 == c1:
                degree_int += 1
                continue
            degree_ext_by_cmty[c1] += 1
        if degree!=0 and degree_int / float(degree) > .5:
            N_int_greater_than_half += 1
        if degree!=0 and \
               (len(degree_ext_by_cmty)==0 or
               degree_int > max(degree_ext_by_cmty.itervalues())):
            N_well_defined += 1
    return dict(frac_apm_defined=(N_well_defined / float(N)),
                frac_raddicci_strong=(N_int_greater_than_half / float(N))
                )
