
import collections
import math
from math import ceil, floor, log, exp
import networkx
import numpy



def quantile(sorted_list, p):
    #N = len(sorted_list)
    #if p < .5/N: return sorted_list[0]
    #if p > (N-.5)/N: return sorted_list[-1]
    #k = int(floor((p-.5) * N))
    #pk = (k+.5)/N
    #return sorted_list[k] + N*(p-pk)*(sorted_list[k+1]-sorted_list[k])
    k = (len(sorted_list)-1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return sorted_list[int(k)]
    d0 = sorted_list[int(f)] * (c-k)
    d1 = sorted_list[int(c)] * (k-f)
    return d0+d1

def log_bin(k):
    b = 10
    decadesize = 20
    decadesize = float(decadesize)
    if k < 10:  return k
    i = log(k)/log(b)*decadesize
    i = round(i)
    k = exp((i/decadesize) * log(b))
    return k

class ScaledLinkDensity(object):
    minsize = 2
    def __init__(self):
        self._data = collections.defaultdict(list)
    def calc(self, g, cmtys):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue

            n_edges_subgraph = 0
            for n in cnodes:
                for nbr in g.adj[n]:
                    if nbr in cnodes:
                        n_edges_subgraph += 1
            assert n_edges_subgraph % 2 == 0
            n_edges_subgraph /= 2

            # subgraph for testing
            #sg = g.subgraph(cnodes)
            #assert networkx.is_connected(sg)
            #assert sg.number_of_edges() == n_edges_subgraph

            # compute sld.
            sld = 2 * n_edges_subgraph / float(n_cmty-1)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, sld

    def accumulate(self, data):
        for n, sld in data:
            self._data[n].append(sld)
    def add(self, g, cmtys):
        self.accumulate(self.calc(g, cmtys))

    def write(self, fname, axopts={}):
        from pcd.support import matplotlibutil
        ax, extra = matplotlibutil.get_axes(fname, **axopts)

        quantiles = [ ]
        p_values = numpy.arange(0.0, 1.1, .1)
        ns = [ ]
        for n, values in sorted(self._data.iteritems()):
            ns.append(n)
            values.sort()
            quantiles.append(tuple(quantile(values, p)
                                   for p in p_values))
        print ns
        p_quantiles = zip(*quantiles)
        # Do the actual plotting
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colormap = cm.get_cmap('jet')
        normmap = mcolors.Normalize(vmin=0, vmax=len(p_quantiles)-1)
        for i,qtl in enumerate(p_quantiles):
            #print colormap(normmap(i))
            qtl = numpy.asarray(qtl) + .02 * (i-len(p_quantiles)/2.)
            ax.loglog(ns, qtl, color='blue')
                    #color=colormap(normmap(i)))

        min_v = min(p_quantiles[0])
        max_v = max(p_quantiles[-1])
        # plot mean
        ax.loglog(ns, [numpy.mean(x) for x in quantiles],
                  color='green', lw='3')
        # plot bounds
        ax.loglog(ns, [2]*len(ns), ls='--', color='green')
        ax.loglog(ns, ns, ls='--', color='green')

        ax.set_xlim(min(ns), max(ns))
        print min_v, max_v
        #ax.set_ylim(floor(min_v), ceil(max_v))
        ax.minorticks_on()
        #ax.set_yticks(range(2, int(ceil(max_v))))
        ax.set_yscale('log', basey=2)

        matplotlibutil.save_axes(ax, extra)
