
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

def log_bin(k, base=10, decadesize=10, minlog=10):
    decadesize = float(decadesize)
    if k < minlog:  return k
    i = log(k)/log(base)*decadesize
    i = round(i)
    k = exp((i/decadesize) * log(base))
    return k

class Statter(object):
    minsize = 2
    log_y = False
    xlabel = '$n_\mathrm{cmty}$'
    ylabel = None
    @property
    def title(self):
        if hasattr(self, '_title'): return self._title
        return self.__class__.__name__
    @title.setter
    def title(self, title):
        self._title = title
    def __init__(self):
        self._data = collections.defaultdict(list)
    def calc(self, g, cmtys):
        raise NotImplemented

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
        p_quantiles = zip(*quantiles)
        # Do the actual plotting
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colormap = cm.get_cmap('jet')
        normmap = mcolors.Normalize(vmin=0, vmax=len(p_quantiles)-1)

        min_v = min(p_quantiles[0])
        max_v = max(p_quantiles[-1])
        epsilon = (max_v - min_v)*.001

        for i,qtl in enumerate(p_quantiles):
            #print colormap(normmap(i))
            qtl = numpy.asarray(qtl) + epsilon * (i-len(p_quantiles)/2.)
            ax.plot(ns, qtl, color='blue')
                    #color=colormap(normmap(i)))

        # plot mean
        ax.plot(ns, [numpy.mean(x) for x in quantiles],
                  color='green', lw='3')
        ylims = ax.get_ylim()

        ax.set_xscale('log')
        ax.set_xlim(min(ns), max(ns))
        ax.set_ylim(*ylims)
        if self.xlabel: ax.set_xlabel(self.xlabel)
        if self.ylabel: ax.set_ylabel(self.ylabel)
        if self.title: ax.set_title(self.title)
        #print min_v, max_v
        #ax.set_ylim(floor(min_v), ceil(max_v))
        #ax.set_yticks(range(2, int(ceil(max_v))))
        if self.log_y:
            ax.minorticks_on()
            ax.set_yscale('log', basey=2)

        self._hook_write_setup_plot(locals())
        matplotlibutil.save_axes(ax, extra)
    def _hook_write_setup_plot(self, lcls):
        pass



class ScaledLinkDensity(Statter):
    log_y = True
    ylabel = 'scaled link density'
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
    def _hook_write_setup_plot(self, lcls):
        ax = lcls['ax']
        ns = lcls['ns']
        # plot bounds
        ax.plot(ns, [2]*len(ns), ls='--', color='green')
        ax.plot(ns, ns, ls='--', color='green')

        ymin, ymax = lcls['ylims']
        ax.set_ylim(min(ymin, 1.9), max(ymax, 5))


class CmtyMinEmbeddednessOne(Statter):
    log_y = False
    ylabel = 'community embeddedness with worst external cmty'
    def calc(self, g, cmtys):
        nodecmtys = cmtys.nodecmtys_onetoone()
        cmtygraph = cmtys.cmty_graph(g, nodecmtys=nodecmtys)
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue

            k_in = cmtygraph.node[cname]['weight'] * 2
            k_outs = [ dat['weight']
                       for neigh,dat in cmtygraph[cname].iteritems() ]
            x = min( k_in/float(k_in+k_out) for k_out in k_outs )
            n_cmty = log_bin(n_cmty)
            yield n_cmty, x
class CmtyEmbeddedness(Statter):
    log_y = False
    ylabel = 'community embeddedness'
    def calc(self, g, cmtys):
        nodecmtys = cmtys.nodecmtys_onetoone()
        cmtygraph = cmtys.cmty_graph(g)
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue

            k_in = cmtygraph.node[cname]['weight'] * 2
            k_outs = [ dat['weight']
                       for neigh,dat in cmtygraph[cname].iteritems() ]
            x = k_in/float(k_in+sum(k_outs))
            n_cmty = log_bin(n_cmty)
            yield n_cmty, x
class NodeEmbeddedness(Statter):
    log_y = False
    ylabel = 'node embeddedness'
    def calc(self, g, cmtys):
        nodecmtys = cmtys.nodecmtys_onetoone()
        cmtygraph = cmtys.cmty_graph(g)
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue
            n_cmty = log_bin(n_cmty)
            for n in cnodes:
                k_tot = g.degree(n)
                k_in = sum(1 for neigh in g.neighbors(n) if neigh in cnodes)
                x = k_in / float(k_tot)
                yield n_cmty, x

#class CmtyExtRatio(CmtyMaxExtRatio):
#    def calc(self, g, cmtys):
#        nodecmtys = cmtys.nodecmtys_onetoone()
#        cmtygraph = cmtys.cmty_graph(g)
#        for cname, cnodes in cmtys.iteritems():
#            n_cmty = len(cnodes)
#            if n_cmty < self.minsize:
#                continue
#
#            e_in = cmtygraph.node[cname]['weight']
#            e_outs = [ dat['weight']
#                       for neigh,dat in cmtygraph[cname].iteritems() ]
#            # Only difference from above:
#            ext_ratio = sum(e_outs) / (2.*e_in)
#            n_cmty = log_bin(n_cmty)
#            yield n_cmty, ext_ratio
#
#class CmtyExtRatio(Statter):
#    def calc(self, g, cmtys):
#        nodecmtys = cmtys.nodecmtys_onetoone()
#        cmtygraph = cmtys.cmty_graph(g)
#        for cname, cnodes in cmtys.iteritems():
#            n_cmty = len(cnodes)
#            if n_cmty < self.minsize:
#                continue
#
#            e_in = cmtygraph.node[cname]['weight']
#            e_outs = [ dat['weight']
#                       for neigh,dat in cmtygraph[cname].iteritems() ]
#            # Only difference from above:
#            ext_ratio = sum(e_outs) / (2.*e_in)
#            n_cmty = log_bin(n_cmty)
#            yield n_cmty, ext_ratio


class NeighborFraction(Statter):
    def calc(self, g, cmtys):
        nodecmtys = cmtys.nodecmtys_onetoone()
        cmtygraph = cmtys.cmty_graph(g)
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue

            neighbors = set(cnodes)  # make a copy
            assert id(neighbors) != id(cnodes)
            for n in cnodes:
                neighbors.update(g.neighbors_iter(n))
            neighbor_ratio = float(n_cmty)/len(neighbors)
            neighbor_ratio = (len(neighbors)-n_cmty)/float(n_cmty)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, neighbor_ratio
