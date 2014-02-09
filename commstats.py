
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

def lin_bin(k, range=1, bins=10):
    #print k, floor(k*bins/float(range)) * range/bins
    return floor(k*bins/float(range)) * range/bins
def lin_bin_width(k, range=1, bins=10):
    return range/float(bins)

def log_bin(k, base=10, decadesize=10, minlog=1):
    """Given a point, return its bin center"""
    decadesize = float(decadesize)
    if k < minlog:  return k
    i = log(k)/log(base)*decadesize
    i = round(i)
    k = exp((i/decadesize) * log(base))
    return k
#def log_bin(k, base=10, decadesize=10, minlog=10):
#    return k
def log_bin_width(k, base=10, decadesize=10, minlog=1):
    """Given a bin center, return its bin width"""
    if k < minlog: return 1
    i = log(k)/log(base)*decadesize
    i = round(i)
    max_ = exp((i+.5)/decadesize * log(base))
    min_ = max(minlog, exp((i-.5)/decadesize * log(base)))
    width = max_ - min_
    return float(width)



class Statter(object):
    """Base statter."""
    minsize = 2
    log_x = True
    log_y = False
    log_p = True
    xlabel = '$n_\mathrm{cmty}$'
    ylabel = None
    legend_loc = None
    @property
    def title(self):
        if hasattr(self, '_title'): return self._title
        return self.__class__.__name__
    @title.setter
    def title(self, title):
        self._title = title
    def __init__(self):
        self._data = collections.defaultdict(lambda: collections.defaultdict(list))
        self.label_order = [ ]
    def calc(self, g, cmtys):
        return self.calc(g, cmtys)

    def add_label(self, label):
        if label not in self._data:
            self.label_order.append(label)
            self._data[label]
    def accumulate(self, data, label):
        if label not in self._data:
            self.label_order.append(label)
        for n, sld in data:
            self._data[label][n].append(sld)
    def add(self, g, cmtys, label):
        self.accumulate(self.calc(g, cmtys, label=label))

    def write(self, fname, axopts={}):
        from pcd.support import matplotlibutil
        ax, extra = matplotlibutil.get_axes(fname, **axopts)

        data = self._data

        # Do the actual plotting
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colormap = cm.get_cmap('jet')
        normmap = mcolors.Normalize(vmin=0, vmax=len(data))

        label_order = self.label_order
        if not label_order:
            label_order = sorted(data.keys())
        for i, label in enumerate(label_order):
            points = data[label]
        #for i, (label, points) in enumerate(data.iteritems()):
            #print label, points.keys()
            xy = sorted((a, numpy.mean(b)) for a,b in points.iteritems())
            xvals, y = zip(*xy)
            ax.plot(xvals, y, color=colormap(normmap(i)),
                    label=label)

        if self.log_x:
            ax.set_xscale('log')
        #ax.set_xlim(min(ns), max(ns))
        #ax.set_ylim(*ylims)
        if self.xlabel: ax.set_xlabel(self.xlabel)
        if self.ylabel: ax.set_ylabel(self.ylabel)
        if self.title: ax.set_title(self.title)
        #ax.set_ylim(floor(min_v), ceil(max_v))
        #ax.set_yticks(range(2, int(ceil(max_v))))
        if self.log_y:
            ax.minorticks_on()
            ax.set_yscale('log', basey=10)

        if len(data) > 1:
            ax.legend(loc=self.legend_loc)

        self._hook_write_setup_plot(locals())
        matplotlibutil.save_axes(ax, extra)
    def _hook_write_setup_plot(self, lcls):
        pass

class QuantileStatter(Statter):
    quantiles = numpy.arange(0.0, 1.1, .1)
    def __init__(self):
        self._data = collections.defaultdict(list)
    def calc(self, g, cmtys):
        return self.calc(g, cmtys)

    def add_label(self, label):
        pass
        #if label not in self._data:
            #self.label_order.append(label)
            #self._data[label]
    def accumulate(self, data):
        for n, sld in data:
            self._data[n].append(sld)
    def add(self, g, cmtys):
        self.accumulate(self.calc(g, cmtys))

    def write(self, fname, axopts={}):
        from pcd.support import matplotlibutil
        ax, extra = matplotlibutil.get_axes(fname, **axopts)

        quantiles = [ ]
        p_values = self.quantiles
        xvals = [ ]
        for n, values in sorted(self._data.iteritems()):
            xvals.append(n)
            values.sort()
            quantiles.append(tuple(quantile(values, p)
                                   for p in p_values))
        p_quantiles = zip(*quantiles)
        # Do the actual plotting
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colormap = cm.get_cmap('jet')
        #normmap = mcolors.Normalize(vmin=0, vmax=len(p_quantiles)-1)
        normmap = mcolors.Normalize(vmin=0, vmax=(len(p_quantiles)-1)/2)
        def quant_norm(x):
            "Divide quantile in half"
            lambda x:  abs(x - (len(p_quantiles)-1)/2)

        min_v = min(p_quantiles[0])
        max_v = max(p_quantiles[-1])
        epsilon = (max_v - min_v)*.001
        def adjust(i):
            return epsilon * (i-len(p_quantiles)/2.)

        for i,qtl in enumerate(p_quantiles):
            #print colormap(normmap(i))
            qtl = numpy.asarray(qtl) + adjust(i)
            ax.plot(xvals, qtl, #color='blue')
                    color=colormap(normmap(quant_norm(i))))

        # plot mean
        ax.plot(xvals, [numpy.mean(x) for x in quantiles],
                  color='green', lw='3')
        ylims = ax.get_ylim()

        if self.log_x:
            ax.set_xscale('log')
        ax.set_xlim(min(xvals), max(xvals))
        ax.set_ylim(*ylims)
        if self.xlabel: ax.set_xlabel(self.xlabel)
        if self.ylabel: ax.set_ylabel(self.ylabel)
        if self.title: ax.set_title(self.title)
        #print min_v, max_v
        #ax.set_ylim(floor(min_v), ceil(max_v))
        #ax.set_yticks(range(2, int(ceil(max_v))))
        if self.log_y:
            ax.minorticks_on()
            ax.set_yscale('log', basey=10)

        self._hook_write_setup_plot(locals())
        matplotlibutil.save_axes(ax, extra)


class DistStatter(Statter):
    bin = True
    def write(self, fname, axopts={}):
        from pcd.support import matplotlibutil
        ax, extra = matplotlibutil.get_axes(fname, **axopts)

        data = self._data

        # Do the actual plotting
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colormap = cm.get_cmap('jet')
        normmap = mcolors.Normalize(vmin=0, vmax=len(data))

        import growsf
        from functools import partial

        if self.bin and self.log_y:
            binparams = { }
            binfunc = log_bin
            binwidth = log_bin_width
        elif self.bin:
            binparams = {'bins':20 }
            binfunc = lin_bin
            binwidth = lin_bin_width
        else:
            binparams = { }
            binfunc = lambda x: x
            binwidth = lambda x: 1

        label_order = self.label_order
        if not label_order:
            label_order = sorted(data.keys())
        for i, label in enumerate(label_order):
            points = data[label]
        #for i, (label, points) in enumerate(data.iteritems()):
            vals = []
            [vals.extend(x) for x in points.values() ]
            vals = map(partial(binfunc, **binparams), vals)
            vals, vals_counts = growsf.counts(vals)
            vals_counts = [ c/binwidth(s, **binparams)
                            for s,c in zip(vals, vals_counts) ]
            p_vals = growsf.norm(vals_counts)
            ax.plot(vals, p_vals, color=colormap(normmap(i)),
                    label=label)

        if self.log_y:
            ax.set_xscale('log')
        #ax.set_xlim(min(ns), max(ns))
        #ax.set_ylim(*ylims)
        if self.xlabel: ax.set_xlabel(self.ylabel)
        if self.ylabel: ax.set_ylabel("P(%s)"%self.ylabel)
        if self.title: ax.set_title(self.title)
        #ax.set_ylim(floor(min_v), ceil(max_v))
        #ax.set_yticks(range(2, int(ceil(max_v))))
        if self.log_p:
            ax.minorticks_on()
            ax.set_yscale('log', basey=10)

        if len(data) > 1:
            ax.legend(loc=self.legend_loc)

        self._hook_write_setup_plot(locals())
        matplotlibutil.save_axes(ax, extra)




class CmtySize(Statter):
    """Community sizes"""
    ylabel = "size"
    log_y = True
    def calc(self, g, cmtys):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            yield n_cmty, n_cmty


class Density(Statter):
    log_y = False
    ylabel = 'edge density'
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

            d = 2 * n_edges_subgraph / float(n_cmty*(n_cmty-1))

            n_cmty = log_bin(n_cmty)
            yield n_cmty, d


class ScaledLinkDensity(Statter):
    log_y = False
    ylabel = 'scaled link density'
    legend_loc = 'lower right'
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
        if isinstance(self, DistStatter):
            return
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
    legend_loc = 'lower right'
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
    legend_loc = 'lower right'
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


class CmtySelfNeighborFraction(Statter):
    """n_cmty / count(union(neighbors(node) for node in cmty)"""
    ylabel = "neighbor fraction"
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
            #neighbor_ratio = (len(neighbors)-n_cmty)/float(n_cmty)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, neighbor_ratio


class CmtyHubness(Statter):
    """max(internal_degree) / (n_cmty-1)."""
    ylabel = "cmty max hubness"
    legend_loc = 'lower right'
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

            sg = g.subgraph(cnodes)
            maxdeg = max(sg.degree().values())
            x = maxdeg / (float(n_cmty)-1)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, x
class CmtyAvgDegree(Statter):
    """sum(degree) / n_cmty.  Average total degree (internal + external)"""
    ylabel = "cmty average degree"
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

            sumdeg = sum(g.degree(cnodes).itervalues())
            x = sumdeg / float(n_cmty)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, x

class MaxDegNodeJaccard(Statter):
    """jaccard(community, neighbors of highest degree node)"""
    ylabel = "cmty jacc with hub node neighbors"
    n = 1
    legend_loc = 'lower right'
    def calc(self, g, cmtys):
        nodecmtys = cmtys.nodecmtys_onetoone()
        cmtygraph = cmtys.cmty_graph(g)
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            sg = g.subgraph(cnodes)
            nodes_by_degree = sorted(sg.degree(cnodes).iteritems(),
                                     reverse=True, key=lambda x: x[1])
            if n_cmty < self.minsize:
                continue

            neighs = set()
            for i in range(self.n):
                n = nodes_by_degree[i][0]
                neighs.update(g.neighbors(n))
            overlap = len(neighs & cnodes)
            jacc =  overlap / float(len(cnodes) + len(neighs) - overlap)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, jacc



#class MaxDegNodeFocusednessLabel(MaxDegNodeFocusedness,Statter):
#    pass


# For every class, make subclasses for Dist and Quantiles
_glbs = list(globals().iteritems())
for name, obj in _glbs:
    #print name, type(obj)
    try:
        if not issubclass(obj, object): continue
    except TypeError:
        continue
    if obj in  (Statter, QuantileStatter, DistStatter): continue
    if issubclass(obj, Statter) and not issubclass(obj, QuantileStatter):
        globals()[name+'Qtl'] = type(name+'Qtl', (obj,QuantileStatter), {})
    if issubclass(obj, Statter) and not issubclass(obj, DistStatter):
        globals()[name+'Dist'] = type(name+'Dist', (obj,DistStatter), {})
