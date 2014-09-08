"""Community statistics functions.

This module contains classes which allow calculation of properties of communities.  There are two base classes:

- Statter - produces <measure(cmty)> vs size(cmty)
- DistStatter - produces histograpms of <measure(cmty)>
where `measure` is some general function of one community.

General use
~~~~~~~~~~~

sttr = CmtyDensity()   # initialize statter

# Add two averages for 'p=0.1'.  That these are two averages of the
# same thing is indicated by the same label.
sttr.add(g, cmtys_p1a, label='p=0.1')
sttr.add(g, cmtys_p1b, label='p=0.1')

# Add two averages for 'p=0.5'.  Since the label is different, there
# will be a different line on the final plot.
sttr.add(g, cmtys_p5a, label='p=0.5')
sttr.add(g, cmtys_p5b, label='p=0.5')

# Write to a file.  This uses options and formatting available in
# pcd.support.matplotlibutil.get_axes.
sttr.write('filename.[pdf,png]')

# Write to another matplotlib.axes.Axes object.  This could be used to
# add the plots to a figure with subfigures.  For doing final layout,
# you may want to set decorate=False to turn off all labels and
# titles, and add those yourself manually outside of the write
# function.
sttr.write(ax)


Subclasses
~~~~~~~~~~

Let's use CmtyDensity as an example.  CmtyDensity is a subclass of
Statter, and only defines thene extra pieces of data:

- calc(self, g, cmtys, ...)
  - The actual logic of calculating whatever this class does.  More on
  - this later.

- log_y = False
  ylabel = 'edge density'
  - These two class attributes provide configuration for the writing
    function, so that it can have good defaults for this particular
    type of plot.  For example, legend location, log scales, and so on
    can be configured, for both Statter and DistStatter writers.

At the end of the module, *every* class defined has some automatic
subclasses made.  For CmtyDensity, these wolud be:

- CmtyDensityQtl   (from QuantileStatter)
  - Like Statter (density vs size), but also draws quantiles
    (10-tiles) of the distribution, so that you can see the standard
    deviation, min, max, and skewness from the mean.

- CmtyDensityDist   (from DistStatter)
  - Made by mixing CmtyDensity with DistStatter

- CmtyDensityCount   (from CountStater)
  - Like CmtyDensityDist, but instead of normalizing, prints raw counts.

- CmtyDensityCcdf   (from CcdfStatter)
  - Community density distribution, as a complimentary cumulative
    distribution function

- CmtyDensityCccf   (from CccfStatter)
  - Community density distribution, as a complimentary cumulative
    counts function

You can look at the definitions of all of these statters, but all the
logic is in Statter or DistStatter.


Configuring output
~~~~~~~~~~~~~~~~~~

All of the writing logic is in the write() methods on Statter or
DistStatter, but these will probably be kind of confusing.  All
configuration is done by instance attributes.  For example, to change
the location of a legend, do::

  sttr.legend_loc = 'lower left'

You can get an idea of what is available by looking at the source for
these classes.

There are generally lots of weird hooks that can be used for
configuration, but for 'professional' layout you should draw to an
axis object and manually configure all labels and items.
"""

import collections
import math
from math import ceil, floor, log, exp, log10
import networkx
import numpy
import random


# Copied from growsf.py.  Move somewhere better sometime.
def counts(it, domain=None):
    counts = collections.defaultdict(int)
    for i in it:
        counts[i] += 1
    if domain is not None:
        keys = sorted(set(counts.iterkeys())|set(domain))
    else:
        keys = sorted(counts.iterkeys())
    values = [ counts[k] for k in keys ]
    return keys, values

def norm(it):
    sum_ = float(sum(it))
    return [ x/sum_ for x in it ]


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

def lin_bin(k, range=1, bins=10, ints=None, pointat=.5, offset=-.5):
    #print k, floor(k*bins/float(range)) * range/bins
    return (floor(k*bins/float(range)-offset)+offset) * range/bins + pointat*(range/float(bins))
def lin_bin_width(k, range=1, bins=10, ints=None, pointat=.5, offset=-.5):
    return range/float(bins)

def log_bin(k, base=10, decadesize=10, minlog=1,
            ints=False):
    """Given a point, return its bin center"""
    decadesize = float(decadesize)
    if k < minlog:  return int(k)
    i = log(k)/log(base)*decadesize
    i = round(i)
    k = exp((i/decadesize) * log(base))
    return k
#def log_bin(k, base=10, decadesize=10, minlog=10):
#    return k
def log_bin_width(k, base=10, decadesize=10, minlog=1,
                  ints=False):
    """Given a bin center, return its bin width"""
    if k < minlog: return 1
    i = log(k)/log(base)*decadesize
    i = round(i)
    max_ = exp((i+.5)/decadesize * log(base))
    min_ = max(minlog, exp((i-.5)/decadesize * log(base)))
    if ints:
        width = floor(max_) - ceil(min_) + 1
    else:
        width = max_ - min_
    return float(width)

# Pointers to original functions, in case lin_bin and log_bin have
# been overridden.
_lin_bin = lin_bin
_log_bin = log_bin
_lin_bin_width = lin_bin_width
_log_bin_width = log_bin_width


def cache_get(cache, name, func):
    """Simple dictionary based.

    Returns func(), or cache of the results.

    cache: dictionary or mappable.
        The cache.  If cache is None, then do not do any caching and
        just return the value.
    name: hashable
        Dictionary key to cache under.
    func: callable
        Callable to evaluate for value, only if value is not already
        in cache.
    """
    if cache is None:
        return func()
    if name not in cache:
        cache[name] = func()
    return cache[name]


class MultiStatter(object):
    def __init__(self,
                 sttrs,
                 layout):
        self.sttrs = sttrs
        self.layout = layout

        def add(self, g, cmtys, label, cache=None):
            for sttr in self.sttrs:
                sttr.accumulate(sttr.calc2(g, cmtys, cache=cache), label=label)
        def calc(self, g, cmtys, label, cache=None):
            results = [ ]
            for sttr in self.sttrs:
                results.append(sttr.calc(g=g, cmtys=cmtys, label=label,
                                         cache=cache))
            return results
        def accumulate(self, data, label):
            for sttr, data in zip(self.sttrs, data):
                sttr.accumulate(data, label)
        def write(self, fname, axopts, title=None):
            from string import lowercase

            ymax = max(len(_) for _ in self.layout)
            xmax = len(self.layout)
            figsize = (y_max*13, x_max*13)
            #figargs=dict(subplotpars=matplotlib.figure.SubplotParams(
            #   hspace=.3)))
            fig, extra = get_axes(fname=fname, ret_fig=True, figsize=figsize)


            for y, ylst in enumerate(self.layout):
                for x, sttr in enumerate(ylst):
                    # Use None to not write anything.
                    if sttr is None: continue
                    # Format for specifying statters:
                    # - integer
                    # - "N.write_function"
                    if isinstance(sttr, str):
                        i = int(sttr)
                    if isinstance(sttr, int):
                        i = sttr
                    sttr = self.layout[i]

                    sax = fig.add_subplot(maxy, maxx, 1+maxy*y+x)
                    sax.set_title(cdname, fontsize=20)
                    sax.text(0, 1.03, '(%s)'%lowercase[i%26], fontsize=20,
                             transform=sax.transAxes)
                    sttr.write(ax=sax)
                    #sttr.data.clear()

        save_axes(fig, extra)




class Statter(object):
    """Base statter."""
    minsize = 2
    log_x = True
    log_y = False
    log_p = True
    xlabel = '$n_\mathrm{cmty}$'
    ylabel = None
    legend_loc = None
    legend_kwargs = { }  # to set fontsize, use prop=dict(size=10)
    bin_ints = False
    xlim = None
    ylim = None
    legend_make = True
    markers = None
    # decorate: put all the title, axis labels, etc, on the plot.
    decorate = True
    # Pre-defined fixed colors for all data.
    colormap = None
    # Plotstyle dictionary
    plotstyle = collections.defaultdict(lambda: '-o')
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
    # Pickle support.  The lambda functions in self._data
    # collections.defaultdict are not pickleable, so convert to and
    # from regular dictionary when pickling.
    def __getstate__(self):
        state = self.__dict__.copy()
        state['_data'] = { }
        state['_data'].update(self._data)
        return state
    def __setstate__(self, state):
        _data = collections.defaultdict(lambda: collections.defaultdict(list))
        _data.update(state['_data'])
        state['_data'] = _data
        self.__dict__ = state

    def calc(self, g, cmtys, cache=None):
        """Calculate values.  Does not change object.

        See add() docstring."""
        raise NotImplementedError("This is a prototype.")
    def calc2(self, g, cmtys, cache=None):
        """Support for picking igraph calc() methods if exists."""
        if 'igraph' in str(type(g)):
            return self.calc_igraph(g, cmtys, cache=cache)
        return self.calc(g, cmtys, cache=cache)

    def add_label(self, label):
        if label not in self._data:
            self.label_order.append(label)
            self._data[label]
    def accumulate(self, data, label):
        """Append calculated values to internal lists.

        This adds the given data to self._data."""
        if label not in self._data:
            self.label_order.append(label)
        for n, sld in data:
            self._data[label][n].append(sld)
    def add(self, g, cmtys, label, cache=None):
        """Add data.  This does calc() and accumulate().

        calc() and accumulate() were separated in order to separate
        calculation from data collection, for possible future
        parallelism."""
        self.accumulate(self.calc2(g, cmtys, cache=cache), label=label)

    quantiles = None
    def write(self, fname, axopts={}, title=None):
        """Write data.

        `fname` can be an matplotlib.axes.Axes object, in which case
        no new files are created and lines are added to this axis.
        """
        from pcd.support import matplotlibutil
        ax, extra = matplotlibutil.get_axes(fname, **axopts)

        data = self._data

        # Do the actual plotting
        # Get colors and styles.
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colormap = cm.get_cmap('jet')
        normmap = mcolors.Normalize(vmin=0, vmax=len(data))

        label_order = self.label_order
        if not label_order:
            label_order = sorted(data.keys())
        for i, label in enumerate(label_order):
            # Colors
            if self.colormap and label in self.colormap:
                color = self.colormap[label]
            else:
                color = colormap(normmap(i))
            # Line type
            if label in self.plotstyle:
                plotstyle = self.plotstyle[label]
            elif self.markers:
                plotstyle = '-%s'%(self.markers[i%len(self.markers)])
            else:
                plotstyle = self.plotstyle[label]
            # Make all points data.  `points` is a list of (cmty_size,
            # [list of <measure> values]) pairs and `xy` is a list of
            # (cmty_size, <mean_measure>) pairs.
            points = data[label]
            if len(points) == 0:
                print "skipping %s with no data"%label
                continue
            xy = sorted((a, numpy.mean(b)) for a,b in points.iteritems())
            xvals, y = zip(*xy)

            # Plot the mean:
            ax.plot(xvals, y, plotstyle, lw=2, color=color, label=label)

            # If we are a subclass of QuantileStatter, plot quantiles,
            # too.
            if self.quantiles is not None and len(self.quantiles):
                quantiles = [ ]
                for x in xvals:
                    values = sorted(points[x])
                    quantiles.append(tuple(quantile(values, p)
                                           for p in self.quantiles))
                p_quantiles = zip(*quantiles)
                # This does a little bit of "offset" to make the
                # quantiles more visible.
                #min_v = min(p_quantiles[0])
                #max_v = max(p_quantiles[-1])
                #epsilon = (max_v - min_v)*.001
                #def adjust(i):
                #    return epsilon * (i-len(p_quantiles)/2.)
                adjust = lambda i: 0
                # Do the actual plotting
                for i, qtl in enumerate(p_quantiles):
                    #print colormap(normmap(i))
                    qtl = numpy.asarray(qtl) + adjust(i)
                    ax.plot(xvals, qtl, color=color,
                            lw=.5)

        # General plot layout and configuration.
        ylims = ax.get_ylim()
        if self.log_x:
            ax.set_xscale('log')
        if self.log_y:
            ax.minorticks_on()
            ax.set_yscale('log', basey=10)
        if self.xlim: ax.set_xlim(*self.xlim)
        if self.ylim: ax.set_ylim(*self.ylim)
        if self.decorate:
            if self.xlabel: ax.set_xlabel(self.xlabel)
            if self.ylabel: ax.set_ylabel(self.ylabel)
            if   title:      ax.set_title(self.title)
            elif self.title: ax.set_title(self.title)
            if len(data) > 1 and self.legend_make:
                ax.legend(loc=self.legend_loc, **self.legend_kwargs)


        self._hook_write_setup_plot(locals())
        matplotlibutil.save_axes(ax, extra)
    def _hook_write_setup_plot(self, lcls):
        pass

class QuantileStatter(Statter):
    quantiles = numpy.arange(0.0, 1.1, .1)


class DistStatter(Statter):
    bin = True
    dist_is_counts = False
    ylabel_dist = None
    dist_xlim = None
    dist_ylim = None
    domain = None
    dist_is_ccdf = False  # plot compl. cumul. dist func
    dist_is_cccf = False  # plot compl. cumul. count func
    def write(self, fname, axopts={}, title=None):
        from pcd.support import matplotlibutil
        ax, extra = matplotlibutil.get_axes(fname, **axopts)

        data = self._data

        # Do the actual plotting
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colormap = cm.get_cmap('jet')
        normmap = mcolors.Normalize(vmin=0, vmax=len(data))

        from functools import partial

        if self.bin and self.log_y:
            binparams = {'ints': self.bin_ints}
            binfunc = log_bin
            binwidth = log_bin_width
        elif self.bin:
            binparams = {'bins':20, 'ints': self.bin_ints}
            binfunc = lin_bin
            binwidth = lin_bin_width
        else:
            binparams = {'ints': self.bin_ints}
            binfunc = lambda x: x
            binwidth = lambda x: 1
        if hasattr(self, 'binparams'):
            binparams.update(self.binparams)

        label_order = self.label_order
        if not label_order:
            label_order = sorted(data.keys())
        for i, label in enumerate(label_order):
            points = data[label]
            if self.colormap and label in self.colormap:
                color = self.colormap[label]
            else:
                color = colormap(normmap(i))
            # Line type
            if label in self.plotstyle:
                plotstyle = self.plotstyle[label]
            elif self.markers:
                plotstyle = '-%s'%(self.markers[i%len(self.markers)])
            else:
                plotstyle = self.plotstyle[label]
            #
            vals = []
            [vals.extend(x) for x in points.values() ]
            vals = map(partial(binfunc, **binparams), vals)
            # Set a domain of all values to include
            if self.domain is not None:
                domain = set(binfunc(x, **binparams) for x in self.domain)
            else:
                domain = set()
            vals, vals_counts = counts(vals, domain=domain)
            if self.dist_is_counts:
                lines = ax.plot(vals, vals_counts, plotstyle, color=color,
                                label=label)
            # Complimentary cumulative distribution function plotting.
            elif self.dist_is_ccdf or self.dist_is_cccf:
                vals, vals_counts = \
                      zip(*sorted(zip(vals, vals_counts), key=lambda x: x[0])
                          )
                total = float(sum(vals_counts))
                import pcd.util
                counts_cumul = list(pcd.util._accumulate(vals_counts))
                counts_cumul = [counts_cumul[-1]-x for x in counts_cumul]
                if self.dist_is_ccdf:
                    counts_cumul = [ x/total for x in counts_cumul ]
                lines = ax.plot(vals, counts_cumul, plotstyle,
                                color=color, label=label)

            else:
                vals_counts = [ c/binwidth(s, **binparams)
                                for s,c in zip(vals, vals_counts) ]
                p_vals = norm(vals_counts)
                lines = ax.plot(vals, p_vals, plotstyle, color=color,
                                label=label)

            if hasattr(self, '_hook_lines'):
                self._hook_lines(locals())

        if self.log_y:
            ax.set_xscale('log')
        if self.dist_xlim: ax.set_xlim(*self.dist_xlim)
        if self.dist_ylim: ax.set_ylim(*self.dist_ylim)
        if self.log_p:
            ax.minorticks_on()
            ax.set_yscale('log', basey=10)
        if self.decorate:
            if self.ylabel: ax.set_xlabel(self.ylabel)
            if self.ylabel_dist:
                ax.set_ylabel(self.ylabel_dist)
            elif self.dist_is_counts:
                if self.ylabel: ax.set_ylabel("Count(%s)"%self.ylabel)
            elif self.dist_is_ccdf:
                if self.ylabel:
                    ax.set_ylabel("CCDF(%s)"%self.ylabel)
            elif self.dist_is_cccf:
                if self.ylabel:
                    ax.set_ylabel("Complimentary cumulative count(%s)"%self.ylabel)
            else:
                if self.ylabel: ax.set_ylabel("P(%s)"%self.ylabel)
            if   title:      ax.set_title(self.title)
            elif self.title: ax.set_title(self.title)
            if len(data) > 1 and self.legend_make:
                ax.legend(loc=self.legend_loc, **self.legend_kwargs)

        self._hook_write_setup_plot(locals())
        matplotlibutil.save_axes(ax, extra)

class CountStatter(DistStatter):
    dist_is_counts = True

class CcdfStatter(DistStatter):
    dist_is_ccdf = True
class CccfStatter(DistStatter):
    dist_is_cccf = True




class Null(Statter):
    def calc(self, g, cmtys, cache=None):
        # Force generation of possible lazy communites, don't do an
        # iteration.
        cmtys.iteritems()
        return ()


class NodeDeg(Statter):
    """Node degree"""
    xlabel = '1'
    ylabel = "$k$"
    log_y = True
    bin_ints = True
    binparams = dict(minlog=10)
    def calc(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            for node in cnodes:
                yield 1, g.degree(node)


class CmtySize(Statter):
    """Community sizes"""
    ylabel = "size"
    log_y = True
    bin_ints = True
    binparams = dict(minlog=10)
    def calc(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            yield n_cmty, n_cmty


class CmtyDensity(Statter):
    log_y = False
    ylabel = 'edge density'
    def calc(self, g, cmtys, cache=None):
        adj = g.adj
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue

            n_edges_subgraph = 0
            for n in cnodes:
                for nbr in adj[n]:
                    if nbr in cnodes and nbr != n:
                        n_edges_subgraph += 1
            assert n_edges_subgraph % 2 == 0
            n_edges_subgraph /= 2

            d = 2 * n_edges_subgraph / float(n_cmty*(n_cmty-1))

            n_cmty = log_bin(n_cmty)
            yield n_cmty, d


class ScaledLinkDensity(Statter):
    log_y = True
    ylabel = 'scaled link density'
    legend_loc = 'lower right'
    def calc(self, g, cmtys, cache=None):
        adj = g.adj
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue

            n_edges_subgraph = 0
            for n in cnodes:
                for nbr in adj[n]:
                    if n == nbr:
                        print 'self-loop'
                        continue
                    if nbr in cnodes:
                        n_edges_subgraph += 1
            assert n_edges_subgraph % 2 == 0, 'non-symmetric or self-loops?'
            n_edges_subgraph /= 2

            # subgraph for testing
            #sg = g.subgraph(cnodes)
            #assert networkx.is_connected(sg)
            #assert sg.number_of_edges() == n_edges_subgraph

            # compute sld.
            sld = 2 * n_edges_subgraph / float(n_cmty-1)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, sld
    def calc_igraph(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue
            sg = g.subgraph(cnodes)
            n_edges_subgraph = sg.ecount()
            sld = 2 * n_edges_subgraph / float(n_cmty-1)
            n_cmty = log_bin(n_cmty)
            yield n_cmty, sld
    def _hook_write_setup_plot(self, lcls):
        if isinstance(self, DistStatter):
            return
        if 'xvals' not in lcls:
            # There was no data, so don't do anything
            return
        ax = lcls['ax']
        ns = lcls['xvals']
        # plot bounds
        ax.plot(ns, [2]*len(ns), ls='--', color='green')
        ax.plot(ns, ns, ls='--', color='green')

        ymin, ymax = lcls['ylims']
        ax.set_ylim(min(ymin, 1.9), max(ymax, 5))

class AvgClusteringCoef(Statter):
    """Average clustering coefficient within community"""
    def calc(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue

            cc = networkx.average_clustering(g, cnodes)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, cc



class CmtyEmbeddedness(Statter):
    log_y = False
    ylabel = 'community embeddedness'
    legend_loc = 'lower right'
    # This allows us to transform embeddedness into conductance by
    # setting _val_map = lambda x: 1-x
    _val_map = staticmethod(lambda x: x)
    def calc(self, g, cmtys, cache=None):
        #cmtygraph = cmtys.cmty_graph(g)
        adj = g.adj
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            # Skip communities below some minimum size.  We must have
            # minsize at least 2, since edge density is not defined
            # for size=1.
            if n_cmty < self.minsize:
                continue
            # Old test code
            #_k_in2 = cmtygraph.node[cname]['weight'] * 2
            #_k_outs = [ dat['weight']
            #           for neigh,dat in cmtygraph[cname].iteritems() ]
            ##x = _k_in2/float(_k_in2+sum(_k_outs))
            k_tot = sum(len(adj[n]) for n in cnodes )
            k_in = sum(1 for n in cnodes
                       for n1 in adj[n] if n1 in cnodes)
            k_out = k_tot - k_in
            #assert k_in2 == _k_in and k_out == sum(_k_outs)
            x = k_in / float(k_in + k_out)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, self._val_map(x)
    def calc_igraph(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue
            sg = g.subgraph(cnodes)
            k_in = 2 * sg.ecount()
            k_tot = sum(g.degree(cnodes, loops=False))
            x = k_in / float(k_tot)
            n_cmty = log_bin(n_cmty)
            yield n_cmty, self._val_map(x)
CmtyEmbeddedness2 = CmtyEmbeddedness
class CmtyEmbeddednessWorst(Statter):
    """k_in / (k_in + k_out_max).

    k_in is every edge that go between two nodes in a community.

    k_out_max is the number of edges that start on a node in a
    community and end on a node that is in the other community.  If
    the two communities overlap, edges in that overlap are counted
    twice for this out degree.

    By this definition, two communities that exactly overlap will have
    an embeddedness of 0.5.  There is also no normalization as for the
    other community size.
    """
    log_y = False
    ylabel = 'community worst embeddedness'
    legend_loc = 'lower right'
    def calc(self, g, cmtys, cache=None):
        nodecmtys = cache_get(cache, 'nodecmtys', lambda: cmtys.nodecmtys())
        adj = g.adj
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            k_map = collections.defaultdict(int)
            k_int = 0
            k_tot = 0
            for n in cnodes:
                for neigh in adj[n]:
                    if neigh not in nodecmtys:
                        continue
                    for c2 in nodecmtys[neigh]:
                        k_tot += 1
                        if c2 == cname:
                            k_int += 1
                            continue
                        k_map[c2] += 1
            n_cmty = log_bin(n_cmty)

            if not k_map:
                yield n_cmty, self._val_map(1.0)
                continue
            if k_int == 0:
                print cname, cnodes, n_cmty, k_int, k_map, k_tot
                yield n_cmty, self._val_map(0)
                continue

            c_best = max(k_map, key=lambda _c: k_map[_c])
            x = k_int / float(k_int + k_map[c_best])

            yield n_cmty, self._val_map(x)

class CmtyConductance(CmtyEmbeddedness):
    """community conductance, 1-embeddedness"""
    ylabel = 'community conductance'
    _val_map = lambda x: 1.0-x
class CmtyConductanceWorst(CmtyEmbeddednessWorst):
    """community worst conductance, 1-worst embeddedness"""
    ylabel = 'community weak conductance'
    _val_map = lambda x: 1.0-x


class NodeEmbeddedness(Statter):
    log_y = False
    ylabel = 'node embeddedness'
    legend_loc = 'lower right'
    def calc(self, g, cmtys, cache=None):
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
#    def calc(self, g, cmtys, cache=None):
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
#    def calc(self, g, cmtys, cache=None):
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
    def calc(self, g, cmtys, cache=None):
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
    def calc(self, g, cmtys, cache=None):
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
    def calc_igraph(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue
            sg = g.subgraph(cnodes)
            maxdeg = g.maxdegree()
            x = maxdeg / (float(n_cmty)-1)
            n_cmty = log_bin(n_cmty)
            yield n_cmty, x
class CmtyAvgDegree(Statter):
    """sum(degree) / n_cmty.  Average total degree (internal + external)"""
    ylabel = "cmty average degree"
    log_y = True
    def calc(self, g, cmtys, cache=None):
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

#
# Statters related to overlaps
#
class CmtyOverlap(Statter):
    """Number of overlaps per community."""
    ylabel = "number of overlaps"
    log_y = True
    def calc(self, g, cmtys, cache=None):
        nodecmtys = cache_get(cache, 'nodecmtys', lambda: cmtys.nodecmtys())
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            overlaps = set.union(*(nodecmtys[n] for n in cnodes))
            yield log_bin(n_cmty), len(overlaps)
class CmtyAvgNodeOverlap(Statter):
    """Average memberships per node in community

    For each community: sum(memberships per node) / n_nodes"""
    ylabel = "avg memberships per node"
    log_y = True
    def calc(self, g, cmtys, cache=None):
        nodecmtys = cache_get(cache, 'nodecmtys', lambda: cmtys.nodecmtys())
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            total_memb = sum(len(nodecmtys[n]) for n in cnodes)
            x = total_memb / float(n_cmty)

            n_cmty = log_bin(n_cmty)
            yield n_cmty, x
class CmtyNodeOverlap(Statter):
    """Average memberships per node in community"""
    ylabel = "memberships per node"
    log_y = True
    def calc(self, g, cmtys, cache=None):
        nodecmtys = cache_get(cache, 'nodecmtys', lambda: cmtys.nodecmtys())
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            n_cmty = log_bin(n_cmty)
            for n in cnodes:
                yield n_cmty, len(nodecmtys[n])
class NodeOverlapByDeg(Statter):
    """Average memberships per node by degree"""
    ylabel = "memberships per node"
    xlabel = "node degree"
    log_y = True
    def calc(self, g, cmtys, cache=None):
        nodecmtys = cache_get(cache, 'nodecmtys', lambda: cmtys.nodecmtys())
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            for n in cnodes:
                yield log_bin(g.degree(n)), len(nodecmtys[n])


class MaxDegNodeJaccard(Statter):
    """jaccard(community, neighbors of highest degree node)"""
    ylabel = "cmty jacc with hub node neighbors"
    n = 1
    legend_loc = 'lower right'
    def calc(self, g, cmtys, cache=None):
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

class CmtyEmbeddednessVsSLD(Statter):
    xlabel = "SLD"
    ylabel = "cmty embeddedness"
    log_x = False
    legend_loc = "upper left"
    def calc(self, g, cmtys, cache=None):
        adj = g.adj
        cmtygraph = cmtys.cmty_graph(g)
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            n_edges_subgraph = 0
            for n in cnodes:
                for nbr in adj[n]:
                    if nbr in cnodes:
                        n_edges_subgraph += 1
            assert n_edges_subgraph % 2 == 0
            n_edges_subgraph /= 2

            # compute sld.
            sld = 2 * n_edges_subgraph / float(n_cmty-1)

            # compute embeddedness:

            k_in = cmtygraph.node[cname]['weight'] * 2
            k_outs = [ dat['weight']
                       for neigh,dat in cmtygraph[cname].iteritems() ]
            cmtyembeddedness = k_in/float(k_in+sum(k_outs))

            assert abs(sld - (k_in/float(n_cmty-1))) < 1e-6,  (sld, k_in/float(n_cmty-1))

            yield lin_bin(sld), cmtyembeddedness


class CmtyCompSize(Statter):
    """Size distribution of all components of communities.

    Note: consider using CmtyCompLCCSize instead, it gives more useful
    answers."""
    ylabel = "cmty comp size fraction"
    legend_loc = "upper right"
    log_y = True
    def calc(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            sg = g.subgraph(cnodes)
            ccs = networkx.connected_components(sg)
            #from fitz import interactnow
            ccs_sizes = [len(x) for x in ccs]
            if not all(x==1 for x in ccs_sizes):
                print ccs_sizes

            n_cmty_binned = log_bin(n_cmty)
            for size in ccs_sizes:
                if size == 1: continue
                yield n_cmty_binned, size/float(n_cmty)

class CmtyLCCSize(Statter):
    """Size of the largest connected component within communities."""
    ylabel = "cmty LCC size fraction"
    legend_loc = "upper right"
    log_y = False
    _which_comp = 0
    def calc(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            sg = g.subgraph(cnodes)
            ccs = networkx.connected_components(sg)

            n_cmty_binned = log_bin(n_cmty)
            # If we don't have the n-th largest component requested,
            # then return zero.  This makes sense since in this case,
            # the size of that component has tended to zero.
            if len(ccs) <= self._which_comp:
                yield n_cmty_binned, 0.0
                continue
            yield n_cmty_binned, len(ccs[self._which_comp])/float(n_cmty)
class CmtySLCCSize(CmtyLCCSize):
    """Second largest component size."""
    ylabel = "cmty second LCC size fraction"
    _which_comp = 1  # second largest component
class CmtySLCCRatio(Statter):
    """Size ratio of second to first largest connected comps."""
    ylabel = "cmty SLCC/LCC size ratio"
    legend_loc = "upper right"
    log_y = False
    _which_comp = 0
    def calc(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            sg = g.subgraph(cnodes)
            ccs = networkx.connected_components(sg)

            n_cmty_binned = log_bin(n_cmty)
            # If we don't have the n-th largest component requested,
            # then return zero.  This makes sense since in this case,
            # the size of that component has tended to zero.
            if len(ccs) <= self._which_comp+1:
                yield n_cmty_binned, 0.0
                continue
            yield n_cmty_binned, len(ccs[self._which_comp+1])/float(len(ccs[self._which_comp]))




class CmtyLinkComm(Statter):
    """Fraction of external links that are part of some community.

    This analyzer was made to support link community analysis work."""
    ylabel = "external edges internal fraction"
    #legend_loc = "upper right"
    _which_comp = 0
    def calc(self, g, cmtys, cache=None):
        adj = g.adj
        nodecmtys = cache_get(cache, 'nodecmtys', lambda: cmtys.nodecmtys())
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            k_ext = 0
            k_ext_int = 0
            for n in cnodes:
                for neigh in adj[n]:
                    if neigh in cnodes: continue
                    if neigh not in nodecmtys: continue
                    k_ext += 1
                    if len(nodecmtys[n]&nodecmtys[neigh]):
                        k_ext_int += 1

            n_cmty_binned = log_bin(n_cmty)
            if k_ext + k_ext_int == 0: continue
            yield n_cmty_binned, k_ext_int/float(k_ext)


class CmtyAvgShortestPath(Statter):
    """Average shortest path within community."""
    ylabel = "avg shortest path"
    # damping: a way to specify how many nodes to check for larger
    # communities.
    damping = 0.35
    def calc(self, g, cmtys, cache=None):
        ssspl = networkx.single_source_shortest_path_length

        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            # Complete or partial shortest path calculation?
            if n_cmty <= 100:
                startnodes = cnodes
            else:
                startnodes = random.sample(list(cnodes),
                               int(n_cmty*10**(-self.damping*(log10(n_cmty)-2))))

            sg = g.subgraph(cnodes)  # subgraph
            sum_l = 0.0              # sum of path lenghts accumulator
            num_l = 0.0              # number of paths accumulator
            #print n_cmty, len(startnodes)
            for node in startnodes:
                path_lengths = ssspl(sg, node)
                # Returns one path ength for every reachable node.  If
                # this is less than n_cmty, the cmty is disconnected.
                if len(path_lengths) < n_cmty:
                    # What should we do when the community is not
                    # connected?  We should never pass in communities
                    # that are unconnected, this problem needs to be
                    # solved at a higher level.  If it ever gets down
                    # to this point, we have a problem and should just
                    # abort.
                    raise RuntimeError("Cmty not connected:"
                        " cname=%s, comps=%s (%s, %s)"%(cname,
                        [len(x) for x in networkx.connected_components(sg)],
                                                        sum_l, num_l))

                sum_l += sum(path_lengths.values())
                num_l += len(path_lengths) - 1  # don't include path to self
            if num_l == 0:
                # No paths found
                print sum_l, num_l, len(cnodes)
                print [len(x) for x in networkx.connected_components(sg) ]
                raise RuntimeError("Singleton community?")
                continue
            # Compute average
            avg = sum_l / float(num_l) #float(len(startnodes)*(n_cmty-1))

            #print avg
            n_cmty_binned = log_bin(n_cmty)
            yield n_cmty_binned, avg
    def calc_igraph(self, g, cmtys, cache=None):
        for cname, cnodes in cmtys.iteritems():
            n_cmty = len(cnodes)
            if n_cmty < self.minsize:
                continue

            # Complete or partial shortest path calculation?
            if n_cmty <= 100:
                startnodes = range(n_cmty)
            else:
                startnodes = random.sample(range(n_cmty),
                               int(n_cmty*10**(-self.damping*(log10(n_cmty)-2))))

            sg = g.subgraph(cnodes)  # subgraph

            all_sp = sg.shortest_paths(source=startnodes)
            # Compute average
            avg = numpy.sum(all_sp) / float((n_cmty-1) * len(startnodes))
            n_cmty_binned = log_bin(n_cmty)
            yield n_cmty_binned, avg


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
    if obj in  (Statter, QuantileStatter, DistStatter, CountStatter,
                CcdfStatter, CccfStatter): continue
    globals()[name+'Base'] = obj
    if issubclass(obj, Statter) and not issubclass(obj, QuantileStatter):
        globals()[name+''] = type(name+'', (obj, ), {})
    if issubclass(obj, Statter) and not issubclass(obj, QuantileStatter):
        globals()[name+'Qtl'] = type(name+'Qtl', (obj,QuantileStatter), {})
    if issubclass(obj, Statter) and not issubclass(obj, DistStatter):
        globals()[name+'Dist'] = type(name+'Dist', (obj,DistStatter), {})
    if issubclass(obj, Statter) and not issubclass(obj, CountStatter):
        globals()[name+'Count'] = type(name+'Count', (obj,CountStatter), {})
    if issubclass(obj, Statter) and not issubclass(obj, CcdfStatter):
        globals()[name+'Ccdf'] = type(name+'Ccdf', (obj,CcdfStatter), {})
    if issubclass(obj, Statter) and not issubclass(obj, CccfStatter):
        globals()[name+'Cccf'] = type(name+'Cccf', (obj,CccfStatter), {})
