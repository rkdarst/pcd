# Richard Darst, July 2011

import ast
import math
from math import log, exp, floor, ceil
import numpy
import os
import random
import subprocess

import networkx

import cmodels

floateq = lambda x,y: abs(x-y) <= 5e-6*max(abs(x),abs(y))

def distance(p1, p2, boxsize=None):
    d = numpy.subtract(p1, p2)
    if not boxsize is None:
        numpy.subtract(d, numpy.round(d/boxsize)*boxsize, d)
    numpy.square(d, d)
    d = numpy.sum(d, axis=-1)
    numpy.sqrt(d, d)
    return d

def cmtyIntersect(G1, c1, G2, c2):
    """Number of nodes in intersect of c1 and c2."""
    return cmodels.cmtyIntersect(G1._struct_p, c1, G2._struct_p, c2)
def cmtyUnion(G1, c1, G2, c2):
    """Number of nodes in union of c1 and c2."""
    return cmodels.cmtyUnion(G1._struct_p, c1, G2._struct_p, c2)


class AttrToDict(object):
    def __init__(self, obj, values=None):
        self.obj = obj
        self.values = values
    def __getitem__(self, attrname):        return getattr(self.obj, attrname)
    def __setitem__(self, attrname, value): setattr(self.obj, attrname, value)
    def __iter__(self):
        if self.values is None:
            raise ValueError("AttrToDict does not have list of value names.")
        for v in self.values:
            yield v
class DictToAttr(object):
    def __init__(self, obj):                self.obj = obj
    def __getattr__(self, attrname):        return self.obj[attrname]
    def __setattr__(self, attrname, value): self.obj[attrname] = value

# This class is copied from fitz.mathutil
class Averager(object):
    """Numerically Stable Averager

    Calculate averages and standard deviations in a numerically stable way.

    From the 'On-Line Algorithm' from
    http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

    The optional initialization argument datatype= should be a
    generator function which returns zeros of the type being averaged.
    For example, to average numpy arrays of ten values, use:
      Averager(datatype=lambda: numpy.zeros(10))
    """
    def __init__(self, datatype=float):
        self._n      = 0
        self._mean   = datatype()   # mean
        self._M2     = datatype()   # variance accumulator
    def __setstate__(self, state):
        if 'n' in state:
            state['_n'] = state['n']
            del state['n']
        self.__dict__.update(state)
    def add(self, value):
        """Add a new number to the dataset.
        """
        if isinstance(value, Averager):
            # Add a sub-averager
            return self._add_child(value)
        n = self.n + 1
        delta = value - self._mean
        mean = self._mean + delta/n
        M2 = self._M2 + delta*(value - mean)
        self._n = n
        self._mean = mean
        self._M2 = M2
        return self
    __iadd__ = add

    def _add_child(self, child):
        if hasattr(self, "_child"):
            self._child.add(child)
        else:
            self._child = child
        return self

    @property
    def n(self):
        """Number of points"""
        if hasattr(self, "_child"):
            return self._n + self._child.n
        return self._n
    @property
    def mean(self):
        """Mean"""
        if hasattr(self, "_child"):
            #if self._n > 1e3 and self._child.n > 1e3:
            #    # Large value one
            #    return (self._n*self._mean + self._child.n*self._child.mean) /\
            #           (self._n + self._child.n)

            delta = self._child.mean - self._mean
            return self._mean + delta * self._child.n/(self._n+self._child.n)
        return self._mean
    @property
    def M2(self):
        """M2 algorithm parameter"""
        if hasattr(self, "_child"):
            delta = self._child.mean - self._mean
            return self._M2 + self._child.M2 + \
                   delta**2 * (self._n*self._child.n)/(self._n+self._child.n)
        return self._M2
    @property
    def std(self):
        """Population standard deviation"""
        if self.n == 0: return float('nan')
        return math.sqrt(self.M2 / self.n)
    @property
    def stdsample(self):
        """Sample standard deviation"""
        if self.n <= 1: return float('nan')
        return math.sqrt(self.M2 / (self.n-1))
    @property
    def stdmean(self):
        """(Sample) standard deviation of the mean.  std / sqrt(n)"""
        if self.n <= 1: return float('nan')
        return math.sqrt(self.M2 / ((self.n-1)*self.n))
    @property
    def var(self):
        """Population variance"""
        if self.n == 0: return float('nan')
        return self.M2 / self.n
    @property
    def varsample(self):
        """Sample variance"""
        if self.n <= 1: return float('nan')
        return self.M2 / (self.n-1)
    def proxy(self, expr):
        """Proxy for computing other"""
        return _AvgProxy(main=self, expr=expr)
    def proxygen(self, expr):
        return lambda: self.proxy(expr=expr)
class _AvgProxy(object):
    def __init__(self, main, expr):
        self.main = main
        self.expr = expr
    def add(self, x):
        raise "You can't add values to proxy objects"
    @property
    def mean(self):
        return eval(self.expr, dict(o=self.main))

class AutoAverager(object):
    def __init__(self, datatype=float, newAverager=Averager,
                 depth=1, auto_proxy=[]):
        """
        if depth=1, then self[name] will be an Averager
        if depth=2, then self[name] will be AutoAverager...
                  ...and self[name2] will be an Averager
        and so on.

        newAverager is what the leaf averager object will be created
        with.  This should be replaced with AutoAverager.

        This will automatically make hierarchical averagers
        """
        self.datatype = datatype
        self.depth = depth
        self.newAverager = newAverager
        self.names = [ ]
        self.data = { }
        self.data_list = [ ]
        self.auto_proxy = auto_proxy
    def __getitem__(self, name, newAverager=None, do_proxy=True):
        # newAverager is used for proxy objects.
        if name not in self.data:
            if self.depth > 1:
                new = self.__class__(datatype=self.datatype,
                                     newAverager=self.newAverager,
                                     depth=self.depth-1,
                                     auto_proxy=self.auto_proxy,
                                     )
            else:
                if newAverager is not None:
                    new = newAverager()
                else:
                    new = self.newAverager(datatype=self.datatype)
            self.names.append(name)
            self.data[name] = new
            self.data_list.append(new)
            # Add auto-proxies (e.g. computing std dev?)
            if self.depth == 1 and not isinstance(new, _AvgProxy) and do_proxy:
                for suffix, expr in self.auto_proxy:
                    self.add_proxy(name, name+suffix, expr=expr)
        return self.data[name]
    get = __getitem__
    def add(self, name, val, do_proxy=True):
        self.get(name, do_proxy=do_proxy).add(val)
    def remove(self, name_or_index):
        if isinstance(name_or_index, int):
            name = self.names[-1]
            del self.data[name]
            del self.names[-1]
        else:
            del self.data[name]
            self.names.remove(name)
    def __iter__(self):
        for key in self.names:
            return (key, self.data[key])
    def add_proxy(self, name_orig, name_new, expr):
        """Add a proxy object.

        Example usage to add a standard deviation column:
        .add_proxy('q', 'q_std', 'o.stdmean')"""
        if name_new not in self.data:
            self.get(name_new, newAverager=self[name_orig].proxygen(expr))
    def column(self, name, attrname='mean'):
        assert self.depth == 2
        return [getattr(self.data[name_][name], attrname) for name_ in self.names]
    def column_names(self):
        return self.data[self.names[0]].names
    def table(self, attrname='mean'):
        assert self.depth == 2
        return zip(*(self.column(name, attrname=attrname)
                     for name in self.column_names()))
class _DerivProxy(_AvgProxy):
    def __init__(self, main, name, t):
        self.main = main
    def add(self, x):
        raise "You can't add values to proxy objects"
    @property
    def mean(self):
        idx = main.names.index(t)
        if idx == 0 or idx == len(main.names)-1:
            return float('nan')
        return (main[t+1][name].mean-main[t-1][name].mean)/2.

# This class is copied from fitz.loginterval
class LogInterval(object):
    """
    interval - how many
    """
    # How many
    interval = 10
    def __init__(self, low=None, high=None,
                 interval=10, density=10, offset=1,
                 ):
        self.interval = interval
        self.density = density
        self.offset = offset
        self.expConstant = exp(log(interval) / density)
        if low:
            self._indexlow  = int(floor(round(self.index(low), 10)))
        if high:
            self._indexhigh = int(ceil(round(self.index(high), 10)))
    def value(self, index):
        return self.offset * self.expConstant**index
    def index(self, value):
        return log(value/self.offset) / (log(self.interval)/self.density)
    def indexes(self):
        return tuple(range(self._indexlow, self._indexhigh+1))
    def values(self):
        return tuple(self.value(i) for i in self.indexes())

    def __iter__(self, maxValue=None):
        # We have an upper bound
        if hasattr(self, '_indexhigh'):
            for v in self.values():
                yield v
            return

        # We have no upper bound, iterate forever.
        index = self._indexlow
        while True:
            value = self.value(index)
            yield value
            index += 1
            if maxValue is not None and value > maxValue:
                break

log2 = lambda x: log(x, 2)

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
    VI = G0.entropy_python + G1.entropy_python - 2*I
    return VI
def In(G0, G1):
    """Normalized Mutual Information"""
    I = mutual_information(G0, G1)
    In = 2.*I / (G0.entropy_python + G1.entropy_python)
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

def mutual_information_overlap(G0, G1, weighted=False):
    """LF overlap-including normalized mutual information."""
    N = 1 - .5 * (HX_Ynorm(G0, G1, weighted) + HX_Ynorm(G1, G0, weighted) )
    return N
N = mutual_information_overlap

def matrix_swap_basis(array, a, b):
    """Swap rows a,b and columns a,b in a array."""
    array[a,:], array[b,:] = array[b,:].copy(), array[a,:].copy()
    array[:,a], array[:,b] = array[:,b].copy(), array[:,a].copy()

def networkx_from_matrix(array, ignore_diagonal=True, ignore_values=()):
    """Convert a interactions matrix to a networkx graph.

    The values in the matrix are used as the graph edge weights.

    `ignore_diagonals` (default True) means to not add graph edges
    pointing back at the same node (self-loops).

    `ignore_values` can be used to ignore certain values in the
    matrix, making them not become edges.  This can be useful to
    exclude the default weights for non-interacting nodes.
    """
    assert len(array.shape) == 2
    assert array.shape[0] == array.shape[1]
    assert (array == array.T).all()
    N = array.shape[0]
    g = networkx.Graph()
    g.add_nodes_from(range(N))
    for i, row in enumerate(array):
        for j, weight in enumerate(row):
            if ignore_diagonal and i == j:
                continue
            if weight in ignore_values:
                continue
            g.add_edge(i, j, weight=weight)
    return g



def color_graph(g):
    """Color a networkX graph.

    This is rather naive right now.
    """
    #For each node,
    colors = set()
    nodes = list(g.nodes())
    while nodes:
        n = random.choice(nodes)
        nodes.remove(n)
        neighbor_colors = set(g.node[n2].get('color', None) for n2 in
                              g.neighbors(n)
                              if n2 != n)
        avail_colors = colors - neighbor_colors
        if avail_colors:
            color = random.choice(list(avail_colors))
        else:
            color = len(colors)
            colors |= set((color,))
        g.node[n]['color'] = color
    return colors

class ColorMapper(object):
    """Helper class for color mapping
    """
    def __init__(self, g, colormap_name='gist_rainbow',
                 colorizer=color_graph):
        """Initialize the system, set up normalization

        g: networkx graph.  In common usage in pcd, the nodes of this
        graph will represent communities (G.supernode_networkx).

        colormap_name: matplotlib colormap name

        colorizer: function to colorize networkx graph g
        """
        # Import only if actually used here.
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colors = colorizer(g)
        self.g = g
        self.colors = colors
        self.colormap = cm.get_cmap(colormap_name)
        self.normmap = mcolors.Normalize(vmin=0, vmax=len(colors))
    def __len__(self):
        return len(self.colors)
    def __call__(self, c):
        """Return matplotlib RGB color."""
        return self.colormap(self.normmap(self.color_id(c)))
    __getitem__ = __call__
    def color_id(self, c):
        """Return an integer representing the color of th given node."""
        return self.g.node[c]['color']

def wrap_coords(coords, boxsize):
    """Wrap coordinates to [0, boxsize)"""
    #return coords - boxsize*numpy.floor(coords/boxsize)
    return coords % boxsize
def wrap_dists(coords, boxsize):
    """Wrap coordinates to [-boxsize/2, boxsize/2)"""
    return coords - boxsize*numpy.floor((coords/boxsize)+.5)

def mean_position(coords, boxsize):
    r0 = coords[0]
    dists = wrap_dists(coords-r0, boxsize)
    mean = dists.mean(axis=0)
    mean += r0
    mean %= boxsize
    return mean

def mean_position2(coords, boxsize):
    """Handles wrapped particles"""
    coords = coords[list(nodelist)]
    return coords.mean(axis=0)
#def uniform_image(coords, boxsize):

def check_sparse_symmetric(G):
    for i, row in enumerate(G.simatrixId):
        for idx, j in enumerate(row):
            weight = G.simatrix[i, idx]
            assert i in G.simatraxId[j, :G.simatrixN[j]]
            invindex = numpy.where(G.simatraxId[j, :G.simatrixN[j]]==i)[0][0]
            difference = G.simatrix[j]

def leval(s):
    try:                              return ast.literal_eval(s)
    except (ValueError,SyntaxError):  return s

def args_to_dict(l):
    kwargs = { }
    for arg in l:
        assert '=' in arg
        k, v = arg.split('=', 1)
        kwargs[k] = leval(v)
    return kwargs

def eval_gamma_str(gammas):
    """Parse a gamma string of the form 'low,high[,density]'.

    gammas is evaluated using ast.literal_eval, so you can use any
    (and must) use python literal syntax such as .1,10 or (.1,10,20).
    In this case, the return is a dict
    dict(low=LOW,high=HIGH[,density=DENSITY]).

    If a single integer or float is given, then only this value is
    returned as an integer."""
    gammas = leval(gammas)
    if isinstance(gammas, (tuple,list)):
        if len(gammas) == 2:
            gammas = dict(low=gammas[0], high=gammas[1])
        elif len(gammas) == 3:
            gammas = dict(low=gammas[0], high=gammas[1], density=gammas[2])
        gammas['start'] = .01
    return gammas

def extremawhere(a, mask, func=numpy.max):
    """Returns location of extrema.

    The mask is *not* applied when returning the location of the mask.
    So, to find the maximum value, you must do a[where], not
    a[mask][where].
    """
    extrema = func(a[mask])

    def approxequal(a,b):
        # True if a,b are within 1%
        return a==b or abs(a-b)/float(max(abs(a),abs(b))) < 1e-2
        return numpy.less_equal(
            numpy.divide((a-b)/float((numpy.add(numpy.abs(a), numpy.abs(b))))),
            1e-2)

    #extrema_where = numpy.where(a == extrema)[0]
    #extrema_where = numpy.where(approxequal(a, extrema))[0]
    extrema_where = [i for (i,x) in enumerate(a)
                     if approxequal(x, extrema)]
    plateaus = [ [ ] ]
    for i in extrema_where:
        if not mask[i]:
            continue
        if len(plateaus[-1]) != 0 and i != plateaus[-1][-1] + 1:
            plateaus.append( [] )
        plateaus[-1].append(i)

    longest_plateau = max(plateaus, key=lambda x: len(x))
    w = longest_plateau[len(longest_plateau)//2]
    return w

def extrema(a, mask, func=numpy.max):
    """Return extrema of an array, subject to mask.

    This function is only considers where 'mask' is true.  It returns
    the middle of the longest continuous run of whatever the extrema
    is.

    mask - only consider these values as being used for the maximum.

    """
    return a[extremawhere(a, mask, func)]


def fraction_defined(G0, G1, ill=False):
    n_nodes = 0
    n_detected = 0
    n_welldef = 0
    n_welldef_detected = 0
    mapping = { }
    for c0 in G0.cmtys():
        maxOverlap = 0
        bestCmty = None
        for c1 in G1.cmtys():
            overlap = cmtyIntersect(G0,c0, G1,c1)
            if overlap > maxOverlap:
                bestCmty = c1
                maxOverlap = overlap
        mapping[c0] = bestCmty
        c1 = bestCmty
        # Considering all nodes: fraction detected
        n_detected += cmtyIntersect(G0,c0, G1,c1)
        # Considering only well-defined nodes:
        if ill:
            #print 'ill defined: ', len([ n for n in G0._graph.nodes() if G0._graph.node[n]['ill'] ])
            #print 'well defined: ', len([ n for n in G0._graph.nodes() if not G0._graph.node[n]['ill'] ])
            plantedNodes = set(G0.cmtyContents(c0))
            wellDefinedPlantedNodes = set(n for n in plantedNodes
                              if not G0._graph.node[G0._nodeLabel[n]]['ill'])
            n_welldef += len(wellDefinedPlantedNodes)
            n_welldef_detected += len(wellDefinedPlantedNodes
                                      & set(G1.cmtyContents(c1)))
    if ill:
        return (
            float(n_detected) / G0.N,
            float(n_welldef_detected) / n_welldef if n_welldef else float('nan'),
            )
    else:
        raise NotImplementedError

def hasrealattr(obj, name):
    print obj.name
    return hasattr(obj, name) and not hasattr(getattr(obj,name), '__get__')
def listAttributes(obj, depth=1, exclude=[], exclude_deep=[]):
    """Return dict of trivial attributes.

    Given an object, return a dict attr:value of all int,long,float,string attributes.  If depth>0, allow tuples"""
    attrs = { }
    #excludeNames = set(('__init__',))
    for name in dir(obj):
        if name in exclude: continue
        value = getattr(obj, name)
        if value in exclude_deep:
            if _isTrivialType(value, depth=0):
                attrs[name] = value
            else:
                attrs[name] = value.__class__
        #if name in excludeNames: continue
        if name.startswith('__'): continue
        if _isTrivialType(value, depth=depth):
            attrs[name] = value
    return attrs
_NoneType = type(None)
def _isTrivialType(t, depth=0):
    if isinstance(t, (int,long,float,str,bool,_NoneType)):
        return True
    if isinstance(t, (tuple,list)) and depth > 0:
        return all(_isTrivialType(x, depth=depth-1) for x in t)
    if isinstance(t, dict) and depth > 0:
        return all(_isTrivialType(x, depth=depth-1) for x in t.iterkeys()) \
               and \
               all(_isTrivialType(x, depth=depth-1) for x in t.itervalues())
    return False

import contextlib
import shutil
import tempfile
@contextlib.contextmanager
def chdir_context(dirname):
    """Context manager for chdir.

    chdir to 'dirname'.  Upon cleanup, chdir to the former working
    directory."""
    olddir = os.getcwd()
    os.chdir(dirname)
    yield
    os.chdir(olddir)

@contextlib.contextmanager
def tmpdir_context(chdir=False, delete=True, suffix='', prefix='tmp', dir=None):
    """Context manager for temporary directories.

    Create a temporary directory using context manager.  If 'chdir' is
    true (default false), then chdir to the directory and to the
    original directory upon cleanup.  If delete is true (default
    true), then delete the temporary directory and all contents upon
    cleanup.  Other arguments (suffix, prefix, dir) are passed to
    tempfile.mkdtemp.

    Theb context argument is the tmpdir name, absolute path."""
    olddir = os.getcwd()
    tmpdir = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
    if chdir:
        os.chdir(tmpdir)
    yield tmpdir
    if chdir:
        os.chdir(olddir)
    if delete:
        shutil.rmtree(tmpdir)

def NMI_LF(cmtys1, cmtys2, check=True, use_existing=False):
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

    with tmpdir_context(chdir=True, dir='.', prefix='tmp-nmi-'):
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
ovIn_LF = NMI_LF

#
# The following class WeightedChoice is used for selecting items from
# a group, with a corresponding weight.
#

# This doesn't exist until 3.2, so redefine here:
import operator
def _accumulate(iterable, func=operator.add):
    'Return running totals'
    # accumulate([1,2,3,4,5]) --> 1 3 6 10 15
    # accumulate([1,2,3,4,5], operator.mul) --> 1 2 6 24 120
    it = iter(iterable)
    total = next(it)
    yield total
    for element in it:
        total = func(total, element)
        yield total

import bisect
class WeightedChoice(object):
    """Class to pick among n items, weighted.

    Initialize class as
    itemsweights = [('a', 5), ('b', 4), ('c', 1)]
    WeightedChoice(itemsweights)

    and then each call to .choice() will return 'a' 50% of the time,
    'b' 40% of the time, and 'c' 10% of the time.  Weights can be any
    numbers with satisfy normal addition and comparison operations.

    rng: instance of random.Random, default built-in RNG
        Use this as the random number generator, instead of the
        automatically created random.Random instance (found via
        random.random.__self__).
    """
    def __init__(self, itemsweights, rng=random.random.__self__):
        self.rng = rng
        items, weights = zip(*itemsweights)
        self.cumdist = list(_accumulate(weights))
        self.norm = self.cumdist[-1]
        self.items = list(items)
    def choice(self):
        """Pick one item according to the weights."""
        x = self.rng.random() * self.norm
        return self.items[bisect.bisect(self.cumdist, x)]
    def add(self, item, weight):
        """Add a given item, with """
        self.cumdist.append(self.cumdist[-1]+weight)
        self.norm = self.cumdist[-1]
        self.items.append(item)
    def __len__(self):
        """Number of items in the chooser"""
        assert len(self.cumdist) == len(self.items)
        return len(self.cumdist)


def mainfunc(ns, defaultmethod='run'):
    import sys
    mode = sys.argv[1]
    kwargs = { }
    if len(sys.argv) > 2:
        if '=' in sys.argv[2]:
            kwargs = args_to_dict(sys.argv[2:])
        else:
            defaultmethod = sys.argv[2]
            if len(sys.argv) > 3:
                kwargs = args_to_dict(sys.argv[3:])
    getattr(ns[mode], defaultmethod)(**kwargs)


def graph_stats(g, prefix='', recurse=1, level=1, _test=False):
    """Compute graph statistics.

    recurse: int
        if true, recurse and compute stats of largest connected
        component, too.
    level: int, default 1
        0: run only O(n) or O(m) tests
        1: also run triangle tests
        2: also calculate diameter
    _test: bool, default False
        Run internal assertions to verify this functions's triangle
        calculations match networkx's. (Some test logic is reproduced
        for efficiency reasons.)
    prefix: str
        internal use, prefix output lines with this.  Do not use.

    Returns: list
        List of strings containing human-readable description of
        graph.  To print it out, use '\n'.join(lines).
    """
    stats = [ ]
    #print "nodes"
    N = g.number_of_nodes()
    stats.append("Nodes: %d"%N)
    #print "edges"
    E = g.number_of_edges()
    stats.append("Edges: %d"%E)
    #print "density"
    stats.append("Density: %f"%(E/float(N*(N-1))))
    stats.append("Mean-Degree: %f"%(2*E/float(N)))
    if not g.is_directed():
        # Undirected graphs:
        if level >= 1:
            ##print "triangles"
            #stats.append("Triangles: %d"%(sum(networkx.triangles(g).itervalues())//3))
            ##print "transitivity"
            #stats.append("Transitivity: %f"%networkx.transitivity(g))
            ##print "acc"
            #stats.append("Average-Clustering-Coefficient: %f"%networkx.average_clustering(g))

            # Implement the above three functions myself.  This saves from
            # extra internal iterations of the
            # cluster._triangles_and_degree_iter function.
            import networkx.algorithms.cluster as cluster
            triangles2 = 0
            possible_triangles2 = 0
            clustering_coef_sum = 0.0
            n_nodes = 0
            # The following function DOUBLE-COUNTS triangles
            for node, d, t in cluster._triangles_and_degree_iter(g):
                # Everywhere here, we have double-counts for t and
                # possible_triangles.
                triangles2 += t
                possible_triangles2 += d*(d-1)
                count_zeros = True
                if t > 0:
                    # exclude nodes with zero triangles from ACC.  This
                    # includes nodes that have no triangles possible.
                    n_nodes += 1
                    clustering_coef_sum += t/float(d*(d-1))
            if triangles2 == 0:   transitivity = 0.0
            else:                 transitivity = triangles2/float(possible_triangles2)
            avg_clustering_coefficient = clustering_coef_sum / float(g.number_of_nodes())
            if n_nodes == 0:
                avg_clustering_coefficient_nonzero = float('nan')
            else:
                avg_clustering_coefficient_nonzero = clustering_coef_sum / float(n_nodes)
            #print n_nodes, clustering_coef_sum
            if _test:
                assert triangles2//6 == sum(networkx.triangles(g).itervalues())//3
                assert abs(transitivity - networkx.transitivity(g)) < 1e-6
                assert abs(avg_clustering_coefficient-networkx.average_clustering(g)) < 1e-6

            stats.append("Triangles: %d"%(triangles2//6))
            #print "transitivity"
            stats.append("Transitivity: %f"%transitivity)
            #print "acc"
            stats.append("Average-Clustering-Coefficient: %f"%avg_clustering_coefficient)
            stats.append("Average-Clustering-Coefficient-Nonzero: %f"%avg_clustering_coefficient_nonzero)

        if level >= 1.5 and transitivity != 1.0:
            # This function raises a warning when run on complete
            # graphs:
            # networkx/algorithms/assortativity/correlation.py:285:
            # RuntimeWarning: invalid value encountered in
            # double_scalars
            stats.append("Degree-Assortativity-Coefficient: %s"%
                         networkx.degree_assortativity_coefficient(g))

        if level >= 1.7:
            comps = networkx.connected_components(g)
            if len(comps)==1 and level >= 2:
                #print "diameter"
                stats.append("Diameter: %d"%networkx.diameter(g))
                #print "radius"
                stats.append("Radius: %d"%networkx.radius(g))
            stats.append("Is-Connected: %d"%(len(comps)==1))
            if len(comps) > 1:
                stats.append("Conn-Components-Number-Of: %d"%len(comps))
                stats.append("Conn-Components-Fraction-Nodes-In-Largest: %g"%(
                    len(comps[0])/float(len(g))))
                stats.append("Conn-Components-Sizes-Top-20: %s"%(
                    " ".join(str(len(x)) for x in comps[:20])))
                stats.append("Conn-Components-Size-Fractions-Top-5: %s"%(
                    " ".join("%g"%(len(x)/float(len(g))) for x in comps[:5])))
                stats.append("Conn-Components-Number-Singletons: %d"%(
                                 sum(1 for x in comps if len(x)==1)))
                # Fails for directed graphs...
                if recurse:
                    #print "recursing"
                    lcc = g.subgraph(comps[0])
                    stats.extend(graph_stats(lcc, prefix='LCC-', recurse=recurse-1))

    if g.is_directed() and level >= 1.7:
        #print "ncc: s/w"
        s_comps = networkx.strongly_connected_components(g)
        stats.append("Is-Strongly-Connected: %d"%(len(s_comps) == 1))
        w_comps = networkx.weakly_connected_components(g)
        stats.append("Is-Weakly-Connected: %d"%(len(w_comps) == 1))
        #
        if len(s_comps) > 1:
            stats.append("Strong-Conn-Components-Number-Of: %d"%len(s_comps))
            stats.append("Strong-Conn-Components-Fraction-Nodes-In-Largest: %g"%(
                             len(s_comps[0])/float(len(g))))
            stats.append("Strong-Conn-Components-Sizes-Top-20: %s"%(
                             " ".join(str(len(x)) for x in s_comps[:20])))
            stats.append("Strong-Conn-Components-Size-Fractions-Top-5: %s"%(
                             " ".join("%g"%(len(x)/float(len(g))) for x in s_comps[:5])))
            stats.append("Strong-Conn-Components-Number-Singletons: %d"%(
                             sum(1 for x in s_comps if len(x)==1)))
            if recurse:
                lcc = g.subgraph(s_comps[0])
                stats.extend(graph_stats(lcc, prefix='LSCC-', recurse=recurse-1))
        if len(w_comps) > 1:
            stats.append("Weak-Conn-Components-Number-Of: %d"%len(w_comps))
            stats.append("Weak-Conn-Components-Fraction-Nodes-In-Largest: %g"%(
                             len(w_comps[0])/float(len(g))))
            stats.append("Weak-Conn-Components-Sizes-Top-20: %s"%(
                             " ".join(str(len(x)) for x in w_comps[:20])))
            stats.append("Weak-Conn-Components-Size-Fractions-Top-5: %s"%(
                             " ".join("%g"%(len(x)/float(len(g))) for x in w_comps[:5])))
            stats.append("Weak-Conn-Components-Number-Singletons: %d"%(
                             sum(1 for x in w_comps if len(x)==1)))
            if recurse:
                lcc = g.subgraph(w_comps[0])
                stats.extend(graph_stats(lcc, prefix='LWCC-', recurse=recurse-1))

    #stats.append(": %d"%networkx.(g))
    if prefix:
        stats = [prefix+line for line in stats]
    return stats

class Proxy(object):
    def __init__(self, gen):
        self.__gen = gen
        self.__obj = None
    def __get(self):
        if self.__obj is None:
            self.__obj = self.__gen()
            print self.__obj
        return self.__obj
    def __getattr__(self, name):
        return getattr(self.__get(), name)
    def __len__(self):
        return self.__get().__len__()



from sqlite import SQLiteDict

if __name__ == "__main__":
    # tests
    print logTime(-10)
    print logTime(logTimeIndex(100))
