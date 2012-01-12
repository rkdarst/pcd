# Richard Darst, July 2011

from taco.mathutil.stats import Averager
from math import log, exp, floor, ceil
import numpy
import random

import networkx

import cmodels

def distance(p1, p2, boxsize=None):
    d = numpy.subtract(p1, p2)
    if not boxsize is None:
        numpy.subtract(d, numpy.round(d/boxsize)*boxsize, d)
    numpy.square(d, d)
    d = numpy.sum(d, axis=-1)
    numpy.sqrt(d, d)
    return d

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
        """Population Variance"""
        if self.n == 0: return float('nan')
        return math.sqrt(self.M2 / self.n)
    @property
    def stdsample(self):
        """Sample Variance"""
        if self.n <= 1: return float('nan')
        return math.sqrt(self.M2 / (self.n-1))
    @property
    def var(self):
        """Population Standard Deviation"""
        if self.n == 0: return float('nan')
        return self.M2 / self.n
    @property
    def varsample(self):
        """Sample Standard Deviation"""
        if self.n <= 1: return float('nan')
        return self.M2 / (self.n-1)


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
    """Calculate mutual information between two graphs."""
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
            s0 = set(G0.cmtyll[c0, :n0])
            s1 = set(G1.cmtyll[c1, :n1])
            n_shared = len(s0 & s1)
            #assert n_shared == len(s0 & s1)

            if n_shared == 0:
                continue

            MI += (n_shared/float(N)) * log2(n_shared*N/float(n0*n1))

    return MI
def mutual_information_c(G0, G1):
    return cmodels.mutual_information(G0._struct_p, G1._struct_p)
mutual_information = mutual_information_c




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
def HX_Ynorm_c(GX, GY):
    return cmodels.HX_Ynorm(GX._struct_p, GY._struct_p)
HX_Ynorm = HX_Ynorm_c

def mutual_information_overlap(G0, G1):
    N = 1 - .5 * (HX_Ynorm(G0, G1) + HX_Ynorm(G1, G0) )
    return N


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


if __name__ == "__main__":
    # tests
    print logTime(-10)
    print logTime(logTimeIndex(100))
