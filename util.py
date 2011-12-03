# Richard Darst, July 2011

from taco.mathutil.stats import Averager
from math import log, exp, floor, ceil
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
            self._indexlow  = int(floor(self.index(low)))
        if high:
            self._indexhigh = int(ceil(self.index(high)))
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


if __name__ == "__main__":
    # tests
    print logTime(-10)
    print logTime(logTimeIndex(100))
