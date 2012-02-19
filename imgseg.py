# Richard Darst, January 2012

import collections
import colorsys
import math
import numpy
#import numpy.fft
import scipy.misc
import scipy.ndimage

from pcd import Graph, MultiResolution
import pcd.util

vrgb_to_hsv = numpy.vectorize(colorsys.rgb_to_hsv)
vhsv_to_rgb = numpy.vectorize(colorsys.hsv_to_rgb)

def imread(fname):
    """Load R,G,B (in [0,1]) channels of image.

    The shape of R,G,B is (Ly,Lx).  When printed, x is the horizontal
    axis and y is the vertical axis."""
    a = scipy.misc.imread(fname)
    if len(a.shape) == 2:
        return a/255.  # for grayscale images, return only main channel
    # If RGBA, return only RGB.
    R, G, B = a[...,0], a[...,1], a[...,2]
    R, G, B = R/255., G/255., B/255.
    return R, G, B

def imsave(fname, bands):
    """Make image (filename) out of R,G,B in [0,1]"""
    #R,G,B = vhsv_to_rgb(H, S, V)
    scipy.misc.imsave(fname, numpy.asarray(bands)*255)

def replace_hue(R, G, B, H):
    shape = R.shape
    oldH,S,V = vrgb_to_hsv(R,G,B)
    #def _h(x,y):
    #        return x/shape[0]
    #H = numpy.fromfunction(_h, shape)
    S = .7  # intensity of the color
    V = numpy.clip(V, a_min=.3, a_max=.7)  # remove pure white and black
    return vhsv_to_rgb(H,S,V)
def replace_channels(R, G, B, channels):
    shape = R.shape
    H,S,V = vrgb_to_hsv(R,G,B)
    if 'H' in channels:
        H = channels['H']
    if 'S' in channels:
        S = channels['S']
        if isinstance(S, tuple) and S[0] == 'clip':
            S = numpy.clip(S,a_min=S[1],a_max=S[2])
    if 'V' in channels:
        V = channels['V']
        if isinstance(S, tuple) and V[0] == 'clip':
            V = numpy.clip(V,a_min=V[1],a_max=V[2])
        #V = numpy.clip(channels['V'],a_min=.9,a_max=.9)
    return vhsv_to_rgb(H,S,V)
def dict_to_array(d, shape=None, keys=None):
    yLow = min(k[0] for k in d.keys())
    xLow = min(k[1] for k in d.keys())
    Ly = max(k[0] for k in d.keys())-yLow+1
    Lx = max(k[1] for k in d.keys())-xLow+1
    if shape is None:
        shape = shape=(Ly, Lx)
    array = numpy.zeros(shape=(Ly, Lx))
    for node in interactions.keys():
        array[tuple(numpy.subtract(node, (yLow,xLow)))] = d[node]
    return array



class Image(object):
    def __init__(self, channels, mode=None):
        self.channels = channels

        if isinstance(channels, tuple):
            R,G,B = channels
            self.R = R
            self.G = G
            self.B = B
            self.mode = 'RGB'
            self.shape = self.R.shape
        else:
            self.L = channels
            self.mode = 'L'
            self.shape = self.L.shape
    def crop(self, crop):
        """Crop image inplace.

        Argument is slice syntax, for example '10:90,20:50' to do
        array[10:90,20:50]"""
        for channel in self.mode:
            chan = getattr(self, channel)
            chan = eval('chan[%s]'%crop)
            setattr(self, channel, chan)
            self.shape = chan.shape
    def resize(self, scale):
        """Resize image to 'scale' fraction of previous size.

        If scale is:
          int: percent of old size
          float: fraction of old size
          tuple: new dimensions
        """
        # Note: see scipy.misc.imresize for docs on how this could be
        # made more flexible than just a fraction.
        for channel in self.mode:
            chan = getattr(self, channel)
            chan = scipy.misc.imresize(chan*255, scale, )/255. #interp='cubic'
            setattr(self, channel, chan)
            self.shape = chan.shape
    def save(self, fname):
        if isinstance(fname, str):
            fname = open(fname, 'wb')
        if self.mode == 'RGB':
            imsave((self.R, self.G, self.B))
        elif self.mode == 'L':
            imsave(fname, self.L)
        else:
            raise TypeError('unknown mode: %s'%self.mode)

    def get(self):
        if self.mode == 'RGB':
            return self.HSV[2]
        elif self.mode == 'L':
            return self.L

    @property
    def HSV(self):
        H,S,V = vrgb_to_hsv(self.R, self.G, self.B)
        return H,S,V

    @property
    def avg(self):
        if self.mode == 'L':
            return self.L
        elif self.mode == 'RGB':
            return (self.R + self.G + self.B) / 3.


    def plotCmtys(self, fname, G, I):

        # Make and print community map.
        cmtys = I._newarr(dtype=int)
        for node in I.iterCoords():
            cmtys[node] = G.cmty[G._nodeIndex[node]]
        #print "plotCmtys:", cmtys.max()
        if cmtys.max() > 0:
            cmtys = cmtys / float(cmtys.max())

        if self.mode == 'RGB':
            R,G,B = replace_hue(self.R, self.G, self.B, cmtys)
            imsave(fname, (R,G,B))
        elif self.mode == 'L':
            S = I._newarr() ; S[:] = .7
            V = I._newarr() ; V = self.L
            V = numpy.clip(V, a_min=.3, a_max=.7) # remove pure white and black
            imsave(fname, vhsv_to_rgb(cmtys, S, V))


class ImgSeg(object):
    def __init__(self, img, width=5, corrlength=5, L=2.0,
                 mode='intensity',
                 defaultweight=0,
                 #cutoff=None,
                 Vbar=0,
                 vartop=False,
                 dist=None):
        self.Ly, self.Lx = img.shape
        self.img = img

        if numpy.max(img) - numpy.min(img) < .5:
            raise Warning("Not much dynamic range in the image...")

        self.L = float(L)
        self.defaultweight = defaultweight
        #gamma = .05
        self.width = width
        self.corrlength = corrlength
        #self.cutoff = cutoff
        self.Vbar = Vbar

        self._inner_dtype = numpy.float32

        if mode == 'A':
            self.width = 0
            self.inner_func = self.inner_intensity
            self.outer_func = self.outer_
        elif mode == 'intensity':
            # Maximum difference: 0(attractive) -- 1(repulsive),
            # adjusted range: -1(attr) -- 0
            self.inner_func = self.inner_intensity
            self.outer_func = self.outer_intensityavg
            self.dist_func  = self.dist_cutoff
            #self.dist_func  = self.expdecay
        elif mode == 'intensity-cutoff':
            # Range: -1(attr) -- 0
            self.inner_func = self.inner_intensity
            self.outer_func = self.inner_intensitycutoff
            self.dist_func  = self.dist_expdecay
        elif mode == 'frequency':
            # range: -x (overlap, attr) -- 0(no overlap, repl)
            self.inner_func = self.inner_FFT_mean
            self.outer_func = self.outer_FFT_conj
            self.dist_func = self.dist_expdecay
            self._inner_dtype = numpy.complex64
        elif mode == 'frequencyb':
            # range: -x (overlap, attr) -- 0(no overlap, repl)
            self.inner_func = self.inner_FFT_paper
            self.outer_func = self.outer_FFT_conj
            self.dist_func = self.dist_expdecay
            self._inner_dtype = numpy.complex64
        elif mode == 'frequency-abs':
            self.inner_func = self.inner_FFT_mean
            self.outer_func = self.outer_FFT_absdiff
            self.dist_func = self.dist_expdecay
            self._inner_dtype = numpy.complex64
        else:
            raise ValueError("Unknown analysis mode: %s"%mode)

        if dist is None:
            pass
        elif dist is 'cutoff':
            self.dist_func  = self.dist_cutoff
        elif dist is 'exp':
            self.dist_func  = self.dist_expdecay
        else:
            raise ValueError("Unknown distance mode: %s"%dist)

    def _newarr(self, extrashape=(), dtype=numpy.float32):
        return numpy.zeros(shape=(self.Ly, self.Lx)+extrashape,
                           dtype=dtype)

    def dist_expdecay(self, d):
        return math.exp(-d/self.L)
    def dist_cutoff(self, d):
        return 1 if d <= self.L else 0

    # All of these outer-functions should have lower = more attractive.

    # Method 'intensity-cutoff': intensity, cutoff.  Set width=0 or
    # else it will use a mean.
    def inner_intensity(self, a):
        return a
    def outer_intensitycutoff(self, a1, a2):
        """-1 if mean of differences is < self.cutoff, else 0"""
        return -(numpy.mean(numpy.abs(a1-a2)) < self.cutoff).astype(numpy.int_)
    # Method 'intensity': Intensity, continuous.  Set width=0 or else it will
    # use a mean of the local block.
    def outer_intensityavg(self, a1, a2):
        """Mean of absolute values of differences."""
        return numpy.mean(numpy.abs(a1-a2)) - 1

    # Method 2: average intensity differences between blocks
    @staticmethod
    def outer_intensityblock(a1, a2):
        # Positive function.
        return numpy.abs(numpy.mean(a1 - a2))

    # Method 3: frequency domair
    @staticmethod
    def inner_FFT_mean(a):
        return numpy.fft.fft2(a - numpy.mean(a))
    #_inner_dtype = numpy.complex64
    @staticmethod
    def outer_FFT_conj(a1, a2):
        """This is the method in the paper."""
        return - numpy.mean(numpy.abs((a1.conj() * a2).flat))
    # Method 3b: frequency domain
    @staticmethod
    def inner_FFT_paper(a):
        return numpy.fft.fft2(a - a.flat[0])
    #_inner_dtype = numpy.complex64
    @staticmethod
    def outer_FFT_absdiff(a1, a2):
        """FFT - difference of absolute values"""
        return numpy.mean(numpy.abs(numpy.abs(a1)-numpy.abs(a2)))


    #@staticmethod
    #def outer_FFT_absdiff(a1, a2):
    #    return numpy.mean(numpy.abs(numpy.abs(a1)-numpy.abs(a2)))

    def Gmask(self, G):
        return (numpy.isnan(G.simatrix) == False)

    def adjust_weights(self, V='mean', G=None, default=None):
        print "Adjusting weights"
        if G is None:
            G = self.G()

        mask = self.Gmask(G)
        print "before:", \
              "min:", G.simatrix[mask].min(), \
              "max:", G.simatrix[mask].max(), \
              "mean:", G.simatrix[mask].mean(), \
              "std:", G.simatrix[mask].std(), \
              "rmatrixDefault:", G.srmatrixDefault
        if V == 'mean':
            V = numpy.mean(G.simatrix[mask])
        elif V == 'median':
            V = numpy.median(G.simatrix[mask])
        elif V == 'max':
            V = numpy.max(G.simatrix[mask])
        elif V == 'half':
            max_ = numpy.max(G.simatrix[mask])
            min_ = numpy.min(G.simatrix[mask])
            V = min_ + .5*(max_ - min_)
        elif '%' in V:
            import re
            percent = float(re.match(r'([-0-9.]+)%', V).group(1))/100.
            max_ = numpy.max(G.simatrix[mask])
            min_ = numpy.min(G.simatrix[mask])
            V = min_ + percent*(max_ - min_)
        #print meanV
        #assert meanV < 0
        G.simatrix[mask] -= V
        #G.srmatrixDefault -= meanV/2.
        if default == 'max':
            G.srmatrixDefault = numpy.max(G.simatrix[mask])
        elif default is not None:
            G.srmatrixDefault = default
        else:
            G.srmatrixDefault -= V
        print "after:", \
              "min:", G.simatrix[mask].min(), \
              "max:", G.simatrix[mask].max(), \
              "mean:", G.simatrix[mask].mean(), \
              "std:", G.simatrix[mask].std(), \
              "rmatrixDefault:", G.srmatrixDefault


    #def make_local0(self):
    #    local = self._newarr(extrashape=(2*self.width+1,
    #                                     2*self.width+1))
    #    width = self.width
    #    inner_func = self.inner_func
    #
    #    i = [0]
    #    pieces = { }
    #    def _f(a):
    #        a.shape = 2*self.width+1, 2*self.width+1
    #        #print a
    #        #print i[0]
    #        n = i[0]
    #        x, y = n//self.Ly, n%self.Ly
    #        local[(y, x)] = inner_func(a)
    #        #print coords
    #        #print a[2,2], array[coords]
    #        i[0] += 1
    #        return i[0]-1
    #    scipy.ndimage.filters.generic_filter(self.img, _f,
    #                             size=(2*self.width+1, 2*self.width+1),
    #                             origin=0,
    #                             mode='constant', cval=-1)
    #    self.local = local
    #    return local
    def make_blocks(self):
        self.blocks = self._newarr(extrashape=(2*self.width+1,
                                              2*self.width+1),
                            dtype=getattr(self, '_inner_dtype'))
        self.mean_blocks = self._newarr()
        width = self.width
        for y, x in self.iterCoords():
            l = self.img[y-width:y+width+1,
                         x-width:x+width+1]
            # center = l[width-1:width+1+1,width-1:width+1+1]
            a = self.inner_func(l)
            self.mean_blocks[y,x] = numpy.mean(numpy.abs(a))
            self.blocks[y,x] = a
            #if (y,x) == (5,5) or (y,x) == (44,44):
            #    from fitz.interact import interact ; interact()

    @property
    def N(self):
        return (self.Ly-width*2) * (self.Lx-width*2)
    def iterCoords(self):
        for y in range(self.width, self.Ly-self.width):
            for x in range(self.width, self.Lx-self.width):
                yield y, x
    def iterAdjCoords(self, y, x):
        for adjy in range(max(y-self.corrlength, self.width),
                          min(y+self.corrlength+1, self.Ly-self.width-1)):
            for adjx in range(max(x-self.corrlength, self.width),
                              min(x+self.corrlength+1, self.Lx-self.width-1)):
                if y==adjy and x==adjx:
                    continue
                yield adjy, adjx
    @property
    def mask(self):
        return (slice(self.width, self.Ly-self.width),
                slice(self.width, self.Lx-self.width))


    def iterweights(self):
        d_map = collections.defaultdict(list)
        self.mean_conn = self._newarr()

        for y, x in self.iterCoords():
            node_conn = [ ]
            for adjy, adjx in self.iterAdjCoords(y, x):
                d = math.sqrt((y-adjy)**2+(x-adjx)**2)
                #print x, y, adjx, adjy
                a1 = self.blocks[y,x]
                a2 = self.blocks[adjy,adjx]
                #print fft1, fft2
                weight = self.outer_func(a1, a2)
                weight = weight * self.dist_func(d)
                #assert weight <= 0
                #print "%5.2f % 5.2f"%(d, weight)
                d_map[d].append(weight)
                #if abs(adjx - x) <= 1 and abs(adjy - y) <= 1:
                node_conn.append(weight)

                #offy, offx = adjy-y+self.corrlength, adjx-+self.corrlength
                yield (y, x), (adjy, adjx), weight

                #self.weights[y][x][offy][offx] = weight
            self.mean_conn[(y,x)] = numpy.mean(node_conn)
        #for k,v in sorted(d_map.iteritems()):
        #    print k, numpy.mean(v), numpy.std(v)

    #def mean_weight(self):
    #    numpy.mean([numpy.mean(connections.values())
    #                for connections in self.weights.values()])
    #print mean_interaction

    def enableVT(self):
        G = self.G()
        G.enableVT2()
        print G.srmatrixDefaultOnlyDefined
        print G.srmatrixDefault
        #import time ; time.sleep(100)

    def G(self):
        if hasattr(self, "_G"):
            return self._G
        maxconn = (2*self.corrlength +1)**2
        self._G = Graph.from_sparseiter(nodes=self.iterCoords(),
                                        weights=self.iterweights(),
                                        default=self.defaultweight,
                                        maxconn=maxconn
                                        )
        #print G.imatrix
        return self._G

    def do_stats(self, basename):
        print "Information:"
        print "inner", self.inner_func
        print "outer", self.outer_func
        print "dist", self.dist_func
        print "dtype", self._inner_dtype
        print "Information on mean_blocks (inner blocks):"
        mask = self.mask
        # Make and print average interaction map
        mean_blocks = self.mean_blocks.copy()[mask]
        print "Blocks stats:", mean_blocks.min(), mean_blocks.max(), \
                              mean_blocks.mean(), mean_blocks.std()
        mean_blocks -= mean_blocks.min()
        mean_blocks /= mean_blocks.max()
        mean_blocks2 = self.mean_blocks.copy()
        mean_blocks2[mask] = mean_blocks
        mean_blocks = mean_blocks2
        print imsave(basename+'-blocksmap.png',
                     mean_blocks)
                 #replace_channels(R,G,B, channels=dict(H=0, S=0,
                 #                                      V=mean_blocks)))

        # Print mean_conn.
        print "Information on mean_conn"
        mean_conn = self.mean_conn.copy()
        print "Connection stats:",mean_conn[mask].min(),mean_conn[mask].max(),\
                                  mean_conn[mask].mean(), mean_conn[mask].std()
        mean_conn[mask] -= mean_conn[mask].min()
        mean_conn[mask] /= mean_conn[mask].max()
        print imsave(basename+'-connmap.png',
                     mean_conn)

        print "G stats:"
        G = self.G()
        mask = self.Gmask(G)
        print "imatrix stats:", G.simatrix[mask].min(),G.simatrix[mask].max(),\
                                G.simatrix[mask].mean(),G.simatrix[mask].std()
        print "simatrixDefault, srmatrixDefault:", G.simatrixDefault,\
              G.srmatrixDefault

    def do_connmap(self, orient=None):
        self._newarr()

    def do_gamma(self, fname, gamma):
        print "Doing minimizaition"
        # Do a minimization
        G = self.G()
        G.trials(gamma, 5, maxrounds=25, minimizer='greedy2')
        print "Plotting cmtys."
        fullimg.plotCmtys(fname, G, self)
        #sys.exit()

    def do(self, basename, gammas=None, replicas=3, trials=2):
        """Do a full MR."""
        print "Making IS blocks neighborhoods"

        print "Making G object"
        G = self.G()

        # Do full mulitresolution
        print "Doing MultiResolution"
        def callback(gamma, data, **kwargs):
            G = data['Gmin']
            #fullimg.plotCmtys(basename+'-mr_gamma%011.7f.png'%gamma, G, self)
            left, right = 10, 5
            fname = basename+'-mr_gamma%%0%d.%df.png'%(left+right+1,right)
            fullimg.plotCmtys(fname%gamma, G, self)
        MR = MultiResolution(overlap=False,
                             minimizerargs=dict(trials=trials,
                                                maxrounds=10,
                                                minimizer='greedy2'),
                             output=basename+'-mrvalues.txt',
                             )
        MR.run(Gs=[G]*replicas,
               gammas=gammas,
               threads=6,
               callback=callback
               )
        MR.plot(basename+'-MR.png')
        MR.write(basename+'-MR.txt')


if __name__ == "__main__":
    import sys
    import optparse
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-g", "--gammas", help="gL,gH[,D] or a dict.")
    parser.add_option("-w", "--width", type=int, default=0)
    parser.add_option("-c", "--corrlength", type=int, default=1)
    parser.add_option("-L", "--length", type=float, default=2.0, help="Length for exp decay or distance cutoff (%default).")
    parser.add_option("-m", "--mode", default='intensity')
    parser.add_option("--dist", help="Override distance function mode")
    parser.add_option("--VT", default=False, action='store_true', help="Use variable topology potts model")
    parser.add_option("-s", "--shift")
    parser.add_option("-T", "--trials", type=int, default=2)
    parser.add_option("-R", "--replicas", type=int, default=3)
    parser.add_option("--crop", type=str, help="Python array slice to crop with, example: 10:50,10:50")
    parser.add_option("--resize", type=float, help="Resize array to this fraction of previous size.  Applied after --crop."  )
    options, args = parser.parse_args()



    fname = args[0]
    if len(args) >= 2:
        basename = args[1]
    else:
        basename = args[0]

    fullimg = Image(imread(fname))

    print "Image shape:", fullimg.shape
    if options.crop:
        print "Cropping",
        fullimg.crop(options.crop)
        print "new shape", fullimg.shape
    if options.resize:
        print "Resizing",
        fullimg.resize(options.resize)
        print "new shape", fullimg.shape

    img = fullimg.get()

    print "Making IS object"
    I = ImgSeg(img, mode=options.mode,
               width=options.width, corrlength=options.corrlength,
               #width=3, corrlength=3,
               L=options.length,
               defaultweight=0,
               #cutoff=.06,
               dist=options.dist
               )
    #I.dist_func = I.dist_expdecay

    I.make_blocks()
    G = I.G()
    #from fitz import interactnow

    I.do_stats(basename=basename)

    if options.VT:
        I.enableVT()

    if options.shift:
        default = None
        shift = pcd.util.leval(options.shift)
        if isinstance(shift, (list,tuple)):
            shift, default = shift
        elif ',' in shift:
            shift, default = shift.split(',')
            shift = pcd.util.leval(shift)
            default = pcd.util.leval(default)

        I.adjust_weights(shift, default=default)
        #I.adjust_weights('max')
        #I.adjust_weights(.06, default=0)
    #I.adjust_weights('half', default=0)

    I.do_stats(basename=basename)

    #sys.exit(2)
    #from fitz import interactnow
    print '>'
    import time
    time.sleep(10)

    if options.gammas:
        gammas = pcd.util.leval(options.gammas)
        if isinstance(gammas, (tuple,list)):
            if len(gammas) == 2:
                gammas = dict(low=gammas[0], high=gammas[1])
            elif len(gammas) == 3:
                gammas = dict(low=gammas[0], high=gammas[1], density=gammas[2])
    else:
        gammas = dict(low=.0000001, high=1, density=2)

    if isinstance(gammas, (float, int)):
        #sys.exit()
        I.do_gamma(basename+'_gamma%09.5f.png'%gammas, gammas)
        #sys.exit()
    else:
        I.do(basename=basename, gammas=gammas, trials=options.trials,
             replicas=options.replicas)

    #I.do(basename=out,
    #     gammas=dict(low=.01, high=.01, density=2))

