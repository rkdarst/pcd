# Richard Darst, January 2012
# -*- encoding: utf-8 -*-

import collections
import colorsys
import math
import numpy
#import numpy.fft
import scipy.misc
import scipy.ndimage
import sys
import time

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
    if a.shape[2] == 4:
        R, G, B, A = a[...,0], a[...,1], a[...,2], a[...,3]
        R, G, B, A = R/255., G/255., B/255., A/255.
        return R, G, B, A
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

        if mode not in (None, 'RGB', 'HSV', 'L', 'RGBA'):
            raise ValueError("Unknown mode")
        if mode == 'HSV':
            channels = vhsv_to_rgb(*channels)
        if mode == 'RGBA' or \
               isinstance(channels, tuple) and len(channels) == 4:
            self.mode = 'RGBA'
            self.R, self.G, self.B, self.A = channels
            self.shape = self.R.shape
            return
        elif isinstance(channels, tuple):
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
    def copy(self, mode=None):
        if mode and mode != self.mode:
            if mode == 'RGB':
                if self.mode == 'L':
                    return Image((self.L, self.L, self.L))
        if self.mode == 'RGB':
            return Image((self.R.copy(), self.G.copy(), self.B.copy()))
        if self.mode == 'RGBA':
            return Image((self.R.copy(), self.G.copy(), self.B.copy(), self.A.copy()))
        elif self.mode == 'L':
            return Image(self.L.copy())
        else:
            raise ValueError("Unknown mode.")
    def RGBchannels(self):
        if self.mode == 'RGB' or self.mode == 'RGBA':
            return numpy.asarray((self.R, self.G, self.B))
        elif self.mode == 'L':
            return numpy.asarray((self.L, self.L, self.L))
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
          '=NN': scale largest dimension to this size
        """
        # Note: see scipy.misc.imresize for docs on how this could be
        # made more flexible than just a fraction.
        if isinstance(scale, str) and scale.startswith("="):
            if scale.endswith('h'):
                scale = int(scale[1:-1]) / float(self.shape[0])
            elif scale.endswith('w'):
                scale = int(scale[1:-1]) / float(self.shape[1])
            else:
                scale = int(scale[1:])/float(max(self.shape[0],self.shape[1]))
        for channel in self.mode:
            chan = getattr(self, channel)
            chan = scipy.misc.imresize(chan*255, scale, )/255. #interp='cubic'
            setattr(self, channel, chan)
            self.shape = chan.shape
    def compose(self, other):
        """Compost other on top of self.  Other's alpha channel is used.

        Other must have an alpha channel."""
        self.R = self.R*(1-other.A) + other.R*other.A
        self.G = self.G*(1-other.A) + other.G*other.A
        self.B = self.B*(1-other.A) + other.B*other.A
    def save(self, fname):
        if isinstance(fname, str):
            fname = open(fname, 'wb')
        if self.mode == 'RGB':
            imsave(fname, (self.R, self.G, self.B))
        elif self.mode == 'RGBA':
            imsave(fname, (self.R, self.G, self.B, self.A))
        elif self.mode == 'L':
            imsave(fname, self.L)
        else:
            raise TypeError('unknown mode: %s'%self.mode)
    @classmethod
    def read(cls, fname):
        return cls(imread(fname))
    def get(self):
        if self.mode =='RGB' or self.mode == 'RGBA':
            return self.HSV[2]
        elif self.mode == 'L':
            return self.L
        else:
            raise RuntimeError('self.mode is unkwown: %s'%self.mode)

    @property
    def HSV(self):
        H,S,V = vrgb_to_hsv(self.R, self.G, self.B)
        return H,S,V

    @property
    def avg(self):
        if self.mode == 'L':
            return self.L
        elif self.mode == 'RGB' or self.mode == 'RGBA':
            return (self.R + self.G + self.B) / 3.


    def getCmtyMap(self, G, I, map_=None, overlay_img=True):

        # Make and print community map.
        cmtys = I._newarr(dtype=int)
        if G is not None:
            for node in I.iterCoords():
                cmtys[node] = G.cmty[G._nodeIndex[node]]
        elif map_ is not None:
            for k, v in map_.iteritems():
                cmtys[k] = v
        #print "plotCmtys:", cmtys.max()
        if cmtys.max() > 0:
            # the +1 is since hue(0) == hue(1)
            cmtys = cmtys / float(cmtys.max()+1)


        print cmtys
        if not overlay_img:
            # We do not want to underlay the original image...
            S = I._newarr() ; S[:] = 1
            V = I._newarr() ; V[:] = 1
            return Image(vhsv_to_rgb(cmtys, S, V))
        elif self.mode == 'RGB' or self.mode == 'RGBA':
            R,G,B = replace_hue(self.R, self.G, self.B, cmtys)
            return Image((R,G,B))
        elif self.mode == 'L':
            S = I._newarr() ; S[:] = .7
            V = I._newarr() ; V = self.L
            V = numpy.clip(V, a_min=.3, a_max=.7) # remove pure white and black
            return Image(vhsv_to_rgb(cmtys, S, V))

    def plotCmtys(self, fname, G, I, overlay=None):
        cmtyimg = self.getCmtyMap(G, I)
        if overlay:
            cmtyimg.compose(overlay)
        cmtyimg.save(fname)
    def getPixbuf(self):
        return img_to_pixbuf(self)


class ImgSeg(object):
    modes = ('A',
             'intensity', 'intensity-cutoff',
             'frequency', 'frequencyb', 'frequency-abs', 'frequency-test',
             )
    distmodes = (None, 'cutoff', 'exp', 'exp-rweights')
    def __init__(self, img, width=5, corrlength=5, L=2.0,
                 mode='intensity',
                 defaultweight=0,
                 #cutoff=None,
                 #Vbar=0,
                 dist=None,
                 VTmode=None,
                 exp_beta=1,
                 blur=0,
                 shift=None,
                 overlap=False,
                 overlay=None):
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
        #self.Vbar = Vbar
        if options.exp_beta is not None:
            self.exp_beta = exp_beta
        else:
            self.exp_beta = 1
        self.blur = blur
        self.shift = shift
        self.overlap = overlap
        self.overlay = overlay

        self.mode = mode
        self.distmode = dist
        self.VTmode = VTmode
        self.setup()

    def setup(self):
        """Setup internal data structures based on config."""
        mode = self.mode
        dist = self.distmode
        VTmode = self.VTmode
        self.rweights = False


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
            self.outer_func = self.outer_FFT_conj_mean
            self.dist_func = self.dist_expdecay
            self._inner_dtype = numpy.complex64
        elif mode == 'frequencyb':
            # range: -x (overlap, attr) -- 0(no overlap, repl)
            self.inner_func = self.inner_FFT_paper
            self.outer_func = self.outer_FFT_conj_mean
            self.dist_func = self.dist_expdecay
            self._inner_dtype = numpy.complex64
        elif mode == 'frequency-abs':
            self.inner_func = self.inner_FFT_mean_norm
            self.outer_func = self.outer_FFT_absdiff
            self.dist_func = self.dist_expdecay
            self._inner_dtype = numpy.complex64
        elif mode == 'frequency-test':
            self.inner_func = self.inner_FFT_mean_norm
            self.outer_func = self.outer_FFT_conj_sum
            self.dist_func = self.dist_expdecay
            self._inner_dtype = numpy.complex64
            if VTmode is None:
                VTmode = 'preDefinedOnly' #dict(mode='preDefined',repelValue=0)
            if dist is None:
                dist = 'exp-rweights'
        else:
            raise ValueError("Unknown analysis mode: %s"%mode)

        if dist is None:
            pass
        elif dist == 'cutoff':
            self.dist_func  = self.dist_cutoff
        elif dist == 'exp':
            self.dist_func  = self.dist_expdecay
        elif dist == 'exp-rweights':
            self.rweights = True
            self.dist_func = self.dist_expdecay
        else:
            raise ValueError("Unknown distance mode: %s"%dist)

        self._realVTmode = VTmode

    def _newarr(self, extrashape=(), dtype=numpy.float32):
        return numpy.zeros(shape=(self.Ly, self.Lx)+extrashape,
                           dtype=dtype)

    def dist_expdecay(self, d):
        return math.exp(-d/self.L)
    def dist_cutoff(self, d):
        return 1. if d <= self.L else 0
    def dist_one(self, d):
        return 1.

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

    ## Method 2: average intensity differences between blocks
    #@staticmethod
    #def outer_intensityblock(a1, a2):
    #    # Positive function.
    #    return numpy.abs(numpy.mean(a1 - a2))

    # Method 3: frequency domair
    @staticmethod
    def inner_FFT_mean(a):
        return numpy.fft.fftn(a - numpy.mean(a))
    @staticmethod
    def inner_FFT_mean_norm(a):
        """Normalized from -1 -- 0.

        -1 = perfect overlap
        0 = no overlap."""

        #a = a - a.mean()
        ##assert numpy.sum(numpy.abs(a)) > 1e-5
        #if numpy.sum(numpy.abs(a)) > 1e-5:
        #    print "warn: no"
        #    a /= math.sqrt(numpy.square(a).sum())
        #a = numpy.fft.fftn(a)
        #return a

        a_norm = a - a.mean()
        if numpy.sum(numpy.abs(a_norm)) > 1e-5:
            a = a_norm
        else:
            a = a.copy()
        if numpy.sum(numpy.abs(a)) > 1e-5:
            a /= math.sqrt(numpy.square(a).sum() * a.size)
        else:
            a[:] = 1./a.size
        a = numpy.fft.fftn(a)


        #a_norm = a - a.mean()
        #if numpy.sum(numpy.abs(a_norm)) > 1e-5:
        #    a = a_norm
        #else:
        #    a = a.copy()
        #a = numpy.fft.fftn(a)
        #if False:
        #    #mask = numpy.zeros(a.shape)
        #    #mask[:int(mask.shape[0]*.5),   :int(mask.shape[0]*.5)] = 1
        #    mask = numpy.zeros(a.shape)
        #    mask[int(mask.shape[0]*0):int(mask.shape[0]*.5),
        #         int(mask.shape[1]*0):int(mask.shape[1]*.5)
        #         ] = 1
        #    #print "applying mask",int(mask.shape[0]*.25),int(mask.shape[0]*.75)
        #
        #    #mask = numpy.ones(a.shape)
        #    #mask[int(mask.shape[0]*.25):int(mask.shape[0]*.75),
        #    #     int(mask.shape[1]*.25):int(mask.shape[1]*.75)
        #    #     ] = 0
        #
        #
        #    a *= mask
        #
        #if numpy.sum(numpy.abs(a)) > 1e-5:
        #    a /= math.sqrt(numpy.square(numpy.abs(a)).sum())
        #else:
        #    pass
        #    #a[:] = 1./a.size



        return a
    #_inner_dtype = numpy.complex64
    @staticmethod
    def outer_FFT_conj_sum(a1, a2):
        """This is the method in the paper."""
        return - numpy.sum(numpy.abs((a1.conj() * a2).flat))
    @staticmethod
    def outer_FFT_conj_mean(a1, a2):
        """This is my mis-implementation of the method in the paper."""
        return - numpy.mean(numpy.abs((a1.conj() * a2).flat))
    # Method 3b: frequency domain
    @staticmethod
    def inner_FFT_paper(a):
        return numpy.fft.fftn(a - a.flat[0])
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

    def do_shift(self, shift):
        # Common arguments: "mean", "median", "max", "half", "NN%"
        print "Shifting weights..."
        default = None
        shift = pcd.util.leval(shift)
        if isinstance(shift, (list,tuple)):
            shift, default = shift
        elif ',' in shift:
            shift, default = shift.split(',')
            shift = pcd.util.leval(shift)
            default = pcd.util.leval(default)

        self.adjust_weights(shift, default=default)
        #self.adjust_weights('max')
        #self.adjust_weights(.06, default=0)
        #self.adjust_weights('half', default=0)

        self.do_stats(basename=basename)

    def adjust_weights(self, V='mean', G=None, default=None):
        print "Adjusting weights"
        if G is None:
            G = self.G()

        mask = self.Gmask(G)
        print "    before:", \
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
        print "    after:", \
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
    def blocks(self):
        if hasattr(self, '_blocks'):
            return self._blocks
        blocks = self._newarr(extrashape=(2*self.width+1,
                                              2*self.width+1),
                            dtype=getattr(self, '_inner_dtype'))
        self._blocks = blocks
        self.mean_blocks = self._newarr()
        width = self.width
        for y, x in self.iterCoords():
            l = self.img[y-width:y+width+1,
                         x-width:x+width+1]
            # center = l[width-1:width+1+1,width-1:width+1+1]
            a = self.inner_func(l)
            if self.blur:
                if numpy.iscomplexobj(a):
                    ar = scipy.ndimage.filters.gaussian_filter(
                        a.real, self.blur, mode='wrap')
                    ai = scipy.ndimage.filters.gaussian_filter(
                        a.imag, self.blur, mode='wrap')
                    a.real = ar
                    a.imag = ai
                else:
                    a = scipy.ndimage.filters.gaussian_filter(
                        a, self.blur, mode='wrap')
                #print a
            self.mean_blocks[y,x] = numpy.mean(numpy.abs(a))
            blocks[y,x] = a
            #if (y,x) == (5,5) or (y,x) == (44,44):
            #    from fitz.interact import interact ; interact()
        return blocks

    @property
    def N(self):
        return (self.Ly-width*2) * (self.Lx-width*2)
    def iterCoords(self):
        for y in range(self.width, self.Ly-self.width):
            for x in range(self.width, self.Lx-self.width):
                yield y, x
    def iterAdjCoords(self, y, x):
        width = self.width
        corrlength = self.corrlength
        Ly, Lx = self.Ly, self.Lx
        for adjy in range(max(y-corrlength, width),
                          min(y+corrlength+1, Ly-width-1)):
            for adjx in range(max(x-corrlength, width),
                              min(x+corrlength+1, Lx-width-1)):
                if y==adjy and x==adjx:
                    continue
                yield adjy, adjx
    @property
    def mask(self):
        """Mask returning the image values within the active range."""
        return (slice(self.width, self.Ly-self.width),
                slice(self.width, self.Lx-self.width))


    def iterweights(self):
        """Iterate over values suitable for Graph.from_sparseiter.

        This is used for """
        d_map = collections.defaultdict(list)
        self.mean_conn = mean_conn = self._newarr()

        blocks = self.blocks()
        outer_func = self.outer_func
        dist_func = self.dist_func
        iterAdjCoords = self.iterAdjCoords

        for coords in self.iterCoords():
            y, x = coords
            print "\r%4s, %4s"%(y,x),
            sys.stdout.flush()
            node_conn = [ ]
            for adjcoords in iterAdjCoords(y, x):
                adjy, adjx = adjcoords
                #d = math.sqrt((y-adjy)**2+(x-adjx)**2)
                d = numpy.subtract(coords, adjcoords)
                numpy.multiply(d, d, d)
                d = numpy.sqrt(d.sum())
                #print x, y, adjx, adjy
                a1 = blocks[y,x]
                a2 = blocks[adjy,adjx]
                #print a1, a2
                fweight = outer_func(a1, a2)
                dist_weight = dist_func(d)
                weight = fweight * dist_weight
                #assert weight <= 0
                #print "%5.2f % 5.2f"%(d, weight)
                d_map[d].append(weight)
                #if abs(adjx - x) <= 1 and abs(adjy - y) <= 1:
                node_conn.append(weight)

                #offy, offx = adjy-y+self.corrlength, adjx-+self.corrlength
                #print (y,x), (adjy, adjx), weight, dist_weight
                if not self.rweights:
                    yield (y, x), (adjy, adjx), weight
                else:
                    #print (y,x), (adjy, adjx), fweight, dist_weight
                    yield (y, x), (adjy, adjx), weight, dist_weight**self.exp_beta

                #self.weights[y][x][offy][offx] = weight
            mean_conn[(y,x)] = numpy.mean(node_conn)
        #for k,v in sorted(d_map.iteritems()):
        #    print k, numpy.mean(v), numpy.std(v)
        print

    #def mean_weight(self):
    #    numpy.mean([numpy.mean(connections.values())
    #                for connections in self.weights.values()])
    #print mean_interaction

    #def enableVT(self, **kwargs):
    #    G = self.G()
    #    if not kwargs:
    #        kwargs = getattr(self, "_VTargs", dict(mode='onlyDefined'))
    #    G.enableVT(**kwargs)

    def G(self):
        if hasattr(self, "_G"):
            return self._G
        maxconn = (2*self.corrlength +1)**2
        kwargs = { }
        if self.rweights:
            kwargs['rmatrix'] = True
        self._G = Graph.from_sparseiter(nodes=self.iterCoords(),
                                        weights=self.iterweights(),
                                        default=self.defaultweight,
                                        maxconn=maxconn,
                                        **kwargs
                                        )
        if self._realVTmode is not None:
            self._G.enableVT(self._realVTmode)
        if self.shift:
            self.do_shift(self.shift)

        #print G.imatrix
        return self._G

    def do_stats(self, basename):
        print "Information on imgseg object:"
        print "    inner", self.inner_func
        print "    outer", self.outer_func
        print "    dist", self.dist_func
        print "    dtype", self._inner_dtype
        print "  Information on mean_blocks (inner blocks):"
        mask = self.mask
        # Make and print average interaction map
        mean_blocks = self.do_blockmean(mask=mask)
        #mean_blocks = self.mean_blocks.copy()[mask]
        print "    Blocks stats:", mean_blocks.min(), mean_blocks.max(), \
                                   mean_blocks.mean(), mean_blocks.std()
        #mean_blocks -= mean_blocks.min()
        #mean_blocks /= mean_blocks.max()
        #mean_blocks2 = self.mean_blocks.copy()
        #mean_blocks2[mask] = mean_blocks
        #mean_blocks = mean_blocks2
        if basename:
            imsave(basename+'-blocksmap.png', mean_blocks)

        # Print mean_conn.
        print "  Information on mean_conn"
        mean_conn = self.mean_conn.copy()
        print "    Connection stats:", \
                               mean_conn[mask].min(),mean_conn[mask].max(),\
                               mean_conn[mask].mean(), mean_conn[mask].std()
        mean_conn[mask] -= mean_conn[mask].min()
        mean_conn[mask] /= mean_conn[mask].max()
        if basename:
            imsave(basename+'-connmap.png', mean_conn)

        print "  G stats:"
        G = self.G()
        mask = self.Gmask(G)
        print "    imatrix stats:", \
                            G.simatrix[mask].min(),G.simatrix[mask].max(),\
                            G.simatrix[mask].mean(),G.simatrix[mask].std()
        if isinstance(G.srmatrix, numpy.ndarray):
            print "    rmatrix stats:", \
                            G.srmatrix[mask].min(),G.srmatrix[mask].max(),\
                            G.srmatrix[mask].mean(),G.srmatrix[mask].std()
        print "    simatrixDefault, srmatrixDefault, " \
              "srmatrixDefaultOnlyDefined:", G.simatrixDefault, \
              G.srmatrixDefault, G.srmatrixDefaultOnlyDefined
    def do_blockmean(self, mask=None):
        if mask is None:
            mask = self.mask
        blocks = self.blocks() # Generate self.mean_blocks in this function.o
        mean_blocks = self.mean_blocks.copy()[mask]
        print "    Blocks stats:", mean_blocks.min(), mean_blocks.max(), \
                                   mean_blocks.mean(), mean_blocks.std()
        mean_blocks -= mean_blocks.min()
        mean_blocks /= mean_blocks.max()
        mean_blocks2 = self.mean_blocks.copy()
        mean_blocks2[mask] = mean_blocks
        mean_blocks = mean_blocks2
        return mean_blocks

    def do_connmap(self, orient=None):
        self._newarr()
        mean_conn = self.mean_conn.copy()
        print "    Connection stats:", \
                               mean_conn[mask].min(),mean_conn[mask].max(),\
                               mean_conn[mask].mean(), mean_conn[mask].std()
        mean_conn[mask] -= mean_conn[mask].min()
        mean_conn[mask] /= mean_conn[mask].max()

    def do_overlapmap(self, basename, coords, use_distweight=False):
        """Local block overlap with respect to one coordinate.

        Given one coordinate, plot the overlap between its local block
        and all other blocks."""
        arr = self._newarr()
        y, x  = coords
        blocks = self.blocks()
        a1 = blocks[y,x]

        for adjy, adjx in self.iterCoords():
            d = math.sqrt((y-adjy)**2+(x-adjx)**2)
            a2 = blocks[adjy,adjx]

            fweight = self.outer_func(a1, a2)
            dist_weight = self.dist_func(d)
            weight = fweight * dist_weight

            #if adjx==56 and adjy==43:
            #    print a1, a2, d, fweight, dist_weight, weight
            if not use_distweight:
                arr[(adjy, adjx)] = fweight
            else:
                arr[(adjy, adjx)] = weight

        arr[self.mask] -= arr[self.mask].min()
        arr[self.mask] /= arr[self.mask].max()
        if basename is not None:
            imsave(basename+'-overlapmap_%d,%d.png'%(y,x), arr)
        else:
            return arr
    def do_saveblocks(self, basename):
        blocks = self.blocks()
        f = open(basename+'-blocks.txt', 'w')
        def fnum(x):
            if numpy.iscomplex(x):
                return '%+4.2f%+4.2fj'%(x.real, x.imag)
            else:
                return '%+4.2f'%x
        for y, x in self.iterCoords():
            block = blocks[y,x]
            size = block.size
            line = [ ]
            line.append('%04d,%04d: '%(y,x))
            line.append(
                '  '.join(
                ' '.join(
                fnum(x) for x in row) for row in block))
            line.append('\n')
            f.write("".join(line))

    def do_gamma(self, fname, gamma):
        print "Doing minimizaition"
        # Do a minimization
        G = self.G()
        G.trials(gamma, 5, maxrounds=25, minimizer='greedy2')
        print "Plotting cmtys."
        fullimg.plotCmtys(fname, G, self, self.overlay)

    def do(self, basename, gammas=None, replicas=3, trials=2,
           threads=6):
        """Do a full MR."""
        print "Making IS blocks neighborhoods"

        print "Making G object"
        G = self.G()

        # Do full mulitresolution
        print "Doing MultiResolution"
        def callback(gamma, data, **kwargs):
            G = data['Gmin']
            #fullimg.plotCmtys(basename+'-mr_gamma%011.7f.png'%gamma, G, self)
            left, right = 5, 5
            fname = basename+'-mr_gamma%%0%d.%df.png'%(left+right+1,right)
            fullimg.plotCmtys(fname%gamma, G, self, self.overlay)
        if isinstance(self.overlap, (int, float)):
            overlap = self.overlap
        elif self.overlap:
            overlap = trials
        else:
            overlap = False
        MR = MultiResolution(overlap=overlap,
                             minimizerargs=dict(trials=trials,
                                                maxrounds=25,
                                                minimizer='greedy2'),
                             output=basename+'-mrvalues.txt',
                             savestate=False,
                             )
        MR.no_N = True
        MR.run(Gs=[G]*replicas,
               gammas=gammas,
               threads=threads,
               callback=callback
               )
        MR.plot(basename+'-MR.png')
        MR.write(basename+'-MR.txt')

def img_to_pixbuf(img):
    import gtk

    channels = img.RGBchannels()
    channels = channels * 255
    channels = channels.astype(numpy.uint8)  # 3, h, w
    channels = channels.swapaxes(0,2)        # w, h, 3
    channels = channels.swapaxes(0,1)        # h, w, 3
    #pb = gtk.gdk.pixbuf_new_from_array(channels,
    #                                   gtk.gdk.COLORSPACE_RGB, 8)
    pb = gtk.gdk.pixbuf_new_from_data(
        channels.tostring(), gtk.gdk.COLORSPACE_RGB,
        has_alpha=False, bits_per_sample=8,
        width=channels.shape[1], height=channels.shape[0],
        rowstride=channels.shape[1]*3)
    return pb

def gtk_main(I, img):
    import pygtk
    #pygtk.require('2.0')
    import gtk
    import logging
    class ImgSegGTK(object):
        log = logging.getLogger('ImgSegGTK')
        log.setLevel(logging.DEBUG)
        log.addHandler(logging.StreamHandler())
        scale = 3
        coords = None
        vizmode = 'mirror'
        vizmodes = ('mirror', 'blockmean', 'overlap', 'cmtys', 'cmty_local')
        VTmodes=(None, "standard", 'onlyDefined', 'preDefined',
                 'preDefinedOnly')
        shifts = (None, "mean", "median", "max", "25%", "half", "75%")
        modes = ImgSeg.modes
        distmodes = ImgSeg.distmodes

        #
        # Interface wrapper.
        #
        def _combo_map_get(self, widget):
            widget = getattr(self, widget)
            active = widget.get_active_text()
            return active
        def _combo_map_set(self, val, widget, options):
            if val == 'None': val = None
            widget = getattr(self, widget)
            options = list(getattr(self, options))
            #print val, options
            widget.set_active(options.index(val))

        vizmode = property(
            lambda self  : self._combo_map_get('ctrl_vizmode', ),
            lambda self,x: self._combo_map_set(x, 'ctrl_vizmode', 'vizmodes'),
            )
        mode = property(
            lambda self, : self._combo_map_get('ctrl_mode', ),
            lambda self,x: self._combo_map_set(x, 'ctrl_mode', 'modes'),
            )
        distmode = property(
            lambda self, : self._combo_map_get('ctrl_distmode', ),
            lambda self,x: self._combo_map_set(x, 'ctrl_distmode','distmodes'),
            )
        VTmode = property(
            lambda self, : self._combo_map_get('ctrl_VTmode'),
            lambda self,x: self._combo_map_set(x, 'ctrl_VTmode', 'VTmodes'))
        shift = property(
            lambda self, : self._combo_map_get('ctrl_shift'),
            lambda self,x: self._combo_map_set(x, 'ctrl_shift', 'shifts'))

        #def _vizmode_get(self):
        #    active = self.ctrl_vizmode.get_active_text()
        #    return active
        #def _vizmode_set(self, vizmode):
        #    self.ctrl_vizmode.set_active(list(self.vizmodes).index(vizmode))
        #vizmode = property(_vizmode_get, _vizmode_set)

        def _entry_get(self, widget):
            widget = getattr(self, widget)
            val = pcd.util.leval(widget.get_text())
            return val

        width   = property(lambda self: self._entry_get('ctrl_width'))
        corrlen = property(lambda self: self._entry_get('ctrl_corrlen'))
        L       = property(lambda self: self._entry_get('ctrl_L'))
        beta    = property(lambda self: self._entry_get('ctrl_beta'))
        gamma   = property(lambda self: self._entry_get('ctrl_gamma'))
        trials  = property(lambda self: self._entry_get('ctrl_trials'))
        blur    = property(lambda self: self._entry_get('ctrl_blur'))


        #@property
        #def gamma(self):
        #    return pcd.util.leval(self.ctrl_gamma.get_text())
        #
        #@property
        #def trials(self):
        #    return pcd.util.leval(self.ctrl_trials.get_text())

        def update_values(self, width=None):
            """Updates the """
            self.log.debug('update_values')
            def del_if_hasattr(obj, attrname):
                if hasattr(obj, attrname):
                    #print "deleting", attrname
                    delattr(obj, attrname)

            #I = self.I

            #width = pcd.util.leval(self.ctrl_width.get_text())
            #if width is not None and width != I.width:
            #    print "New width:", width
            #    I.width = width
            #    del_if_hasattr(I, '_G')
            #    del_if_hasattr(I, '_blocks')
            #self.update_img()

            stuff_to_update = (
                ('width', self.width, ('_G', '_blocks')),
                ('corrlength', self.corrlen, ('_G', '_blocks')),
                ('L', self.L, ('_G', '_blocks')),
                ('exp_beta', self.beta, ('_G', '_blocks')),
                ('blur', self.blur, ('_G', '_blocks')),

                ('mode', self.mode, ('_G', '_blocks', 'setup')),
                ('distmode', self.distmode, ('_G', '_blocks', 'setup')),
                ('VTmode', self.VTmode, ('_G', '_blocks', 'setup')),
                ('shift', self.shift, ('_G', '_blocks', 'setup')),
                )

            needs_setup = False
            for name, val, to_del in stuff_to_update:
                if val == 'None':
                    val = None
                if val != getattr(self.I, name):
                    self.log.debug("Changed value: %s=%r"%(name, val))
                    setattr(self.I, name, val)
                    for name in to_del:
                        if name == 'setup':
                            needs_setup = True
                            continue
                        del_if_hasattr(I, name)
            if needs_setup:
                I.setup()

        def do_cd(self, dummy=None):
            self.update_values()
            if self.vizmode not in ('cmtys', ):
                self.vizmode = 'cmtys'
            G = self.I.G()
            self.I.do_stats(None)
            G.trials(gamma=self.gamma, trials=self.trials, minimizer='greedy2',
                     threads=6)
            img = self.img.getCmtyMap(G, self.I)
            self.pb_cmtys = img.getPixbuf()
            self.update_img()
        def do_cdlocal(self):
            self.update_values()
            if not hasattr(self, 'coords'):
                return None
            G = self.I.G()
            G.cmtyListClear()
            #y, x = self.coords
            G.cmty[:] = 0
            G.cmty[G._nodeIndex[tuple(self.coords)]] = 1
            G.cmtyListInit()
            #G.cmtyListAdd(0, G._nodeIndex[tuple(self.coords)])
            G.trials(gamma=self.gamma, trials=self.trials,
                     minimizer='ovexpand', initial='current',
                     threads=6)
            map_ = { }
            for coord in self.I.iterCoords():
                n = G._nodeIndex[coord]
                if G.cmtyContains(1, n):
                    map_[coord] = 1
                else:
                    map_[coord] = 0

            img = self.img.getCmtyMap(G=None, I=self.I, map_=map_)
            self.pb_cmty_local = img.getPixbuf()
            self.update_img()


        #
        # Utility
        #
        def scale_pb(self, pb):
            return pb.scale_simple(int(pb.get_width()*self.scale),
                                   int(pb.get_height()*self.scale),
                                   gtk.gdk.INTERP_NEAREST)
        # gtk.gdk.INTERP_HYPER

        #
        # Events
        #
        def event_scale(self, scale_adj):
            """Update event"""
            self.log.debug
            self.scale = scale_adj.get_value()
            self.update_view()

        def event_vizmode(self, vizmodebox):
            self.update_img()

        def event_imageclick(self, widget, event):
            self.log.debug("event_imageclick %r", event)

            if event is None:
                self.log.debug("event_imageclick -> None")
                self.update_values()
                self.update_img()
                return

            if event.button == 1:
                x, y = int(event.x), int(event.y)
                x //= self.scale ; y //= self.scale
                coords = self.coords = y, x
                self.log.info("Clicked: %s, %s", x, y)
                print "block:"
                print self.I.blocks()[y,x]
            # Different event types.
            if event.button==3 and event.type==gtk.gdk.BUTTON_RELEASE:
                self.log.debug("event_imageclick -> 3 release")
                self.pb_bottom = self._old_pb_bottom ; del self._old_pb_bottom
                self.update_view()
                return
            elif event.button == 3:
                if hasattr(self, '_old_pb_bottom'): return
                self.log.debug("event_imageclick -> 3 press")
                self._old_pb_bottom = self.pb_bottom
                self.pb_bottom = self.pb_orig
                self.update_view()
                return
            elif event.button==1 and self.vizmode == 'cmty_local':
                self.do_cdlocal()
            elif event.button==1:
                self.update_values()
                self.log.debug("event_imageclick -> other")
                self.vizmode = 'overlap'
                x, y = int(event.x), int(event.y)
                x //= self.scale ; y //= self.scale
                coords = self.coords = y, x
            self.update_img()


        def update_img(self):
            """Update the graphics for the two sides."""
            self.log.debug("update imgs")
            if self.coords:
                coords = self.coords
                y, x = coords
            else:
                coords = None
            # Top image creation:
            # draw the box
            img = self.img.copy()
            if coords is not None:
                w = I.width
                # top, bottom:
                for x_ in (x-w, x+w):
                    img.R[y-w:y+w+1, x_     ] = 0
                    img.G[y-w:y+w+1, x_     ] = 0
                    img.B[y-w:y+w+1, x_     ] = 0
                for y_ in (y-w, y+w):
                    img.R[y_     , x-w:x+w+1] = 0
                    img.G[y_     , x-w:x+w+1] = 0
                    img.B[y_     , x-w:x+w+1] = 0

            pb = img_to_pixbuf(img)
            self.pb = pb

            if self.vizmode == 'mirror':
                self.pb_bottom = pb
            elif self.vizmode == 'blockmean':
                arr = I.do_blockmean()
                pb = img_to_pixbuf(Image((arr, arr, arr)))
                self.pb_bottom = pb
            elif self.vizmode == 'overlap':
                if coords is not None:
                    arr = I.do_overlapmap(basename=None, coords=(y, x))
                    pb = img_to_pixbuf(Image((arr, arr, arr)))
                self.pb_bottom = pb
            elif self.vizmode == 'cmtys':
                if hasattr(self, 'pb_cmtys'):
                    self.pb_bottom = self.pb_cmtys
                else:
                    self.pb_bottom = pb
            elif self.vizmode == "cmty_local":
                if hasattr(self, 'pb_cmty_local'):
                    self.pb_bottom = self.pb_cmty_local
                else:
                    self.pb_bottom = pb

            self.update_view()


        def update_view(self):
            """Updates scaling/paning of display."""
            self.log.debug("update view")
            self.image.set_from_pixbuf(self.scale_pb(self.pb))
            self.image_overlap.set_from_pixbuf(self.scale_pb(self.pb_bottom))


        #
        # Beans and rice
        #
        def do_saveimg(self, button):
            """Save the bottom image to a file"""
            chooser = gtk.FileChooserDialog(
                parent=self.window,
                action=gtk.FILE_CHOOSER_ACTION_SAVE,
                buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                         gtk.STOCK_ZOOM_IN, gtk.RESPONSE_ACCEPT,
                         gtk.STOCK_OPEN, gtk.RESPONSE_OK))
            #chooser.set_default_response(gtk.RESPONSE_OK)
            response = chooser.run()
            if response == gtk.RESPONSE_OK:
                fname = chooser.get_filename()
                print fname, 'selected'
                pb = self.pb_bottom
                pb.save(fname, type=fname.split('.')[-1])
            elif response == gtk.RESPONSE_ACCEPT:
                fname = chooser.get_filename()
                print fname, 'selected'
                pb = self.scale_pb(self.pb_bottom)
                pb.save(fname, type=fname.split('.')[-1])
            chooser.destroy()


        def __init__(self, I, img):
            self.I = I
            self.img = img.copy(mode='RGB')
            self.pb_orig = self.img.getPixbuf()

            # The window
            self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
            #self.window.set_default_size(1000, 1000)
            self.window.set_size_request(750, 750)
            self.window.set_border_width(5)
            self.window.resize(min(img.shape[1]+50, 1e3),
                               min(img.shape[0]*2+100, 1e3))
            self.window.connect("destroy_event", lambda w: gtk.main_quit())
            self.window.set_title("%s - %s"%(globals().get('fname', '???'),
                                             sys.argv[0]))

            # The image
            self.image = image = gtk.Image()
            self.image.show()
            # The eventbox
            eb = gtk.EventBox()
            eb.connect('button-press-event', self.event_imageclick)
            eb.connect('button-release-event', self.event_imageclick)
            eb.add(image)
            eb.show()
            image = eb

            a = gtk.Alignment(xalign=.5, yalign=.5)
            a.add(image)
            a.show()
            image = a

            sw1 = gtk.ScrolledWindow(hadjustment=None, vadjustment=None)
            sw1.add_with_viewport(image)
            sw1.show()


            # The bottom image
            self.image_overlap = image = gtk.Image()
            self.image_overlap.show()
            # The eventbox
            eb = gtk.EventBox()
            eb.connect('button-press-event', self.event_imageclick)
            eb.connect('button-release-event', self.event_imageclick)
            eb.add(image)
            eb.show()
            image = eb

            a = gtk.Alignment(xalign=.5, yalign=.5)
            a.add(image)
            a.show()
            image = a

            sw2 = gtk.ScrolledWindow(hadjustment=sw1.get_hadjustment(),
                                     vadjustment=sw1.get_vadjustment())
            sw2.add_with_viewport(image)
            sw2.show()

            imgbox = gtk.VBox(homogeneous=False, spacing=5)
            imgbox.pack_start(sw1)
            imgbox.pack_start(sw2)
            imgbox.show()


            # The controls
            controls = gtk.VBox()
            controls.set_size_request(150, -1)
            controls.show()
            def add_label(widget, label):
                hbox = gtk.HBox()
                hbox.show()
                l = gtk.Label(label)
                l.show()
                hbox.pack_start(l, False)
                hbox.pack_start(widget, False)
                return hbox
            def make_combo(options, default=0, c=controls, label=None):
                combobox = gtk.combo_box_new_text()
                for m in options:
                    if m is None: m = 'None'
                    combobox.append_text(m)
                #combobox.set_active(default)
                combobox.show()
                if label:
                    c.pack_start(add_label(combobox, label), False)
                else:
                    c.pack_start(combobox, False)
                return combobox
            def make_entry(label, default, width=None, c=controls):
                """Given label and default value, make a text field."""
                #hbox = gtk.HBox()
                #hbox.show()
                #l = gtk.Label(label)
                #l.show()
                e = gtk.Entry()
                if width: e.set_width_chars(width)
                e.set_text(str(default))
                e.show()
                #hbox.pack_start(l, False)
                #hbox.pack_start(e, False)
                hbox = add_label(e, label)
                c.pack_start(hbox, False)
                return e
            def make_button(label, c=controls):
                b = gtk.Button(label)
                b.show()
                c.pack_start(b, False)
                return b

            scale_adj = gtk.Adjustment(lower=1, upper=10)
            scale_adj.set_value(self.scale)
            scale_adj.connect('value-changed', self.event_scale)
            hscale = gtk.HScale(scale_adj)
            hscale.show()
            controls.pack_start(hscale, False)

            self.ctrl_vizmode = make_combo(self.vizmodes)
            self.ctrl_vizmode.connect("changed", self.event_vizmode)
            self.vizmode = 'mirror'

            update = make_button('Update images')
            update.connect('clicked', self.event_imageclick, None)

            imgseg = gtk.Frame('ImgSeg parameters')
            imgseg.show()
            controls.pack_start(imgseg, False)
            box = gtk.VBox()
            box.show()
            imgseg.add(box)
            imgseg = box

            potts = gtk.Frame('Potts parameters')
            potts.show()
            controls.pack_start(potts, False)
            box = gtk.VBox()
            box.show()
            potts.add(box)
            potts = box

            self.ctrl_mode = make_combo(self.modes, c=imgseg, label="Mode:")
            #self.ctrl_mode.connect("changed", self.event_vizmode)
            self.mode = I.mode

            self.ctrl_distmode = make_combo(self.distmodes, c=imgseg, label="Dist:")
            #self.ctrl_distmode.connect("changed", self.event_vizmode)
            self.distmode = I.distmode

            self.ctrl_VTmode = make_combo(self.VTmodes, c=imgseg, label="VT:")
            self.VTmode = I.VTmode

            self.ctrl_shift = make_combo(self.shifts, c=imgseg, label="Shift:")
            self.shift = I.shift



            self.ctrl_width   = make_entry("width:", self.I.width, c=imgseg, width=6)
            self.ctrl_corrlen = make_entry("corrlen:", self.I.corrlength, c=imgseg, width=6)
            self.ctrl_L               = make_entry("L:", self.I.L, c=imgseg, width=6)
            self.ctrl_beta   = beta   = make_entry(u"β:", 1, c=imgseg, width=6)
            self.ctrl_blur   = blur   = make_entry(u"blur:", self.I.blur, c=imgseg, width=6)

            self.ctrl_gamma  = gamma  = make_entry(u"γ:", .9, c=potts, width=10)
            self.ctrl_trials = trials = make_entry(u"trials:", 3, c=potts, width=6)



            cd_button = make_button("Minimize", c=potts)
            cd_button.connect("clicked", self.do_cd)
            save_button = make_button("Save")
            save_button.connect("clicked", self.do_saveimg)



            #a = gtk.Alignment()
            #a.add(controls)
            #a.show()
            #controls = a

            # Control widgets
            hbox = gtk.HBox(False, False)
            hbox.pack_start(imgbox)
            hbox.pack_start(controls, False)
            hbox.show()

            self.window.add(hbox)
            self.update_img()
            self.window.show()

        def main(self):
            gtk.main()
    W = ImgSegGTK(I, img)
    W.main()

if __name__ == "__main__":
    import sys
    import optparse
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("--crop", type=str, help="Python array slice to crop with, example: 10:50,10:50 .")
    parser.add_option("--resize", type=str, help="Resize array to this fraction of previous size.  Applied after --crop.")
    parser.add_option("-m", "--mode", default='intensity')
    parser.add_option("--dist", help="Override distance function mode")
    parser.add_option("--wait", action='store_true', help="Wait before main loop.")
    parser.add_option("--gtk", action='store_true', help="GTK interface.")
    parser.add_option("--out-suffix", default="", help="Output suffix four automatic output files.")
    parser.add_option("--overlay", default=None, help="Overlay this image on top of the output image.  Should have an alpha channel.")
    parser.add_option("--dry-run", action='store_true', help="Abort after initial set-up (useful for testing parameters).")
    parser.add_option("--threads", default=6, type=int)

    parser.add_option("-g", "--gammas", help="gL,gH[,D] or a dict.")
    parser.add_option("-T", "--trials", type=int, default=3)
    parser.add_option("-R", "--replicas", type=int, default=3)
    parser.add_option("-w", "--width", type=int, default=0, help="Block width (L=w*2+1) (%default)")
    parser.add_option("-c", "--corrlength", type=int, default=1, help="max distance between pixels (%default)")
    parser.add_option("-L", "--length", type=float, default=1.0, help="Length for exp decay or distance cutoff (%default).")
    parser.add_option("--VT", default=None, action='store_true', help="Use variable topology potts model")
    parser.add_option("-s", "--shift")
    parser.add_option("--exp_beta", type=float, help="beta<1 biases towards closer.")
    parser.add_option("--blur", type=float, default=0, help="guassian blur block matrix.")
    parser.add_option("--overlap", action="store_true", help="Allow community overlaps.")
    options, args = parser.parse_args()



    fname = args[0]
    if len(args) >= 2:
        basename = args[1]
    else:
        basename = args[0]
    basename = basename + options.out_suffix

    fullimg = Image(imread(fname))

    overlay=None
    if options.overlay:
        overlay = Image.read(options.overlay)
        color = numpy.divide((232, 118, 0), 255.) # safety orange
        overlay.R[:], overlay.G[:], overlay.B[:] = color

    print "Image shape:", fullimg.shape
    if options.crop:
        print "Cropping",
        fullimg.crop(options.crop)
        if overlay: overlay.crop(options.crop)
        print "new shape", fullimg.shape
    if options.resize:
        print "Resizing",
        fullimg.resize(pcd.util.leval(options.resize))
        if overlay: overlay.resize(pcd.util.leval(options.resize))
        print "new shape", fullimg.shape

    img = fullimg.get()

    print "Making IS object"
    I = ImgSeg(img, mode=options.mode,
               width=options.width, corrlength=options.corrlength,
               #width=3, corrlength=3,
               L=options.length,
               defaultweight=0,
               #cutoff=.06,
               dist=options.dist,
               VTmode='preDefinedOnly' if options.VT else None,
               exp_beta=options.exp_beta,
               blur=options.blur,
               shift=options.shift,
               overlap=options.overlap,
               overlay=overlay,
               )
    #I.dist_func = I.dist_expdecay
    if options.exp_beta:
        I.exp_beta = 1

    if options.gtk:
        gtk_main(I, fullimg)
        sys.exit()

    #I.blocks()

    G = I.G()
    #from fitz import interactnow

    I.do_stats(basename=basename)

    #if options.VT:
    #    print "Enabling variable topology..."
    #    I.enableVT()
    #    #I.enableVT(mode="onlyDefined")
    #    #I.enableVT(mode="preDefined", repelValue=0)
    #    I.do_stats(basename=basename)

    #if options.shift:
    #    I.do_shift(options.shift)

    if options.gammas:
        gammas = pcd.util.eval_gamma_str(options.gammas)
    else:
        gammas = dict(low=.0000001, high=1, density=2)

    # Coordinates: (y,x) (note order).  (0,0) is upper left corner.  y
    # goes down, x goes across.
    #I.do_saveblocks(basename)
    #I.do_connmapfull(basename, (90,10))
    #I.do_connmapfull(basename, (75,38))
    #I.do_connmapfull(basename, (25,25))



    #raw_input('>')
    if options.wait:
        from code import interact ; interact(local=locals(), banner="")
    if options.dry_run:
        print "Aborting, --dry-run enabled"
        sys.exit(0)

    if isinstance(gammas, (float, int)):
        I.do_gamma(basename+'_gamma%09.5f.png'%gammas, gammas)
    else:
        I.do(basename=basename, gammas=gammas, trials=options.trials,
             replicas=options.replicas, threads=options.threads)

    #I.do(basename=out,
    #     gammas=dict(low=.01, high=.01, density=2))

