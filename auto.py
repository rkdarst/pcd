# Richard Darst, May 2012

import numpy

from pcd import Graph, MultiResolution
import pcd.util

class Auto(object):
    trials = 6
    replicas = 5
    gammas = dict(low=.001, high=10, density=4)
    threads = 2
    VT = None
    sparse = True
    overlap = None
    G0callback = None
    MRkwargs = { }
    initial = 'random'  # initial state for minimizations.  If you want
                        # to refine some other communities, add them here.
                        # Format is pcd.Graph.getcmtystate()
    quiet = False

    plot = True
    plotIntermediate = False
    basename = None
    leaf = False
    #leaf_every = False
    modularity = None

    def __init__(self, g=None, options={}):
        self.setOptions(options)
        if self.basename is None:
            raise ValueError("required option 'basename' is missing.")
        if g:
            self.run_g(g)

    def setOptions(self, options):
        for attrname, value in options.items():
            if not hasattr(self, attrname):
                raise ValueError("%s doesn't have attribute '%s'"%(
                    self.__class__.__name__, attrname))
            setattr(self, attrname, value)
    def run_g(self, g, extrafields=(), options={}):
        self.setOptions(options)
        opts = { }

        for a,b,d in g.edges(data=True):
            if getattr(self, 'weight_attr', None):
                d['weight'] = d[self.weight_attr]
            if getattr(self, 'weight', None):
                #g.edge[a][b]['weight'] = options.weight
                d['weight'] = self.weight
            if getattr(self, 'weight_negate', None):
                #g.edge[a][b]['weight'] = -d['weight']
                d['weight'] = -d['weight']


        defaultweight = 1
        if self.VT is not None:
            defaultweight = 0
        G = Graph.fromNetworkX(g, defaultweight=defaultweight)
        if self.quiet:
            G.verbosity = -1
        self.G = G.copy()
        del G._graph
        if self.VT is not None:
            G.make_sparse(default=0.0, cutoff_op=numpy.not_equal)
            G.enableVT(mode='standard')
        elif self.modularity is not None:
            G.enableModularity(mode=self.modularity)
        elif self.sparse:
            G.make_sparse(default='auto')
        if self.G0callback:
            G0 = G.copy()
            G0.cmtyListClear()
            G0 = self.G0callback(G0, g)
        else:
            G0 = None
        #else:
        #    G0=None
        self.run_G(G, G0=G0, extrafields=extrafields)

    def MR(self):
        if hasattr(self, '_MR'):
            return self._MR

        basename = self.basename

        _overlap = self.overlap
        if isinstance(_overlap, (int, float)):
            overlap = _overlap
        elif _overlap:
            overlap = self.trials
        else:
            overlap = False

        self._MR = MultiResolution(overlap=overlap,
                                  minimizerargs=dict(trials=self.trials,
                                                     maxrounds=25,
                                                     minimizer='greedy2',
                                                     initial=self.initial),
                                  output=basename+'-mrvalues.txt',
                #savefigargs=dict(fname=basename+'-MR-gamma%(gamma)10.5f.png')
                                   **self.MRkwargs
                                  )
        return self._MR


    def run_G(self, G, G0=None, options={}, extrafields=()):
        self.G = G
        basename = self.basename
        self.setOptions(options)

        #if self.VT is not None:
        #    G.make_sparse(default=0.0, cutoff_op=numpy.not_equal)
        #    G.enableVT(mode='standard')
        #from fitz import interactnow

        MR = self.MR()
        MR.enable('F1')

        if G0:
            MR.G0 = G0
            MR.enable('mi_0')

        import threading
        lock = threading.Lock()

        def callback(gamma, data, state, MR, MRR):
            G = data['Gmin']
            #if hasattr(self, 'coords'):
            #    G.coords = self.coords
            G.coords = None
            if self.plot or self.plotIntermediate:
                lock.acquire()
                G.viz(show=False,
                      fname=basename+'-MR-gamma%010.5f.png'%gamma,
                      )
                lock.release()
        MR.no_N = True
        MR.run(Gs=[G]*self.replicas,
               gammas=self.gammas,
               threads=self.threads,
               callback=callback,
               extradata={'custom':extrafields},
               )
        self.write()
    def best_communities(self, mode=None, **kwargs):
        """Return the best results from this detection.

        Uses MR.extrema() (extrema of the similarity measures) to make
        a guess at the best communities.  Return list of
        pcd.cmty.Communities object of the results.  These communities
        objects have these extra attributes:

        cmtys.varname: which variable was optimized
        cmtys.maxval: extremum of the varname
        cmtys.name: human-readable description of what was optimized
        cmtys.gamma: the corresponding gamma
        """
        MR = self._MR
        cmtys = MR.best_communities_mapping([self.G]*self.replicas, leaf=self.leaf, **kwargs)
        if mode is not None:
            return cmtys[mode]
        return cmtys

    def write(self):
        if not self.plot:
            return
        basename = self.basename
        MR = self.MR()
        MR.write(basename+'-MR.txt')
        MR.plot(basename+'-MR-plot-g-VI-I-q.png',
                ax1items=['VI', 'In', 'F1', 's_F1'], y1lim=(0,1))
        MR.plot(basename+'-MR-plot-g-VI-I-n.png', ax2items=['n_mean'],
                y1lim=(0,1))
        MR.plot(basename+'-MR-plot-n-VI-I.png', xaxis='n_mean',
                ax1items=['VI',], ax2items=['In', 'F1', 's_F1'],
                xlim=(0,50), y2lim=(0,1))
        MR.plot(basename+'-MR-plot-q-VI-I.png', xaxis='q',
                ax1items=['VI'], ax2items=['In', 'F1', 's_F1'], y2lim=(0,1))
        MR.save(basename+'-state.pickle.gz')



        #ax1items=("VI", "In", "ov_N", "s_F1_ov"),
        #ax2items=("q","n_mean","ov_n_mean"),
        #plotstyles=dict(
        #    ov_N=dict(scale=6),
        #    In=dict(scale=6),
        #    n_mean=dict(scale=1))



if __name__ == "__main__":
    import sys
    import optparse
    from optparse import OptionParser

    parser = OptionParser()
    #parser.add_option("-m", "--mode", default='intensity')
    #parser.add_option("--wait", action='store_true', help="Wait before main loop.")
    parser.add_option("--out-suffix", default="", help="Output suffix four automatic output files.")
    parser.add_option("--dry-run", action='store_true', help="Abort after initial set-up (useful for testing parameters).")
    parser.add_option("--threads", default=6, type=int)
    parser.add_option("--weight", default=-1., type=float, help="Set all graph edges to this weight.")
    parser.add_option("--weight_negate", default=-1., type=float, help="Set all graph edges to this weight.")
    parser.add_option("--weight_attr", help="Name of graph attribute holding weights.")

    parser.add_option("-g", "--gammas", default=".001,100,4", help="gL,gH[,D] or a dict.")
    parser.add_option("-T", "--trials", type=int, default=6)
    parser.add_option("-R", "--replicas", type=int, default=6)
    parser.add_option("--VT", default=None, action='store_true', help="Use variable topology potts model")
    parser.add_option("-s", "--shift")
    parser.add_option("--overlap", action="store_true", help="Allow community overlaps.")
    options, args = parser.parse_args()

    fname = args[0]
    if len(args) >= 2:
        basename = args[1]
    else:
        basename = args[0]
    basename = basename + options.out_suffix


    import networkx
    if fname.endswith(".gml"):
        g = networkx.read_gml(fname)

    coords = networkx.spring_layout(g)

    #print coords
    #for n, c in coords.iteritems():
    #    print n, c, g.nodes()
    #    g.node[n]['pos'] = c

    #if options.wait:
    #    from code import interact ; interact(local=locals(), banner="")

    alloptions = [ o.dest for o in parser._long_opt.itervalues()
                   if o.dest is not None ]
    #for n, o in parser._long_opt.iteritems():
    #    print n, o.dest, o
    #print alloptions

    A = Auto(pcd.util.AttrToDict(options, values=alloptions))
    A.setOptions(dict(basename=basename))

    A.gammas = pcd.util.eval_gamma_str(options.gammas)

    A.run_g(g, dict(coords=coords,
                    overlap=options.overlap))
