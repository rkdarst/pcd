# Richard Darst, September 2011

from math import log, exp, floor, ceil
import numpy
import time

import util
from util import LogInterval

def recursive_dict_update(d, dnew):
    for k,v in dnew.iteritems():
        if k in d and isinstance(d[k],dict) and isinstance(v,dict):
            recursive_dict_update(d[k], v)
        else:
            d[k] = v
    return d



class MultiResolutionCorrelation(object):
    def __init__(self, gamma, Gs, trials=None,
                 nhistbins=50, **kwargs):
        self.calc(gamma, Gs, nhistbins=nhistbins)
        self.replicas = len(Gs)
        self.trials = trials
        for k,v in kwargs.iteritems():
            setattr(self, k, v)

    def calc(self, gamma, Gs, nhistbins=50):
        pairs = [ ]
        pairsOverlap = [ ]

        overlapGs = [ G.copy() for G in Gs ]
        [ G.trials(gamma, trials=10, initial='current',
                   minimizer='overlapMinimize')
          for G in overlapGs ]

        Gmin = min(Gs, key=lambda G: G.energy(gamma))
        self.cmtyState = tuple(G.getcmtystate() for G in Gs)
        self.cmtyStateBest = Gmin.getcmtystate()

        GminOverlap = Gmin.copy()
        #GminOverlap.setOverlap(True)
        GminOverlap.trials(gamma, trials=10, initial='current',
                           minimizer='overlapMinimize')
        # Do analysis using GminOverlap


        for i, G0 in enumerate(Gs):
            for j, G1 in enumerate(Gs[i+1:]):
                pairs.append((G0, G1))
                pairsOverlap.append((overlapGs[i], overlapGs[j]))

        Is = [ util.mutual_information(G0,G1) for G0,G1 in pairs ]

        VI = numpy.mean([G0.entropy + G1.entropy - 2*mi
                         for ((G0,G1), mi) in zip(pairs, Is)])
        In = numpy.mean([2*mi / (G0.entropy + G1.entropy)
                         for ((G0,G1), mi) in zip(pairs, Is)
                         if G0.q!=1 or G0.q!=1])

        NmiO = numpy.mean([ util.mutual_information_overlap(G1, G2)
                            for G1,G2 in pairsOverlap])


        self.fieldnames = ("gamma", "q", "q_std", "qmin",
                           "E", "entropy",
                           "I", "VI", "In", "NmiO",
                           "n_mean")
        self.gamma   = gamma
        self.q       = numpy.mean(tuple(G.q for G in Gs))
        self.q_std   = numpy.std(tuple(G.q for G in Gs), ddof=1)
        self.qmin    = Gmin.q

        self.E       = numpy.mean(tuple(G.energy(gamma) for G in Gs))
        self.entropy = numpy.mean(tuple(G.entropy for G in Gs))

        self.I       = numpy.mean(Is)
        self.VI      = VI
        self.In      = In
        self.NmiO    = NmiO

        self.N       = N = Gs[0].N
        self.n_mean  = sum(tuple(item[0]*item[1]/float(G.q) for G in Gs for
                                 item in G.n_counts().items())) / ( float(len(Gs)) )

        n_hists=[]
        edges=numpy.logspace(0,numpy.log10(N+1),num=nhistbins,endpoint=True)
        for G in Gs:
            n_counts = []
            histdata = tuple(item for item in G.n_counts().items())
            for number,count in histdata:
                n_counts = n_counts+[number]*count*number
            n_counts = sorted(n_counts)
            n_hist,self.n_hist_edges = numpy.histogram(n_counts, bins=edges)
            n_hists.append(n_hist/float(sum(n_hist)))
        n_hists_array = numpy.array(n_hists)
        self.n_hist = n_hists_array.mean(axis=0)

    def getGs(self, Gs):
        """Return the list of all G replicas.

        This recreates the list of all Gs, from stored states.  You
        need to provide a graph instance G, which has the same node
        arrangement as the original graph(s) used in the multi-replica
        analysis.
        """
        Gs = [ G.copy() for G in Gs ]
        for G, state in zip(Gs, self.cmtyState):
            G.setcmtystate(state)
        return Gs

class MultiResolution(object):
    """Class to do full multi-resolution analysis.
    """
    _callback = None
    _output = None
    _lock = None
    def __init__(self, output=None, savefigargs=None,
                 minimizer='trials',
                 minkwargs=dict(minimizer='minimize', trials=10),
                 plotargs=None,
                 calckwargs={}):
        """Create a multiresolution helper system.

        output: save table of gamma, q, H, etc, to this file.

        minimizer: name of minimizer to use.  This should be a method
        name on the Graph object.

        minkwargs: keyword arguments to be passed to the minimizer.

        savefigargs: keyword arguments (dictionary) to pass to
        .savefig() of the minimum G of each gamma.  To save to a
        particular filename, use
        savefigargs={'fname':'filename%(gamma)f.png'}


        Simple greedy minimization:
        minimiser=default, minkwargs=dict(trials=NTRIALS)

        Annealing minimization:
        minimizer=default, minkwargs=dict(minimizer='anneal', trials=NTRIALS)

        Annealing minimization, no trials:
        minimizer='anneal', minkwargs={args to anneal, ...}


        """
        self._output = output
        self._savefigargs = savefigargs
        self._plotargs = plotargs
        self._minimizer = minimizer
        self._minimizerkwargs = minkwargs
        self._calckwargs = calckwargs

        self._data = { } #collections.defaultdict(dict)
    def do_gamma(self, gamma, callback=None):
        """Worker thread for gammas.

        For one gamma, look up the respective gamma and do that analysis.
        """
        if callback is None:
            callback = self._callback
        minGs = [ ]
        for G in self._Gs:
            if self._lock:  self._lock.acquire()
            G = G.copy()
            if self._lock:  self._lock.release()
            getattr(G, self._minimizer)(gamma, **self._minimizerkwargs)
            minGs.append(G)

        # Do the information theory VI, MI, In, etc, stuff on our
        # minimized replicas.
        self._data[gamma] = MultiResolutionCorrelation(
            gamma, minGs, trials=self.trials)
        # Save output to a file.
        if self._output is not None and \
               (not hasattr(self, '_writelock' or self._writelock.acquire(False))):
            self.write(self._output)
            if hasattr(self, '_writelock'):
                self._writelock.release()
        if self._savefigargs is not None:
            kwargs = self._savefigargs.copy()
            kwargs['fname'] = kwargs['fname']%{'gamma':gamma}
            G.remapCommunities(check=False)
            G.savefig(**kwargs)
        if self._plotargs is not None:
            #kwargs = self._savefigargs.copy()
            #kwargs['fname'] = kwargs['fname']%{'gamma':gamma}
            self.plot(**self._plotargs)

        # Run callback if we have it.
        if callback:
            totalminimum = min(minGs, key=lambda x: G.energy(gamma))
            replica_number = minGs.index(totalminimum)
            callback(G=totalminimum, gamma=gamma,
                     mrc=self._data[gamma],
                     mr=self,
                     replica_number=replica_number)
    def _thread(self):
        """Thread worker - do gammas until all are exhausted.

        One of these are spawned for each thread, they run until all
        gammas have been finished."""
        # The queue.join() may return by the time we get to the
        # except: clause, where self._EmptyException has already been
        # deleted.  Thus, we need to make a local copy before the
        # try/except block.
        EmptyException = self._EmptyException
        try:
            while True:
                gamma = self._queue.get_nowait()
                self.do_gamma(gamma)
                self._queue.task_done()
        except EmptyException:
            return
    def do(self, Gs, gammas=None, logGammaArgs=None,
           trials=None, threads=1, callback=None):
        """Do multi-resolution analysis on replicas Gs with `trials` each."""
        # If multiple threads requested, do that version instead:
        if gammas is None and logGammaArgs is not None:
            gammas = LogInterval(**logGammaArgs).values()
        if gammas is None:
            raise ValueError("Either gammas or logGammaArgs must be given.")
        if trials is not None:
            self._minimizerkwargs['trials'] = trials
        if threads > 1:
            return self.do_mt(Gs, gammas, trials=trials, threads=threads,
                              callback=callback)

        # Non threaded version:
        self._callback = callback
        self._Gs = Gs
        self.replicas = len(Gs)
        self.trials = trials
        for gamma in gammas:
            self.do_gamma(gamma)
            if self._output:
                self.write(self._output)
        del self._Gs
        del self._callback # unmasks the class object's _callback=None

    def do_mt(self, Gs, gammas, trials=10, threads=2, callback=None):
        """Do multi-resolution analysis on replicas Gs with `trials` each.

        Multi-threaded version."""
        self._Gs = Gs
        self._callback = callback
        self.replicas = len(Gs)
        self.trials = trials

        import threading
        import Queue
        self._EmptyException = Queue.Empty
        self._lock = threading.Lock()
        self._writelock = threading.Lock()
        self._queue = Queue.Queue()

        for gamma in gammas:
            self._queue.put(gamma)

        threadlst = [ ]
        for tid in range(threads):
            t = threading.Thread(target=self._thread)
            #t.daemon = True
            t.start()
            threadlst.append(t)

        self._queue.join()
        del (self._Gs, self._EmptyException, self._lock, self._writelock,
             self._queue)
        del self._callback # unmasks class object's _callback=None

    _calcDone = False
    def calc(self):
        """Put information from all of the MRCs on self.

        This takes
        """
        # Cache to only do it once.  Need a way to override this.
        #if self._calcDone:
        #    return
        #self._calcDone = True


        MRCs = [ MRC for i, MRC in sorted(self._data.iteritems()) ]

        self.fieldnames = MRCs[0].fieldnames
        for fieldname in self.fieldnames:
            setattr(self, fieldname,
                    numpy.asarray([getattr(mrc, fieldname) for mrc in MRCs]))

        self.N         = N         = MRCs[0].N

        self.n_hist       = n_hist       = [mrc.n_hist       for mrc in MRCs]
        self.n_hist_edges = n_hist_edges = [mrc.n_hist_edges for mrc in MRCs]

    def write(self, fname):
        """Save multi-resolution data to a file."""
        self.calc()
        f = open(fname, 'w')
        print >> f, "#", time.ctime()
        print >> f, "# replicas:", self.replicas
        print >> f, "# trials per replica:", self.trials
        print >> f, "#" + " ".join(self.fieldnames)
        for i in range(len(self.gamma)):
            print >> f, self.gamma[i],
            for name in self.fieldnames:
                if name == "gamma": continue
                print >> f, getattr(self, name)[i],
            print >> f

    def plot_nhists(self, fname):
        from matplotlib.backends.backend_agg import Figure, FigureCanvasAgg
        from matplotlib.backends.backend_pdf import PdfPages
        f = Figure()
        c = FigureCanvasAgg(f)

        ax  = f.add_subplot(111)
        pp = PdfPages(fname)

        for i in range(len(self.n_hist)):
            histdata=self.n_hist[i]
            histedges=self.n_hist_edges[i]
            ax.cla()
            ax.set_title("$\gamma=%f$"%self.gamma[i])
            ax.bar( histedges[:-1], histdata, width=(histedges[1:]-histedges[:-1]) )
            ax.set_xscale('log')
            ax.set_xlim( 0.1 , self.N )
            ax.set_ylim( 0   , 1      )
            pp.savefig(f)

        pp.close()

    def plot_nmean(self, fname=None):
        if fname:
            from matplotlib.backends.backend_agg import Figure, FigureCanvasAgg
            f = Figure()
            c = FigureCanvasAgg(f)
        else:
            from matplotlib import pyplot
            f = pyplot.gcf()

        ax_vi  = f.add_subplot(111)
        ax_vi.set_xscale('log')
        ax_n = ax_vi.twinx()

        l2 = ax_vi.plot(self.gamma, self.VI,
                        color='blue')
        l3 = ax_vi.plot(self.gamma, self.In, '--',
                        color='red')
        l = ax_n.plot(self.gamma, self.n_mean,':',c='black')

        ax_n.legend((l2, l3,l), ("$VI$", "$I_n$","$\langle n \\rangle$"), loc=0)
        ax_n.set_ylabel('$\langle n \\rangle$')
        ax_n.set_yscale('log')
#        ax_n.set_ylim(bottom=0)

        ax_vi.set_ylabel('$VI$ and $I_n$')
        ax_vi.set_xlabel('$\gamma$')

        if fname:
            c.print_figure(fname, bbox_inches='tight')

    def plot(self, fname=None, ax1items=['VI', 'In', 'NmiO'], ax2items=['q'],
             plotstyles={}):
        """Plot things from a multiresolution plot.

        ax1items: List of items to plot on axis 1 (left)
        ax1items: List of items to plot on axis 2 (right)

        Each item should be a name of something on the
        MultiResolutionCorrelation object.
        """
        # Recalculate all the stuff
        self.calc()
        # Prepare plot figure
        if fname:
            from matplotlib.backends.backend_agg import Figure, FigureCanvasAgg
            f = Figure()
            c = FigureCanvasAgg(f)
        else:
            from matplotlib import pyplot
            f = pyplot.gcf()

        # Add axes
        ax1  = f.add_subplot(111)
        ax1.set_xscale('log')
        ax2 = ax1.twinx()

        # Defaults (move outside of the function eventually)
        legenditems = [ ]
        defaultplotstyles = {
            'VI':     dict(label='$VI$',color='blue', linestyle='-'),
            'In':     dict(label='$I_n$',color='red', linestyle='--'),
            'NmiO':   dict(label='$N$',color='green', linestyle='-.'),
            'q':      dict(label='$q$',color='black', linestyle=':'),
            '<n>':    dict(label='$<n>$',color='black', linestyle=':'),
            'entropy':dict(label='$H$',color='blue', linestyle='--'),
            None:   dict(),  # defaults
            }
        plotstyles = recursive_dict_update(defaultplotstyles, plotstyles)
        # Fill axes
        for item in ax1items:
            plotstyle = plotstyles.get(item, {}).copy()
            scale = plotstyle.pop('scale', 1)
            l = ax1.plot(self.gamma, getattr(self, item)*scale,
                         **plotstyle)[0]
            legenditems.append((l, plotstyle.get('label', item)))
        for item in ax2items:
            plotstyle = plotstyles.get(item, {}).copy()
            scale = plotstyle.pop('scale', 1)
            #plotstyle2 = plotstyle.copy()
            #plotstyle2.pop('scale', None)
            l = ax2.plot(self.gamma, getattr(self, item)*scale,
                         **plotstyle)[0]
            legenditems.append((l, plotstyle.get('label', item)))

        # Legends and labels
        ax2.legend(*zip(*legenditems), loc=0)
        ax2.set_ylim(bottom=0)

        def axLabels(axItems):
            """Look up axis labels for a list of items."""
            for item in axItems:
                plotstyle = plotstyles.get(item, {})
                label = plotstyle.get('label', item)
                if 'scale' in plotstyle:
                    label = (r'$%s\times$'%plotstyle['scale'])+label
                yield label
        ax1labels = ', '
        ax1.set_ylabel(', '.join(axLabels(ax1items)))
        ax2.set_ylabel(', '.join(axLabels(ax2items)))
        ax1.set_xlabel('$\gamma$')

        if fname:
            c.print_figure(fname, bbox_inches='tight')

    def viz(self):
        import matplotlib.pyplot as pyplot
        gamma     = self.gamma
        q         = self.q
        E         = self.E
        entropy   = self.entropy

        I         = self.I
        VI        = self.VI
        In        = self.In
        pyplot.xlim(gamma[0], gamma[-1])
        pyplot.semilogx(gamma, q)
        pyplot.xlabel('\gamma')
        pyplot.ylabel('q')
        #pyplot.show()
        pyplot.ion()
        from fitz import interactnow
