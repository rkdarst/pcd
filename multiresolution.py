# Richard Darst, September 2011

from math import log, exp, floor, ceil
import numpy
import time
import types

import logging
logger = logging.getLogger('pcd.mr')

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
    calcMethods = [ ]
    def __init__(self, gamma, Gs,
                 nhistbins=50, calcSettings={}):
        self.gamma = gamma
        self.replicas = len(Gs)
        self.calcSettings = calcSettings

        Gmin = min(Gs, key=lambda G: G.energy(gamma))
        GminIndex = Gs.index(Gmin)
        self.GminIndex = GminIndex
        self.cmtyState = tuple(G.getcmtystate() for G in Gs)
        self.cmtyStateBest = Gmin.getcmtystate()

        if calcSettings.get('overlap'):
            overlapTrials = calcSettings.get('overlap')
            #if overlapTrials is True:
            #    overlapTrials = self.trials
            #self.overlapTrials = overlapTrials

            overlapGs = [ G.copy() for G in Gs ]
            [ G.trials(gamma, trials=overlapTrials, initial='current',
                       minimizer=calcSettings.get('overlapMinimizer',
                                                'overlapMinimize'))
              for G in overlapGs ]
            self.cmtyStateOverlap = tuple(G.getcmtystate() for G in overlapGs)
        else:
            overlapGs = None

        self.calc(gamma, Gs=Gs, nhistbins=nhistbins,
                  Gmin=Gmin, overlapGs=overlapGs,
                  calcSettings=calcSettings)


    def calc(self, gamma, Gs=None, nhistbins=50,
             Gmin=None, overlapGs=None, graph_list=None,
             calcSettings=None):
        """Calculate properties of minimized replicas.

        You need to give one of these sets of arguments to this method:

        Gs, Gmin, overlapGs: these are the *minimized* configurations.

        graph_list: if this is given, then recalculate everything with
        this initial list.  Reload state from self.cmtyState*.
        graph_list is *not* input minimized configurations.

        """
        # In all of these functions: graph_list is the original
        # replica Graph objects, but *not* necessarily the minimized
        # versions.  This variable needs to have its state reloaded.
        # The Gs variables are the minimized configurations.
        if Gs is None:
            Gs = self.getGs(graph_list)
        if Gmin is None:
            Gmin = self.getGmin(graph_list)
        if calcSettings.get('overlap') and overlapGs is None:
            overlapGs = self.getOverlapGs(graph_list)

        self.fieldnames = ('gamma', )
        self.gamma = gamma

        data = { }
        data['Gs'] = Gs
        data['Gmin'] = Gmin
        data['overlapGs'] = overlapGs
        data['gamma'] = gamma

        pairIndexes = [ ]
        pairStyle = calcSettings.get('pairStyle', 'all')
        if pairStyle == "all":
            for i in range(len(Gs)):
                for j in range(i+1, len(Gs)):
                    pairIndexes.append((i, j))
        elif pairStyle == "adjacent":
            pairIndexes.extend(zip(range(len(Gs)-1), range(1, len(Gs))))
        else:
            raise ValueError("pairStyle %s not known"%pairStyle)
        data['pairIndexes'] = pairIndexes


        # Not moved to a method yet, this is special.
        N = Gs[0].N
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


        # Calculate things defined in other methods:
        for func in self.calcMethods:
            if isinstance(func, types.MethodType):
                ret = func(data=data, calcSettings=calcSettings)
            ret = func(self, data=data, calcSettings=calcSettings)
            for name, value in ret:
                self.fieldnames += (name, )
                setattr(self, name, value)

    def calc_static(self, data, calcSettings):
        Gs = data['Gs']
        return [('n_mean', sum(G.n_mean() for G in Gs)/float(len(Gs))),
                ('q',      numpy.mean(tuple(G.q for G in Gs))),
                ('q_std',  numpy.std(tuple(G.q for G in Gs), ddof=1)),
                ('qmin',   data['Gmin'].q),
                ('nodes',  Gs[0].N),
                ]
    calcMethods.append(calc_static)
    def calc_thermo(self, data, calcSettings):
        Gs = data['Gs']
        E = numpy.mean(tuple(G.energy(data['gamma']) for G in Gs))
        if Gs[0].oneToOne:
            entropy = numpy.mean(tuple(G.entropy for G in Gs))
        else:
            entropy = 0
        return [('E', E),
                   ('entropy', entropy)
                ]

    calcMethods.append(calc_thermo)
    def calc_changes(self, data, calcSettings):
        Gs = data['Gs']
        returns = [ ]
        # This records the number of rounds needed to perform the
        # minimization.
        if hasattr(Gs[0], 'nChanges'):
            nChanges = numpy.mean([G.nChanges for G in Gs], axis=0)
            for i, value in enumerate(nChanges):
                name = "nChanges%d"%i
                returns.append((name, value))
        return returns
    calcMethods.append(calc_changes)
    def calc_multualInformation(self, data, calcSettings):
        Gs = data['Gs']
        pairIndexes = data['pairIndexes']
        Is = [ util.mutual_information(Gs[i],Gs[j])
               for i,j in pairIndexes ]
        VI = numpy.mean([Gs[i].entropy + Gs[j].entropy - 2*mi
                         for ((i,j), mi) in zip(pairIndexes, Is)])
        In = numpy.mean([2*mi / (Gs[i].entropy + Gs[j].entropy)
                         for ((i,j), mi) in zip(pairIndexes, Is)
                         if Gs[i].q!=1 or Gs[j].q!=1])
        returns = [('I',  numpy.mean(Is)),
                   ('VI', VI),
                   ('In', In),
                   ]

        Nmi = numpy.mean(
            [ util.mutual_information_overlap(Gs[i], Gs[j])
              for i,j in data['pairIndexes'] ])
        returns.append(('N',  Nmi))

        return returns
    calcMethods.append(calc_multualInformation)
    def calc_overlap(self, data, calcSettings):
        overlapGs = data.get('overlapGs', None)
        if overlapGs is None:
            return [ ]
        NmiO = numpy.mean(
            [ util.mutual_information_overlap(overlapGs[i], overlapGs[j])
              for i,j in data['pairIndexes'] ])
        n_mean_ov = sum(G.n_mean() for G in overlapGs)/float(len(overlapGs))

        returns = [('ov_N',      NmiO),
                   ('ov_n_mean', n_mean_ov),
                   ('ov_E',    numpy.mean(tuple(G.energy(data['gamma'])
                                                for G in overlapGs))),
                   ]

        # This records the number of rounds needed to perform the
        # minimization.
        data2 = data.copy()
        data2['Gs'] = data['overlapGs']
        ov_ret = self.calc_changes(data=data2, calcSettings=calcSettings)
        returns.extend(('ov_'+name, value)
                       for name, value in ov_ret)

        return returns
    calcMethods.append(calc_overlap)

    def getGs(self, graph_list):
        """Return the list of all G replicas.

        This recreates the list of all Gs, from stored states.  You
        need to provide a graph list Gs, which has the same node
        arrangement as the original graph list used in the
        multi-replica analysis.
        """
        Gs = [ G.copy() for G in graph_list ] # copy is assumed in callers
        for G, state in zip(Gs, self.cmtyState):
            G.setcmtystate(state)
        return Gs
    def getOverlapGs(self, graph_list):
        overlapGs = self.getGs(graph_list) # this method ensures we have a copy
        if hasattr(self, 'cmtyStateOverlap'):
            for G, state in zip(overlapGs, self.cmtyStateOverlap):
                G.setcmtystate(state)
        else:
            [ G.trials(self.gamma, trials=self.overlapTrials, initial='current',
                       minimizer='overlapMinimize')
              for G in overlapGs ]
            self.cmtyStateOverlap = tuple(G.getcmtystate() for G in overlapGs)
        return overlapGs

    def getGmin(self, graph_list):
        if hasattr(self, 'GminIndex'):
            Gmin = graph_list[self.GminIndex]
            Gmin.setcmtystate(self.cmtyStateBest)
            return Gmin
        else:
            Gs = self.getGs(graph_list)
            Gmin = min(Gs, key=lambda G: G.energy(self.gamma))
            return Gmin
    def getOverlapGmin(self, graph_list):
        Gs = self.getOverlapGs(graph_list)
        #minG = min(Gs, key=lambda G: G.energy(self.gamma))
        minG = Gs[self.GminIndex]
        del Gs
        import gc ; gc.collect()
        return minG

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
                 #overlap=False,
                 #pairstyle='all',
                 calcSettings={}):
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
        self._calcSettings = calcSettings

        self._data = { } #collections.defaultdict(dict)
    def do_gamma(self, gamma, callback=None):
        """Worker thread for gammas.

        For one gamma, look up the respective gamma and do that analysis.
        """
        if callback is None:
            callback = self._callback
        logger.info("Thread %d MR-minimizing g=%f"%(self.thread_id(), gamma))
        minGs = [ ]
        for G in self._Gs:
            if self._lock:  self._lock.acquire()
            G = G.copy()
            if self._lock:  self._lock.release()
            getattr(G, self._minimizer)(gamma, **self._minimizerkwargs)
            minGs.append(G)
        logger.debug("Thread %d MR-minimizing g=%f: done"%(
                                                      self.thread_id(), gamma))
        # Do the information theory VI, MI, In, etc, stuff on our
        # minimized replicas.
        logger.info("Thread %d initializing MRC g=%f"%(
                                              self.thread_id(), gamma))
        self._data[gamma] = MultiResolutionCorrelation(
            gamma, minGs,
            calcSettings=self._calcSettings)
        logger.debug("Thread %d initializing MRC g=%f: done"%(
                                              self.thread_id(), gamma))
        # Save output to a file.
        if self._output is not None and self.lockAcquire(blocking=False):
            logger.debug("Thread %d write"%self.thread_id())
            self.write(self._output)
            self.lockRelease()
        if self._savefigargs is not None and self.lockAcquire(blocking=True):
            logger.debug("Thread %d savefig"%self.thread_id())
            minG = min(minGs, key=lambda x: G.energy(gamma))
            kwargs = self._savefigargs.copy()
            kwargs['fname'] = kwargs['fname']%{'gamma':gamma}
            G.remapCommunities(check=False)
            G.savefig(**kwargs)
            self.lockRelease()
        if self._plotargs is not None and self.lockAcquire(blocking=False):
            logger.debug("Thread %d plot"%self.thread_id())
            self.plot(**self._plotargs)
            self.lockRelease()
        # Run callback if we have it.
        if callback:
            logger.debug("Thread %d callback"%self.thread_id())
            minG = min(minGs, key=lambda x: G.energy(gamma))
            replica_number = minGs.index(minG)
            callback(G=minG, gamma=gamma,
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
                logger.debug("Thread %d begin gamma=%f"%(
                    self.thread_id(), gamma))
                self.do_gamma(gamma)
                logger.debug("Thread %d end gamma=%f"%(
                    self.thread_id(), gamma))
                self._queue.task_done()
        except EmptyException:
            logger.debug("Thread terminated: %d", self.thread_id())
            return
    def thread_id(self):
        import threading
        return threading.current_thread().ident
    def lockAcquire(self, name='writelock', blocking=False):
        """Acquire lock self._NAME.

        blocking: if true, block for lock to acquire.  If false, don't
        block, return fals if lock not acquired."""
        logger.debug('Thread %d attempting lock acquisition: %s',
                     self.thread_id(), name)
        if hasattr(self, '_'+name):
            locked = getattr(self, '_'+name).acquire(blocking)
            logger.debug('Thread %d attempting lock acquisition: %s: %s',
                         self.thread_id(), name, locked)
            return locked
        # If we don't have self._LOCKNAME, we don't lock so always return true.
        return True
    def lockRelease(self, name='writelock'):
        """Release lock self.NAME."""
        logger.debug('Thread %d attempting release: %s',
                     self.thread_id(), name)
        if hasattr(self, '_'+name):
            logger.debug('Thread %d attempting release (phase 2): %s',
                         self.thread_id(), name)
            return getattr(self, '_'+name).release()
    def do(self, Gs, gammas=None, logGammaArgs=None,
           threads=1, callback=None):
        """Do multi-resolution analysis on replicas Gs with `trials` each."""
        # If multiple threads requested, do that version instead:
        if gammas is None and logGammaArgs is not None:
            gammas = LogInterval(**logGammaArgs).values()
        if gammas is None:
            raise ValueError("Either gammas or logGammaArgs must be given.")
        if threads > 1:
            return self.do_mt(Gs, gammas, threads=threads,
                              callback=callback)

        # Non threaded version:
        self._callback = callback
        self._Gs = Gs
        self.replicas = len(Gs)
        for gamma in gammas:
            self.do_gamma(gamma)
            if self._output:
                self.write(self._output)
        del self._Gs
        del self._callback # unmasks the class object's _callback=None

    def do_mt(self, Gs, gammas, threads=2, callback=None):
        """Do multi-resolution analysis on replicas Gs with `trials` each.

        Multi-threaded version."""
        self._Gs = Gs
        self._callback = callback
        self.replicas = len(Gs)

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
    def recalc(self, graph_list):
        """Rerun .calc() method an all MRCs.

        This is a higher-cost method (reruns various calculations, but
        doesn't need to rerun minimizations.)"""
        for gamma, MRC in self._data.iteritems():
            MRC.calc(graph_list)
    def calc(self):
        """Compiles data from all MRCs onto self.

        This takes the data from the MRCs on self, and puts it in
        lists by the field names.  This is the first step to plotting,
        writing, etc.  This is a low-cost method.
        """
        MRCs = [ MRC for i, MRC in sorted(self._data.iteritems()) ]

        self.fieldnames = MRCs[0].fieldnames
        for fieldname in self.fieldnames:
            setattr(self, fieldname,
                    numpy.asarray([getattr(mrc, fieldname) for mrc in MRCs]))

        #self.N         = N         = MRCs[0].N

        self.n_hist       = n_hist       = [mrc.n_hist       for mrc in MRCs]
        self.n_hist_edges = n_hist_edges = [mrc.n_hist_edges for mrc in MRCs]

    def write(self, fname):
        """Save multi-resolution data to a file."""
        self.calc()
        f = open(fname, 'w')
        print >> f, "#", time.ctime()
        print >> f, "# replicas:", self.replicas
        print >> f, "# minimizer and kwargs:", self._minimizer, \
                                               self._minimizerkwargs

        print >> f, "#", " ".join("%d:%s"%(i+1,x)
                                  for i,x in enumerate(self.fieldnames))
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
             xaxis='gamma',
             plotstyles={}):
        """Plot things from a multiresolution plot.

        ax1items: List of items to plot on axis 1 (left)
        ax1items: List of items to plot on axis 2 (right)

        Each item should be a name of something on the
        MultiResolutionCorrelation object.
        """
        # Recalculate all the stuff
        self.calc()
        # remove NmiO if we don't have that attribute
        if ax1items==['VI', 'In', 'NmiO'] and 'NmiO' in ax1items:
            if not hasattr(self, 'NmiO'):
                ax1items.remove('NmiO')
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
            'gamma':  dict(label='$\gamma$',color='blue', linestyle='-'),
            'VI':     dict(label='$VI$',color='blue', linestyle='-'),
            'In':     dict(label='$I_n$',color='red', linestyle='--'),
            'N':      dict(label='$N$',color='green', linestyle='-.'),
            'ov_N':   dict(label='$N_{ov}$',color='green', linestyle='-.'),
            'q':      dict(label='$q$',color='black', linestyle=':'),
            'n_mean': dict(label='$<n>$',color='blue',linestyle=':'),
            'n_mean_ov':dict(label='$<n_{ov}>$',color='green',linestyle=':'),
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
        ax1.set_xlabel(defaultplotstyles.get('gamma',
                                             dict(label=xaxis))['label'])

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
