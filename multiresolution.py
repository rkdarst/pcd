# Richard Darst, September 2011

from math import log, exp, floor, ceil
import numpy
import threading
import time
import types

import logging
logger = logging.getLogger('pcd.mr')

import util
from util import LogInterval
from fitz.mathutil import Averager

def recursive_dict_update(d, dnew):
    for k,v in dnew.iteritems():
        if k in d and isinstance(d[k],dict) and isinstance(v,dict):
            recursive_dict_update(d[k], v)
        else:
            d[k] = v
    return d


class GammaData(object):
    def __init__(self, gamma):
        self.gamma = gamma
        self.data = { }
    def add(self, name, value):
        if not hasattr(self, name):
            self.data[name] = Averager()
        self.data[name].add(value)

class MultiResolution(object):
    calcMethods = [ ]
    nhistbins = 50
    pairStyle = 'all'

    def __init__(self, overlap=5,
                 minimizer='trials',
                 minimizerargs=dict(minimizer='greedy', trials=10),
                 output=None, savefigargs=None,plotargs=None,
                 ):
        self.fieldnames = [ 'gamma', ]
        self._data = { }

        self.overlap = overlap
        self.minimizer = minimizer
        self.minimizerargs = minimizerargs
        self.output = output
        self.savefigargs = savefigargs
        self.plotargs = plotargs
        self._lock = threading.Lock()
    def __getstate__(self):
        state = self.__dict__
        del state['_lock']
        return state
    def __setstate__(self, state):
        self.__dict__ = state
        self._lock = threading.Lock()

    def runner(self):
        return MRRunner(self)
    def run(self, *args, **kwargs):
        runner = self.runner()
        return runner.do(*args, **kwargs)

    def add(self, gamma, data, state, ):
        """Calculate properties of minimized replicas.

        gamma: our gamma

        data: the minimized Gs

        state:
        """
        # Store number of nodes
        if not hasattr(self, 'N'):
            self.N = data['Gs'][0].N
        assert all(self.N == G.N for G in data['Gs'])


        pairIndexes = [ ]
        pairStyle = getattr(self, 'pairStyle', 'all')
        if pairStyle == "all":
            for i in range(len(data['Gs'])):
                for j in range(i+1, len(data['Gs'])):
                    pairIndexes.append((i, j))
        elif pairStyle == "adjacent":
            pairIndexes.extend(zip(range(len(Gs)-1), range(1,len(data['Gs']))))
        else:
            raise ValueError("pairStyle %s not known"%pairStyle)
        data['pairIndexes'] = pairIndexes


        ## Not moved to a method yet, this is special.
        #N = Gs[0].N
        #n_hists=[]
        #edges=numpy.logspace(0,numpy.log10(N+1),num=nhistbins,endpoint=True)
        #for G in Gs:
        #    n_counts = []
        #    histdata = tuple(item for item in G.n_counts().items())
        #    for number,count in histdata:
        #        n_counts = n_counts+[number]*count*number
        #    n_counts = sorted(n_counts)
        #    n_hist,self.n_hist_edges = numpy.histogram(n_counts, bins=edges)
        #    n_hists.append(n_hist/float(sum(n_hist)))
        #n_hists_array = numpy.array(n_hists)
        #self.n_hist = n_hists_array.mean(axis=0)


        # Calculate things defined in other methods:
        returns = [ ]
        for func in self.calcMethods:
            if isinstance(func, types.MethodType):
                ret = func(data=data, settings=self)
            else:
                ret = func(self, data=data, settings=self)
            returns.extend(ret)
        self._addValues(gamma, returns)

    def _addValues(self, gamma, namevals):
        """Add a list of (name,value) pairs to the corresponding gamma data"""
        with self._lock:
            for name, val in namevals:
                if name not in self.fieldnames:
                    self.fieldnames.append(name)
                if gamma not in self._data:
                    self._data[gamma] = GammaData(gamma=gamma)
                self._data[gamma].add(name, val)

    #
    # Calculation methods
    #
    def calc_static(self, data, settings):
        Gs = data['Gs']
        return [('n_mean', sum(G.n_mean() for G in Gs)/float(len(Gs))),
                ('q',      numpy.mean(tuple(G.q for G in Gs))),
                ('q_std',  numpy.std(tuple(G.q for G in Gs), ddof=1)),
                ('qmin',   data['Gmin'].q),
                ('nodes',  Gs[0].N),
                ]
    calcMethods.append(calc_static)
    def calc_thermo(self, data, settings):
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
    def calc_changes(self, data, settings):
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
    def calc_multualInformation(self, data, settings):
        Gs = data['Gs']
        pairIndexes = data['pairIndexes']
        def _entropy(G):
            if G.oneToOne: return G.entropy
            else: return 1
        Is = [ util.mutual_information(Gs[i],Gs[j])
               for i,j in pairIndexes ]
        VI = numpy.mean([_entropy(Gs[i]) + _entropy(Gs[j]) - 2*mi
                         for ((i,j), mi) in zip(pairIndexes, Is)])
        In = numpy.mean([2*mi / (_entropy(Gs[i]) + _entropy(Gs[j]))
                         for ((i,j), mi) in zip(pairIndexes, Is)
                         if Gs[i].q!=1 or Gs[j].q!=1])
        returns = [('I',  numpy.mean(Is)),
                   ('VI', VI),
                   ('In', In),
                   ]

        # Avoid calculating N when n_mean is <= 5.  This is because it
        # takes too long due to the hard recursive and pair-to-pair
        # calculations it requires.
        if Gs[0].n_mean() > 5 and not getattr(settings, 'no_N', False):
            Nmi = numpy.mean(
                [ util.mutual_information_overlap(Gs[i], Gs[j])
                  for i,j in data['pairIndexes'] ])
        else:
            Nmi = float('nan')
        returns.append(('N',  Nmi))

        return returns
    calcMethods.append(calc_multualInformation)
    def calc_overlap(self, data, settings):
        Gs = data['Gs']
        overlapGs = data.get('ovGs', None)
        if overlapGs is None:
            return [ ]
        # Avoid calculating N when n_mean is <= 5.  This is because it
        # takes too long due to the hard recursive and pair-to-pair
        # calculations it requires (same as above).
        if Gs[0].n_mean() > 5 and not getattr(settings, 'no_N', False):
            NmiO = numpy.mean(
                [ util.mutual_information_overlap(overlapGs[i], overlapGs[j])
                  for i,j in data['pairIndexes'] ])
        else:
            NmiO = float('nan')
        n_mean_ov = sum(G.n_mean() for G in overlapGs)/float(len(overlapGs))

        returns = [('ov_N',      NmiO),
                   ('ov_n_mean', n_mean_ov),
                   ('ov_E',    numpy.mean(tuple(G.energy(data['gamma'])
                                                for G in overlapGs))),
                   ]

        # This records the number of rounds needed to perform the
        # minimization.
        data2 = data.copy()
        data2['Gs'] = data['ovGs']
        ov_ret = self.calc_changes(data=data2, settings=settings)
        returns.extend(('ov_'+name, value)
                       for name, value in ov_ret)

        return returns
    calcMethods.append(calc_overlap)

    #
    # Analysis methods
    #
    def table(self):
        """Compiles data from all MRCs onto self.

        This takes the data from the MRCs on self, and puts it in
        lists by the field names.  This is the first step to plotting,
        writing, etc.  This is a low-cost method.
        """
        with self._lock:
            table= { }
            gammas = numpy.asarray(sorted(self._data.keys()))
            table['gamma'] = gammas
            for fieldname in self.fieldnames:
                # gamma is special cased
                if fieldname == 'gamma':
                    continue
                array = numpy.asarray([self._data[gamma].data[fieldname].mean
                                       for gamma in gammas])
                table[fieldname] = array

        #self.n_hist       = n_hist       = [mrc.n_hist       for mrc in MRCs]
        #self.n_hist_edges = n_hist_edges = [mrc.n_hist_edges for mrc in MRCs]
        return table

    def write(self, fname):
        """Save multi-resolution data to a file."""
        table = self.table()
        f = open(fname, 'w')
        print >> f, "#", time.ctime()
        #print >> f, "# replicas:", self.replicas
        #print >> f, "# minimizer and kwargs:", self._minimizer, \
        #                                       self._minimizerargs

        print >> f, "#", " ".join("%d:%s"%(i+1,x)
                                  for i,x in enumerate(self.fieldnames))
        for i in range(len(table['gamma'])):
            print >> f, table['gamma'][i],
            for name in self.fieldnames:
                if name == "gamma": continue
                print >> f, table[name][i],
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
            ax.set_xlim( 0.1 , self.Nnodes )
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

    def plot(self, fname=None, ax1items=['VI', 'In', 'ov_N'], ax2items=['q'],
             xaxis='gamma',
             plotstyles={}):
        """Plot things from a multiresolution plot.

        ax1items: List of items to plot on axis 1 (left)
        ax1items: List of items to plot on axis 2 (right)

        Each item should be a name of something on the
        MultiResolutionCorrelation object.
        """
        # Recalculate all the stuff
        table = self.table()
        # remove NmiO if we don't have that attribute
        if ax1items==['VI', 'In', 'ov_N'] and 'ov_N' in ax1items:
            if 'ov_N' not in table:
                ax1items = ax1items[:]
                ax1items.remove('ov_N')
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
            '_s_F1':   dict(label='$_sF_1$',color='black', linestyle='-.'),
            'n_mean': dict(label='$<n>$',color='blue',linestyle=':'),
            'ov_n_mean':dict(label='$<n_{ov}>$',color='green',linestyle=':'),
            'entropy':dict(label='$H$',color='blue', linestyle='--'),
            's_F1_ov':dict(label='$F1_{s,ov}$'),
            None:   dict(),  # defaults
            }
        plotstyles = recursive_dict_update(defaultplotstyles, plotstyles)
        # Fill axes
        for item in ax1items:
            plotstyle = plotstyles.get(item, {}).copy()
            scale = plotstyle.pop('scale', 1)
            l = ax1.plot(table['gamma'], table[item]*scale,
                         **plotstyle)[0]
            legenditems.append((l, plotstyle.get('label', item)))
        for item in ax2items:
            plotstyle = plotstyles.get(item, {}).copy()
            scale = plotstyle.pop('scale', 1)
            #plotstyle2 = plotstyle.copy()
            #plotstyle2.pop('scale', None)
            l = ax2.plot(table['gamma'], table[item]*scale,
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
        return f

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





class MRRunner(object):
    """Class to do full multi-resolution analysis.
    """
    _callback = None
    _output = None
    _lock = None
    def __init__(self, MR, #output=None, savefigargs=None,
                 #minimizer='trials',
                 #minimizerargs=dict(minimizer='greedy', trials=10),
                 #plotargs=None,
                 #overlap=False,
                 #pairstyle='all',
                 #settings={}
                 ):
        """Create a multiresolution helper system.

        output: save table of gamma, q, H, etc, to this file.

        minimizer: name of minimizer to use.  This should be a method
        name on the Graph object.

        minimizerargs: keyword arguments (dict) to be passed to the minimizer.

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
        #self._output = output
        #self._savefigargs = savefigargs
        #self._plotargs = plotargs
        #self._minimizer = minimizer
        #self._minimizerargs = minimizerargs
        #self._settings = settings
        self.MR = MR

    def do_gamma(self, gamma, callback=None):
        """Worker thread for gammas.

        For one gamma, look up the respective gamma and do that analysis.
        """
        if callback is None:
            callback = self._callback
        logger.info("Thread %s MR-minimizing g=%f"%(self.thread_id(), gamma))

        data = { }
        state = { }
        data['gamma'] = gamma
        state['gamma'] = gamma
        # Minimize the main systems
        minGs = [ ]
        for G in self._Gs:
            if self._lock:  self._lock.acquire()
            G = G.copy()
            if self._lock:  self._lock.release()
            getattr(G, self.MR.minimizer)(gamma, **self.MR.minimizerargs)
            minGs.append(G)
        data['Gs'] = minGs
        state['Gs'] = [ G.getcmtystate() for G in minGs ]

        # Get the total minimum system:
        Gmin_index, Gmin = min(enumerate(minGs),
                               key=lambda x: x[1].energy(gamma))
        data['Gmin'] = Gmin
        state['Gmin'] = Gmin.getcmtystate()
        state['Gmin_index'] = Gmin_index

        if self.MR.overlap:
            overlapTrials = self.MR.overlap
            overlapMinimizer = getattr(self.MR, 'ovMinimizer', 'ovGreedy')
            overlapGs = [ ]
            for G in minGs:
                G = G.copy()
                G.trials(gamma, trials=overlapTrials, initial='current',
                         minimizer=overlapMinimizer)
                overlapGs.append(G)
                data['ovGs'] = overlapGs
                state['ovGs'] = [ G.getcmtystate() for G in overlapGs ]

                # Get the total minimum system:
                Gmin_index, Gmin = min(enumerate(overlapGs),
                                       key=lambda x: x[1].energy(gamma))
                data['ovGmin'] = Gmin
                state['ovGmin'] = Gmin.getcmtystate()
                state['ovGmin_index'] = Gmin_index

        logger.info("Thread %s MR-minimizing g=%f: done"%(
                                                      self.thread_id(), gamma))
        # Do the information theory VI, MI, In, etc, stuff on our
        # minimized replicas.
        logger.info("Thread %s initializing MRD g=%f"%(
                                              self.thread_id(), gamma))
        self.MR.add(gamma=gamma, data=data, state=state)

        logger.debug("Thread %s initializing MRD g=%f: done"%(
                                              self.thread_id(), gamma))
        # Save output to a file.
        if self.MR.output is not None and self.lockAcquire(blocking=False):
            logger.debug("Thread %s write"%self.thread_id())
            self.MR.write(self.MR.output)
            self.lockRelease()
        # do .savefig() on the minimum E config
        if self.MR.savefigargs is not None \
               and self.lockAcquire(blocking=True):
            logger.debug("Thread %s savefig"%self.thread_id())
            minG = data['Gmin']
            kwargs = self.MR.savefigargs.copy()
            kwargs['fname'] = kwargs['fname']%{'gamma':gamma}
            minG.remap(check=False)
            minG.savefig(**kwargs)
            self.lockRelease()
        # do .plot() on the MultiResolutionData
        if self.MR.plotargs is not None and self.lockAcquire(blocking=False):
            logger.debug("Thread %s plot"%self.thread_id())
            #self.plot(**self._plotargs)
            self.MR.plot(**self.MR.plotargs)
            self.lockRelease()
        # Run callback if we have it.
        if callback:
            logger.debug("Thread %s callback"%self.thread_id())
            callback(gamma=gamma, data=data, state=state,
                     MR=self.MR, MRR=self)
    def _thread(self):
        """Thread worker - do gammas until all are exhausted.

        One of these are spawned for each thread, they run until all
        gammas have been finished."""
        # The queue.join() may return by the time we get to the
        # except: clause, where self._EmptyException has already been
        # deleted.  Thus, we need to make a local copy before the
        # try/except block.
        while True:
            gamma = self._getGamma()
            if gamma is None:
                return
            logger.info("Thread %s begin gamma=%f"%(
                self.thread_id(), gamma))
            self.do_gamma(gamma)
            logger.info("Thread %s end gamma=%f"%(
                self.thread_id(), gamma))
        logger.info("Thread terminated: %s", self.thread_id())

    def thread_id(self):
        import threading
        return threading.current_thread().name
    def lockAcquire(self, name='writelock', blocking=False):
        """Acquire lock self._NAME.

        blocking: if true, block for lock to acquire.  If false, don't
        block, return fals if lock not acquired."""
        # Someday both of these could be made into context managers -
        # the only reason they aren't is that sometimes, we *won't*
        # use the locks (and need to abstract out blocking or not)...
        logger.debug('Thread %s attempting lock acquisition: %s',
                     self.thread_id(), name)
        if hasattr(self, '_'+name):
            locked = getattr(self, '_'+name).acquire(blocking)
            logger.debug('Thread %s attempting lock acquisition: %s: %s',
                         self.thread_id(), name, locked)
            return locked
        # If we don't have self._LOCKNAME, we don't lock so always return true.
        return True
    def lockRelease(self, name='writelock'):
        """Release lock self.NAME."""
        logger.debug('Thread %s attempting release: %s',
                     self.thread_id(), name)
        if hasattr(self, '_'+name):
            logger.debug('Thread %s attempting release (phase 2): %s',
                         self.thread_id(), name)
            return getattr(self, '_'+name).release()
    def _getGamma(self):
        self.lockAcquire(name='gammalock', blocking=True)
        if not hasattr(self, '_seenGammas'):
            self._seenIndexes = set()
            self._seenGammas = set()
        if isinstance(self._gammas, (list, tuple)):
            gamma = self._getGamma_list()
        else:
            gamma = self._getGamma_dict()
        self._seenGammas.add(gamma)
        self.lockRelease(name='gammalock')
        #from fitz import interactnow
        return gamma

    def _getGamma_list(self):
        """Return next gamma, if self._gammas is a list"""
        if isinstance(self._gammas, tuple):
            self._gammas = list(self._gammas)
        if len(self._gammas) == 0:
            return None
        i = len(self._gammas) // 2
        return self._gammas.pop(i)

    def _getGamma_dict(self):
        """Return next gamma, if self._gammas is a dict for LogInterval"""
        p = self._gammas
        if not hasattr(self, '_logGammas'):
            args = p.copy()
            if args.get('high') == 'auto':  args.pop('high')
            if args.get('low')  == 'auto':  args.pop('low')
            if args.get('start'):           args.pop('start')
            self._logGammas = LogInterval(**args)
        logGammas = self._logGammas

        start = p.get('start', 1)
        if isinstance(p.get('high'), (int,float)) and p.get('high', 1) < start:
            start = p['high']
        if isinstance(p.get('low'), (int,float)) and p.get('low', 1) > start:
            start = p['low']
        high = self._gammas.get('high', 'auto')
        low  = self._gammas.get('low',  'auto')

        i_start = logGammas.index(start)
        i_high, i_low = 'auto', 'auto'
        if high!='auto': i_high = logGammas.index(high)
        if low !='auto': i_low  = logGammas.index(low)

        i_next = self._getGamma_dict_nextindex(i_start, i_low, i_high)
        if i_next is None:
            return None
        self._seenIndexes.add(i_next)
        next_gamma = self._logGammas.value(i_next)
        return next_gamma

    def _getGamma_dict_nextindex(self, i_start, i_low, i_high):
        # try low:
        def trylow():
            """Return the next lowest gamma if there are any left."""
            next_i_low = min(self._seenIndexes) - 1
            next_gamma = self._logGammas.value(next_i_low)

            # If we are *not* automatic.
            if i_low != 'auto':
                if next_i_low <= i_low:
                    return None
                return next_i_low

            # Automatic detection after this point:

            # If we've already run that gamma, then always run it
            # again (for multiplet MRRunner runs):
            if next_gamma in self.MR._data:
                return next_i_low
            # Emergency override, don't allow auto to go less
            # than .0001 (since for disconnected graphs, we
            # won't every be able to get to q == 1.  FIXME:
            # improve sometime.
            if min(self._seenIndexes) <= self._logGammas.index(.000001):
                return None
            doneGammas = sorted(self._seenGammas)
            qs = [ self.MR._data[g].data['q'].mean
                   for g in doneGammas]
            n_means = [ self.MR._data[g].data['n_mean'].mean
                        for g in doneGammas ]
            #lowestgamma = min(self.MR._data.keys())
            # First see if num clusters == 1, if so we are done
            if qs[0] == 1:
                return None
            elif (len(doneGammas) > 10
                and all(q == qs[0] for q in qs[:10])
                and all(nm==n_means[0] for nm in n_means[:10])
                ):
                return None
            # Next see if cluster size == N, if so we are done
            elif (self.MR._data[doneGammas[0]].data['n_mean'].mean
                  == self.MR._data[doneGammas[0]].data['nodes'].mean ):
                return None

            return next_i_low

        def tryhigh():
            """Return the next highest gamma if there are any more left."""
            # This function has problems, since _i_high should be
            # i_high, but that causes global/local scope problems
            # for i_high.
            next_i_high = max(self._seenIndexes) + 1
            next_gamma = self._logGammas.value(next_i_high)
            # Manual detection
            if i_high != 'auto':
                if next_i_high >= i_high:
                    return None
                return next_i_high
            # Automatic detection after this point:
            # This isn't intelligently detected yet, so just
            # default to 1000.  FIXME: improve sometime.
            default_i_high = self._logGammas.index(1000)
            if next_i_high >= default_i_high:
                return None
            return next_i_high


        # Do the actual logic:
        # If we haven't done anything at all, start with 'start'
        if len(self._seenIndexes) == 0:
            i_start
            return i_start
        if max(self._seenIndexes)-i_start > i_start-min(self._seenIndexes):
            # next_i_high farther away than next_i_low
            next_i = trylow()
            if next_i is not None:
                return next_i
            return tryhigh()
        else:
            next_i = tryhigh()
            if next_i is not None:
                return next_i
            return trylow()

    def do(self, Gs, gammas=None,
           threads=1, callback=None):
        """Do multi-resolution analysis on replicas Gs with `trials` each."""
        # If multiple threads requested, do that version instead:
        #if gammas is None and logGammaArgs is not None:
        #    gammas = LogInterval(**logGammaArgs).values()
        #if gammas is None:
        #    raise ValueError("Either gammas or logGammaArgs must be given.")
        self._gammas = gammas
        if threads > 1:
            return self.do_mt(Gs, gammas, threads=threads,
                              callback=callback)

        # Non threaded version:
        self._callback = callback
        self._Gs = Gs
        self.replicas = len(Gs)
        while True:
            gamma = self._getGamma()
            if gamma is None:
                break
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
        self._lock = threading.Lock()
        self._gammalock = threading.Lock()
        self._writelock = threading.Lock()

        # Start the threads
        threadlst = [ ]
        for tid in range(threads):
            t = threading.Thread(target=self._thread)
            t.daemon = True
            t.start()
            threadlst.append(t)

        # Wait for all threads to join:
        for t in threadlst:
            t.join()

        #self._queue.join()
        del (self._Gs, self._gammalock, self._lock, self._writelock, )
        del self._callback # unmasks class object's _callback=None



def write_table(fname, MRs, values, Vname='V'):
    f = open(fname, 'w')

    print >> f, "#", " ".join("%d:%s"%(i+1,x)
                     for i,x in enumerate((Vname,)+tuple(MRs[0].fieldnames)))

    for v, MR in zip(values, MRs):
        table = MR.table()

        for i in range(len(table['gamma'])):
            print >> f, v, table['gamma'][i],
            for name in MR.fieldnames:
                if name == "gamma": continue
                print >> f, table[name][i],
            print >> f

        print >> f

