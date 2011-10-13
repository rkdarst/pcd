# Richard Darst, October 2011

import collections
import itertools
import math
from math import exp, log, floor, ceil
import numpy
import random

import networkx

import cmodels
import util

log2 = lambda x: math.log(x, 2)

class _GraphAnneal(object):
    """Helper class to hold anneal-related methods
    """

    def anneal(self, gamma, Escale=None, stablelength=10,
               attempts=1000, betafactor=1.001):
        """An annealing round

        Escale: initial energy scale.  Set this high enough for your
        system in order to get good minimization.

        attempts: self.N*attempts is the number of moves we attempt on
        each cycle.

        stablelength: Run until we go through this many cycles with
        all changes rejected.
        """
        if Escale is None:
            beta = 1./self._find_average_E(gamma)
        else:
            beta = 1./Escale
        running_changes = collections.deque((None,)*stablelength)

        for i in itertools.count():
            changes = self._anneal(gamma, beta=beta, steps=self.N*attempts)
            E = self.energy(gamma)
            if self.verbosity >= 2:
                print "  (a%6d) %9.5f %7.2f %6d"%(i, beta, E, changes)
            beta *= betafactor

            running_changes.append(changes)
            del running_changes[0]
            if all(x==0 for x in running_changes):
                break
        return i


    def _anneal(self, gamma, beta, steps=None, deltabeta=0):
        if steps is None:
            steps = self.N*1000
        return cmodels.anneal(self._struct_p, gamma, beta,
                              steps, deltabeta)

    def _find_average_E(self, gamma, trials=None):
        trials = self.N * 100
        Es = [ ]
        for i in range(trials):
            n = int(floor(random.uniform(0, self.N)))
            c = random.choice(tuple(self.cmtys()))
            Es.append(self.energy_cmty_n(gamma, c, n))
        #return numpy.mean(Es)
        #print Es
        return max(abs(x) for x in Es)
