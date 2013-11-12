# Richard Darst, October 2011

import collections
import ctypes
import itertools
import math
from math import exp, log, floor, ceil
import numpy
import random

import networkx

import cmodels
import util

log2 = lambda x: math.log(x, 2)
def approxeq(a, b):
    if a is None or b is None: return False
    if a == 0 and b == 0: return abs(a-b) < 1e-6
    return a==b or 2*abs(a-b)/float(max(abs(a),abs(b))) < 1e-6

class _GraphAnneal(object):
    """Helper class to hold anneal-related methods
    """

    def anneal(self, gamma, Escale=None, stablelength=10,
               attempts=1, betafactor=1.001,
               new_cmty_prob=.01, p_binary=.05, p_collective=.05,
               const_q_SA=False, constqEcoupling=1.,
               min_n_SA=False, minnEcoupling=.1,
               mode=None,
               maxrounds=None):
        """An annealing round

        Escale: initial energy scale.  Set this high enough for your
        system in order to get good minimization.

        attempts: self.N*attempts is the number of moves we attempt on
        each cycle.

        stablelength: Run until we go through this many cycles with
        all changes rejected.
        """
        if Escale is None:
            beta = .1/self._find_average_E(gamma)
        else:
            beta = 1./Escale
        totChanges = 0
        running_changes = collections.deque((None,)*stablelength)

        steps = self.N*attempts

        if mode == 'guimera':
            steps = self.N**2 + self.N
            new_cmty_prob = .001
            p_collective = 0.0
            p_binary = self.N / float(steps)
            print "steps=%s,p_bin=%s"%(steps, p_binary)
        elif mode:
            raise ValueError("Unknown mode: %s"%mode)

        best_E = float('inf')
        best_state = None

        if self.verbosity >= 2:
            E = self.energy(gamma)
            print "  (a------) %9.5fb %7.2fE ------ %4dq"%(
                beta, E, self.q)

        for nRounds in itertools.count():
            move_info = self._anneal(gamma, beta=beta, steps=steps,
                                   #new_cmty_prob=.01, combine_prob=.01
                                   new_cmty_prob=new_cmty_prob, p_binary=p_binary,
                                   p_collective=p_collective,
                                   const_q_SA=const_q_SA, constqEcoupling=constqEcoupling,
                                   min_n_SA=min_n_SA, minnEcoupling=minnEcoupling,
                                   )
            changes = sum(move_info[1::3])
            self.remap()
            totChanges += changes
            E = self.energy(gamma)
            if nRounds > 10 or nRounds%10 == 0:
                if self.verbosity >= 2:
                    print "  (a%6d) %9.3eb %7.2fE %6dc %4dq"%(
                        nRounds, beta, E, changes, self.q),
                    # verbose status info
                    print ":: "+"  ".join( (
                        "%4d %8.2e %9.2e"%(move_info[i*3+0],
                                           move_info[i*3+1]/move_info[i*3+0],
                                           move_info[i*3+2]/move_info[i*3+0])
                        for i in range(3)) )

            beta *= betafactor
            if E < best_E:
                best_state = self.getcmtystate()

            running_changes.append((changes, E))
            del running_changes[0]
            if nRounds >= stablelength-1:
                # If no more changes
                #if all(x[0]==0 for x in running_changes):
                #    break
                # if E doesn't change for the last stable_length rounds.
                if all(approxeq(x[1],running_changes[0][1]) for x in running_changes):
                    break
            if maxrounds and nRounds > maxrounds:
                break
        self.setcmtystate(best_state)
        return (nRounds, totChanges)


    def _anneal(self, gamma, beta, steps=None, deltabeta=0,
                new_cmty_prob=.01, p_binary=.05, p_collective=.05,
                const_q_SA=False, constqEcoupling=1.,
                min_n_SA=False, minnEcoupling=.1,
                ):
        if steps is None:
            steps = self.N*1000
        #new_cmty_prob = .01 # P of shifting to new community.  .01 is former default.
        #combine_prob = .01  # prob to try *either* split or combine move.  .01 is former default.
        move_info = (ctypes.c_double*12)()
        changes = cmodels.anneal(self._struct_p, gamma, beta,
                              steps, deltabeta,
                              new_cmty_prob, p_binary, p_collective,
                              const_q_SA, constqEcoupling,
                              min_n_SA, minnEcoupling,
                                 move_info,
                              )
        move_info = numpy.array(move_info, dtype=float)
        return move_info

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
