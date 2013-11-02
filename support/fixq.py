# Richard Darst, October 2013.

# This module takes a partition and merges/splits further things until
# it achieves the correct number of q.

from pcd import anneal
from pcd import cmodels

import itertools

def fixq_shift_nodes(G, gamma):
    """Shift individual nodes to their best communities, keeping q
    constant throughout."""
    G.const_q = G.q
    changes_cumul = 0
    #changes = G._anneal(gamma, steps=G.N*10, beta=1e3, steps=100,
    #                    new_cmty_prob=0, p_binary=0, p_collective=0)
    while True:
        changes = cmodels.greedy(G, gamma)
        changes_cumul += changes
        #print 'greedy round', changes
        if changes == 0:
            break
    G.const_q = False
    return changes_cumul

def fixq_merge(G, gamma):
    return cmodels.fixq_merge(G, gamma)
def fixq_merge_density(G, gamma):
    return cmodels.fixq_merge_density(G, gamma)

import pcd.support.algorithms

import pcd.support.spectral
import pcd
class FixQ(pcd.support.algorithms.CDMethod):
    _input_format = 'null'
    true_q = None
    q_hint = None
    gamma = 1.0
    parent_method = pcd.support.algorithms.Louvain
    max_louvain_rounds = 10000
    Emode = 'mod'
    merge_mode = 'E'
    shift_every_round = True
    def run(self):
        g = self.g

        true_q = self.true_q
        if true_q is None:
            spec = pcd.support.spectral.Flow(g, num_ev=self.q_hint+1)
            true_q = spec.q()
            self.time_ev = spec.time_ev

        for i in range(self.max_louvain_rounds):
            lv = self.parent_method(g, verbosity=-1)
            print "%s result sizes:"%self.parent_method.name(), \
                  [ r.q for r in lv.results ], ", wanted size %s"%true_q
            # lv.results should have greatest q first.
            results = [ r for r in lv.results if r.q >= true_q ]
            if len(results) > 0:
                break
        self.num_parent_method_trials = i+1
        if len(results) < 1:
            raise ValueError("All %s results return fewer than q communities."%self.parent_method.name())
        cmtys = results[-1]

        G = pcd.Graph.fromNetworkX(g)
        cmtys.load_pcd(G)
        #from fitz import interactnow
        if self.Emode == 'mod':
            G.enableModularity(mode=None)
        elif self.Emode == 'APM':
            pass
        else:        raise ValueError("self.Emode=%s not recognized"%self.Emode)

        self.fixq(G, gamma=self.gamma, q=true_q)

        self.cmtys = cmtys.from_pcd(G)
        self.results = [ self.cmtys ]
    def read_cmtys(self):
        pass

    def fixq(self, G, gamma, q):
        assert q is not None, "q is none"
        assert G.q >= q, "G.q is too low: G.q=%s, q=%s"%(G.q, q)

        for i in itertools.count():
            if G.q <= q:
                break

            # Merge community pairs, using either E or edge density as a
            # criteria.
            if self.merge_mode == 'E':
                deltaE = fixq_merge(G, gamma)
            elif self.merge_mode == 'density':
                deltaE = fixq_merge_density(G, gamma)
            #print "merged two communities:", deltaE

            # Take each node, in random order, and move it to the adjacent
            # community that it has the lowest energy with.
            changes_cumul = None
            if self.shift_every_round:
                changes_cumul = fixq_shift_nodes(G, gamma=gamma)

            print i, G.q, q, deltaE, changes_cumul
        # Shift nodes in place if we don't do it after every merge round.
        if not self.shift_every_round:
            changes_cumul = fixq_shift_nodes(G, gamma=gamma)
        assert G.q == q

        # No return, mutate G in place
