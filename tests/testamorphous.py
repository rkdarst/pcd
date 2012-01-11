# Richard Darst/Glen Hocky, August 2011

import os
import cPickle as pickle
import sys

import numpy

import pcd
import graphs
fast = globals().get('fast', False)


# All of the loading stuff has been moved to support.gromacs
G = graphs.bss2d_n240_T050()
#G = graphs.bss2d_n5760_T040()
#G = graphs.bss2d_n23040_T040()
#pickle.dump(G, open('G.pickle', 'w'), pickle.HIGHEST_PROTOCOL) ; exit()
#G = pickle.load(open('G.pickle'))

dirname = 'tests-output/amorphous/'
if not os.path.exists(dirname):
    os.makedirs(dirname)

G.trials(gamma=0, trials=10)
#exit()
G.savefig(os.path.join(dirname, 'amorphous_gamma0.png'),
          radii=G.radii,base_radius=0.4)

#def callback(G, gamma, **kwargs):
#    G.remapCommunities(check=False)
#    fname = os.path.join(dirname, 'amorphous_gamma%011.5f.svg'%gamma)
#    G.savefig(fname,radii=G.radii,base_radius=0.4)

#nreplicas, ntrials = 10, 20
nreplicas, ntrials = 3, 5
if fast:
    nreplicas, ntrials = 5, 5

from pcd import LogInterval
gammas = LogInterval(.001, 100, density=10).values()

MR = pcd.MultiResolution(
    savefigargs=dict(
        fname=os.path.join(dirname, 'amorphous_gamma%(gamma)011.5f.svg'),
        base_radius=0.4),
    minimizerargs=dict(trials=ntrials),
    )
MR.do([G]*nreplicas, gammas, threads=1,
      )
#MR = pcd.MultiResolution(0.1, 1, callback=callback, number=10)
#MR.do([G]*10, trials=10, threads=2)
MR.calc()
MR.write(os.path.join(dirname, "mr.txt"))
#MR.plot(os.path.join(dirname, amorphous_mr.png"))
#MR.plot_nmean(os.path.join(dirname, "amorphous_mr_nmean.png"))
#MR.plot_nhists(os.path.join(dirname, "nhists.pdf"))
