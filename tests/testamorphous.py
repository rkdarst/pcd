# Richard Darst/Glen Hocky, August 2011

import os
import sys

import numpy

import pcd
import graphs
fast = globals().get('fast', False)


# All of the loading stuff has been moved to support.gromacs
G = graphs.bss2d_n240()

dirname = 'tests-output/amorphous/'
if not os.path.exists(dirname):
    os.makedirs(dirname)

G.minimize(gamma=0)
G.savefig(os.path.join(dirname, 'amorphous_gamma0.svg'),
          radii=G.radii,base_radius=0.4)

def callback(G, gamma, **kwargs):
    G.remapCommunities(check=False)
    fname = os.path.join(dirname, 'imgs/amorphous_gamma%011.5f.svg'%gamma)
    G.savefig(fname,radii=G.radii,base_radius=0.4)

nreplicas, ntrials = 10, 20
if fast:
    nreplicas, ntrials = 5, 5

from pcd import LogInterval
gammas = LogInterval(.001, 100, density=10).values()
MR = pcd.MultiResolution()
MR.do([G]*nreplicas, gammas, trials=ntrials, threads=1,
      #callback=callback
      )
#MR = pcd.MultiResolution(0.1, 1, callback=callback, number=10)
#MR.do([G]*10, trials=10, threads=2)
MR.calc()
MR.write(os.path.join(dirname, "test_amorphous_values.txt"))
#MR.plot(os.path.join(dirname, amorphous_mr.png"))
#MR.plot_nmean(os.path.join(dirname, "amorphous_mr_nmean.png"))
#MR.plot_nhists(os.path.join(dirname, "nhists.pdf"))
