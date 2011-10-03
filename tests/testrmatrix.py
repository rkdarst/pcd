# Richard Darst, October 2011

import os
import sys

import numpy

import pcd
from pcd.support.fractalsquare import fractalsquare1, fractalsquare2


def e_inverse_6(d):
    energy = -(1./d**6)
    energy *= 729
    #energy += 10
    return energy
    #intn = numpy.round(energy)
    #return intn


sep = 4
coords, L = fractalsquare1(sep=sep)
# refunc is the rmatrix efunc.  Here, it is a contstant so we can just
# make it trivial, otherwise it runs the same as the "efunc" argument
# to this method.
G = pcd.Graph.from_coords_and_efunc(coords, e_inverse_6,
                                    refunc=lambda d: 10)


if not os.path.exists('tests-output/rmatrix'):
    os.makedirs('tests-output/rmatrix')

def callback(G, gamma, **kwargs):
    G.remapCommunities(check=False)
    fname = 'tests-output/rmatrix/gamma%09.4f.png'%gamma
    Es = [ G.energy_n(gamma, n) for n in range(G.N) ]
    G.savefig(fname, coords=coords, energies=Es,
              nodes='circles')

nreplicas, ntrials, ldensity = 10, 50, 10
if globals().get('fast', False):
    nreplicas, ntrials, ldensity = 5, 5, 10

MR = pcd.MultiResolution()
MR.do([G]*nreplicas, logGammaArgs=dict(low=.00001, high=100, density=ldensity),
      trials=ntrials, threads=2, callback=callback)
MR.calc()
MR.plot("tests-output/rmatrix/multiresolution.png")
