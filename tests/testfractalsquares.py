# Richard Darst/Glen Hocky, August 2011

import os
import sys

import numpy

import pcd

from support.fractalsquare import fractalsquare1, fractalsquare2


def e_lj(d):
    d *= 2**(1/6.)
    energy = 4 * (1/d**12. - 1/d**6.)
    energy *= 1000
    energy += 1
    intn = numpy.round(energy)
    return intn
def e_inverse_6(d):
    energy = -(1./d**6)
    energy *= 729
    energy += 10
    intn = numpy.round(energy)
    return intn

# Mode 1:
#coords = glencoords()
coords, L = fractalsquare1(4)
G = pcd.Graph.from_coords_and_efunc(coords,
                                     #lambda x: e_lj(x)*100,
                                     e_lj,
                                     boxsize=L)
G = pcd.Graph(N=len(coords))

# Mode 2
#sep = 16
imatrix, coords = fractalsquare2(L=16)
print imatrix
print coords
G = G.from_imatrix(imatrix)

G.minimize(gamma=0)
G.savefig('img.png', coords=coords)

if not os.path.exists('imgs'):
    print "Creating dir:",os.path.join(os.getcwd(),'imgs')
    os.mkdir('imgs')

def callback(G, gamma, **kwargs):
    G.remapCommunities(check=False)
    fname = 'imgs/gamma%011.5f.png'%gamma
    Es = [ G.energy_n(gamma, n) for n in range(G.N) ]
    G.savefig(fname, coords=coords, energies=Es,
              nodes='circles')

nreplicas, ntrials, ldensity = 10, 250, 20
if globals().get('fast', True):
    nreplicas, ntrials, ldensity = 5, 5, 10

MR = pcd.MultiResolution()
MR.do([G]*nreplicas, logGammaArgs=dict(low=.001, high=100, density=ldensity),
      trials=ntrials, threads=2, callback=callback)
MR.calc()
MR.plot("imgs.png")
