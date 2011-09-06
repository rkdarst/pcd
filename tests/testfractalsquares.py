# Richard Darst/Glen Hocky, August 2011

import os

import numpy

import models

def glencoords():
    # By Glen Hocky, August 2011

    natoms=64
    #import pylab
    #fig=pylab.figure()
    #ax = fig.add_subplot(111, aspect='equal')
    positions=[[i,j] for i in range(natoms) for j in range(natoms)]
    values=[[i,j] for i in range(natoms) for j in range(natoms)]
    for idx,positions in enumerate(positions):
        for i in range(1,natoms+2,2):
            if positions[0] > i: values[idx][0]=values[idx][0]+1
            if positions[1] > i: values[idx][1]=values[idx][1]+1
        for i in range(3,natoms+4,4):
            if positions[0] > i: values[idx][0]=values[idx][0]+2
            if positions[1] > i: values[idx][1]=values[idx][1]+2
        for i in range(7,natoms+8,8):
            if positions[0] > i: values[idx][0]=values[idx][0]+3
            if positions[1] > i: values[idx][1]=values[idx][1]+3
        for i in range(15,natoms+16,16):
            if positions[0] > i: values[idx][0]=values[idx][0]+4
            if positions[1] > i: values[idx][1]=values[idx][1]+4

    xpoints=[x[0] for x in values]
    ypoints=[x[1] for x in values]
    #pylab.scatter(xpoints,ypoints)
    #pylab.show()
    return numpy.asarray(values)

def fractalsquare(sep):
    """Make a fractal square with parameter `sep`

    Number of nodes per side is 2^sep.

    Largest separation between nodes is sep.
    """
    if sep == 0:
        return numpy.asarray([[0,0]]), 1
    else:
        sq, L = fractalsquare(sep=sep-1)
        # Offset it by size of smaller square plus the new spacing
        offset = L
        L = L + L
        new = numpy.concatenate((sq,
                                 sq+(0, offset),
                                 sq+(offset, 0),
                                 sq+(offset, offset),
                                 ))
        return new, L+1


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


def set_imatrix(imatrix, coords, periodic=None):
    for n1 in range(len(coords)):
        delta = coords[n1] - coords
        if periodic:
            delta -=  numpy.round(delta/float(periodic))*periodic
        delta = delta**2
        dist = numpy.sum(delta, axis=1)
        dist = numpy.sqrt(dist)
        e = e_lj(dist)
        imatrix[n1, :] = e * 100


#coords = glencoords()
coords, L = fractalsquare(4)

G = models.Graph(N=len(coords))

set_imatrix(G.interactions, coords, periodic=L)

import matplotlib.cm as cm
import matplotlib.figure
import matplotlib.backends.backend_agg

def callback(G, gamma, **kwargs):
    G.remapCommunities(check=False)
    fname = 'imgs/gamma%07.3f.png'%gamma

    G.savefig(fname, coords=coords)

MR = models.MultiResolution(.00001, 10000, callback=callback)
MR.do([G]*10, trials=10, threads=2)
MR.calc()
MR.plot("imgs.png")
