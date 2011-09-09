# Richard Darst/Glen Hocky, August 2011

import os
import sys

import numpy

import models

def readgro(gro_file,dimensions=2):
    configuration_lines=open(gro_file,'r').readlines()
    natoms=int(configuration_lines[1])
    configuration=numpy.array( [ line.split()[3:3+dimensions] for line in configuration_lines[2:2+natoms] ],dtype=float)
    box=[float(x) for x in configuration_lines[-1].split()]
    atomtypes=[line.split()[1] for line in configuration_lines[2:2+natoms] ]
    ntypea=atomtypes.index(atomtypes[-1])
    length=box[0]
    return length,ntypea,configuration

def e_r6_i(d,sigma,epsilon=1):
    energy = -epsilon*(sigma/d)**6
    energy += 0.5*epsilon
    energy *= 1000
    intn = numpy.round(energy)
    return intn

def e_softsphere_i(d,sigma,epsilon=1):
    energy = epsilon*(sigma/d)**12
    energy += 0.5*epsilon
    energy *= 1000
    intn = numpy.round(energy)
    return intn

def e_lj_i(d,sigma,epsilon=1):
    energy = 4*epsilon*( (sigma/d)**12 - (sigma/d)**6 )
    energy += 0.5*epsilon
    energy *= 1000
    intn = numpy.round(energy)
    return intn

def get_imatrix(coords,sigmas,ntypea,periodic=None):
    imatrix=numpy.zeros((len(coords),len(coords)),dtype=int)
    orig_settings = numpy.seterr()
    numpy.seterr(all="ignore")
    for n1 in range(len(coords)):
        atomtype=int(n1>=ntypea) #typa a = 0 type b = 1
        delta = coords[n1] - coords
        if periodic:
            delta -=  numpy.round(delta/float(periodic))*periodic
        delta = delta**2
        dist = numpy.sum(delta, axis=1)
        dist = numpy.sqrt(dist)
        sigma=sigmas[atomtype] 
#        e = e_lj_i(dist,sigma)
        e = e_r6_i(dist,sigma)
        imatrix[n1, :] = e
    numpy.seterr(**orig_settings)
    return imatrix

testsdir=os.path.dirname(sys.argv[0])
datadir=os.path.abspath(os.path.join(testsdir,'../data'))
grofile=os.path.join(datadir,'2dss32_n240_T0.5_1.gro')
L,ntypea,coords=readgro(grofile)
ntypeb=len(coords)-ntypea

sigma1=numpy.concatenate(( numpy.ones(ntypea),1.1*numpy.ones(ntypeb) ))
sigma2=numpy.concatenate(( 1.1*numpy.ones(ntypea),1.4*numpy.ones(ntypeb) ))
sigmas=numpy.array((sigma1,sigma2))
radii=numpy.concatenate(( numpy.ones(ntypea),1.4*numpy.ones(ntypeb) ))

imatrix=get_imatrix(coords,sigmas,ntypea,periodic=L)

G = models.Graph(N=len(coords))
G = G.from_imatrix(imatrix)

if not os.path.exists('imgs'):
    print "Creating dir:",os.path.join(os.getcwd(),'imgs')
    os.mkdir('imgs')

G.minimize(gamma=0)
G.savefig('imgs/amorphous_gamma0.svg', coords=coords,radii=radii,base_radius=0.4)

def callback(G, gamma, **kwargs):
    G.remapCommunities(check=False)
    fname = 'imgs/amorphous_gamma%011.5f.svg'%gamma
    G.savefig(fname, coords=coords,radii=radii,base_radius=0.4)

MR = models.MultiResolution(.001, 100, callback=callback, number=10)
MR.do([G]*10, trials=20, threads=2)
MR.calc()
MR.plot("imgs/amorphous_mr.png")
