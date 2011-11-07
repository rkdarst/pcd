# Glen Hocky, September 2011

# Note: the actual get_test_G could be moved to graphs.py

import numpy
import os

import pcd

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


def readgro(gro_file,dimensions=2):
    configuration_lines=open(gro_file,'r').readlines()
    natoms=int(configuration_lines[1])
    configuration=numpy.array( [ line.split()[3:3+dimensions] for line in configuration_lines[2:2+natoms] ],dtype=float)
    box=[float(x) for x in configuration_lines[-1].split()]
    atomtypes=[line.split()[1] for line in configuration_lines[2:2+natoms] ]
    ntypea=atomtypes.index(atomtypes[-1])
    length=box[0]
    return length,ntypea,configuration

def get_imatrix(coords,sigmas, atomtypes, boxsize=None):
    imatrix=numpy.zeros((len(coords),len(coords)),dtype=int)
    orig_settings = numpy.seterr()
    numpy.seterr(all="ignore")
    for n1 in range(len(coords)):
        #atomtype=int(n1>=ntypea) #typa a = 0 type b = 1
        atomtype = atomtypes[n1] #typa a = 0 type b = 1
        delta = coords[n1] - coords
        if boxsize:
            delta -=  numpy.round(delta/float(boxsize))*boxsize
        delta = delta**2
        dist = numpy.sum(delta, axis=1)
        dist = numpy.sqrt(dist)
        sigma=sigmas[atomtype]
#        e = e_lj_i(dist,sigma)
        e = e_r6_i(dist,sigma)
        imatrix[n1, :] = e
    numpy.seterr(**orig_settings)
    return imatrix


def load_bss2d(fname='2dss32_n240_T0.5_1.gro'):
    if '/' not in fname:
        testsdir = os.path.dirname(__file__)
        datadir = (os.path.join(testsdir,'./../data'))
        grofile = os.path.join(datadir,fname)
    else:
        grofile = fname
    L,ntypea,coords=readgro(grofile)
    ntypeb=len(coords)-ntypea

    atomtypes = numpy.concatenate(((0,)*ntypea, (1,)*ntypeb ))

    sigma1=numpy.concatenate(( numpy.ones(ntypea),1.1*numpy.ones(ntypeb) ))
    sigma2=numpy.concatenate(( 1.1*numpy.ones(ntypea),1.4*numpy.ones(ntypeb) ))
    sigmas=numpy.array((sigma1,sigma2))
    radii=numpy.concatenate(( numpy.ones(ntypea),1.4*numpy.ones(ntypeb) ))

    imatrix=get_imatrix(coords,sigmas,atomtypes,boxsize=L)

    G = pcd.Graph.from_imatrix(imatrix)
    G.coords = coords
    G.boxsize = L
    G.radii = radii
    return G


def load_pysim(fname="2dss_n23040_T0.40.atomdump", openargs={}):
    if '/' not in fname:
        testsdir = os.path.dirname(__file__)
        datadir = (os.path.join(testsdir,'./../data'))
        fname = os.path.join(datadir,fname)

    import pysim.models
    frames = pysim.models.openFname(fname, **openargs)
    frame = frames[0]
    N = frame.N
    coords = numpy.asarray((frame.pos[0], frame.pos[1])).transpose()
    boxsize = frame.L

    #from fitz import interactnow


    atomtypes = frame.atomtype
    sigmas = numpy.zeros(shape=(2,N), dtype=float)
    sigmas[0, atomtypes==1] = 1.0
    sigmas[0, atomtypes==2] = 1.1
    sigmas[1, atomtypes==1] = 1.1
    sigmas[1, atomtypes==2] = 1.4
    radii = numpy.zeros(N, dtype=float)
    radii[atomtypes==1] = 1.0
    radii[atomtypes==2] = 1.4
    # Check it:
    assert (sigmas != 0).all()
    assert (radii  != 0).all()

    imatrix = get_imatrix(coords, sigmas, boxsize=boxsize,
                          atomtypes=atomtypes-1) # reindex atomtypes, -1


    G = pcd.Graph.from_imatrix(imatrix)
    G.coords = coords
    G.boxsize = boxsize
    G.radii = radii
    return G
