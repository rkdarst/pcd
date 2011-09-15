# Glen Hocky, September 2011
# Richard Darst, September 2011

import numpy

from math import log, exp, floor, ceil
log2 = lambda x: log(x, 2)
exp2 = lambda x: 2.**x

# This is by Glen Hocky
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


# Below here is by Richard Darst
def fractalsquare1(sep):
    """Make a fractal square with a side length of 2^sep
    """
    if sep == 0:
        return numpy.asarray([[0,0]]), 1
    else:
        sq, L = fractalsquare1(sep=sep-1)
        # Offset it by size of smaller square plus the new spacing
        offset = L #max(c[0] for c in sq) + sep
        L = L + L
        new = numpy.concatenate((sq,
                                 sq+(0, offset),
                                 sq+(offset, 0),
                                 sq+(offset, offset),
                                 ))
        return new, L+1



def tick(x):
    s = log2(x+1)
    if s == floor(s):
        return int(s)
    return tick(x - exp2(floor(s)))
#for x in range(0, 16):
#    print x, tick(x)

def fractalsquare2(L):
    N = L*L
    a = numpy.arange(N)
    a.shape = L,L

    imatrix = numpy.zeros(shape=(N, N))
    imatrix[:] = 1

    coords = numpy.arange(N)
    coords = numpy.concatenate(((coords//L,), (coords%L,)),
                               axis=0)
    coords = coords.transpose()

    lowest_t = int(ceil(log2(L)))
    lowest_V = int(ceil(log2(L)))+1

    def efunc(t):
        return -lowest_V + (t+1)
        #return (-lowest_V + (t+1)) * 3
        #return -2**(lowest_t+1)  +  2**(t+1)  -  1

    for i in range(L*L):
        y, x = divmod(i, L)
#        coords.append((x,y))

        # left:
        V = efunc(tick(x))
        xprime = (x + 1) % L
        iprime = y*L + xprime
        imatrix[i, iprime] = V
        imatrix[iprime, i] = V
        print x, y, V, xprime, y, i, iprime

        # bottom:
        V = efunc(tick(y))
        yprime = ((y + 1)%L)
        iprime = yprime*L + x
#        print x, y, V, x, yprime, i, iprime
        imatrix[i, iprime] = V
        imatrix[iprime, i] = V

    return imatrix, numpy.asarray(coords)

if __name__ == "__main__":
    import pylab
    fig=pylab.figure()
    ax = fig.add_subplot(111, aspect='equal')

    sq, L = fractalsquare1(6)
    #print sq

    pylab.scatter(sq[:,0], sq[:, 1])
    #pylab.show()

    imatrix, coords = fractalsquare2(4)
    print imatrix
    print coords
