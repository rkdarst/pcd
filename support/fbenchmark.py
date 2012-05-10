
import os
from os.path import join
import shutil
import subprocess
import tempfile
import textwrap
import threading
lock = threading.Lock()

import networkx
import scipy.stats

import pcd.util




#scipy.stats.powerlaw
#def powlaw(a):
#    return scipy.stats.powerlaw(a).rvs()
#def A(gamma, beta, N, avg_deg, mu):
#    g = networkx.Graph()
#    for n in range(N):
#        g.add_node(n)
#        g.node[n]['k'] = powlaw(gamma)
#
#    print g.summary()
#    print sum(d['k'] for n,d in g.nodes(1))


progdir = '/home/richard/research/cd/fortunato_benchmarks'


def read_file(fname):
    for line in open(fname):
        line = line.strip().split()
        yield tuple(pcd.util.leval(x) for x in line)
def makeargs(dict):
    args = [ ]
    for k,v in dict.items():
        args.extend(['-'+k, str(v)])
    return args



def benchmark(N, k, maxk, gamma, beta, mu,
              minc="", maxc=""):
    """Unweighted, Undirected, Non-overlapping benchmark graph.

    This is the graph supporting Lancichenetti, Fortunato, Radicci,
    PRE 78 046110 (2008).

    Arguments are:

    N     # number of nodes
    k     # average degree
    maxk  # maximum degree
    gamma # exponent for the degree distribution
    beta  # exponent for the community size distribution
    mu    # mixing parameter
    minc  # minimum for the community sizes (optional)
    maxc  # maximum for the community sizes (optional)

    Example parameters:
    N=1000, k=15, maxk=100, gamma=2, beta=1, mu=.1
    """

    if minc is None: minc = ""
    if maxc is None: maxc = ""

    params = textwrap.dedent("""\
    %(N)s
    %(k)s
    %(maxk)s
    %(gamma)s
    %(beta)s
    %(mu)s
    %(minc)s
    %(maxc)s
    """%locals())

    tmpdir = tempfile.mkdtemp(dir=".", prefix="benchmark")
    open(join(tmpdir, 'parameters.dat'), 'w').write(params)

    prog = join(progdir, 'benchmark_2_2/benchmark')
    kwargs = { }
    args = [ prog ] + makeargs(kwargs)
    print "Arguments are: ", " ".join(args)

    retcode = subprocess.call(args, cwd=tmpdir)
    assert retcode == 0

    g = networkx.Graph()
    for n, c in read_file(join(tmpdir, 'community.dat')):
        g.add_node(n-1, cmty=c-1)
    for n1, n2 in read_file(join(tmpdir, 'network.dat')):
        g.add_edge(n1-1, n2-1)
    g.graph['statistics'] = open(join(tmpdir, 'statistics.dat')).read()


    shutil.rmtree(tmpdir)
    return g


def binary_graph(**kwargs):
    """
    -N              [number of nodes]
    -k              [average degree]
    -maxk           [maximum degree]
    -mu             [mixing parameter]
    -t1             [minus exponent for the degree sequence]
    -t2             [minus exponent for the community size distribution]
    -minc           [minimum for the community sizes]
    -maxc           [maximum for the community sizes]
    -on             [number of overlapping nodes]
    -om             [number of memberships of the overlapping nodes]

    -N, -k, -maxk, -mu have to be specified. For the others, the
    program can use default values:
    t1=2, t2=1, on=0, om=0, minc and maxc will be chosen close to the
    degree sequence extremes.
    If you set a parameter twice, the latter one will be taken.
    """
    prog = join(progdir, '/binary_networks/benchmark')
    args = [ prog ] + makeargs(kwargs)
    print "Arguments are: ", " ".join(args)

    retcode = subprocess.call(args)
    assert retcode == 0

    g = networkx.Graph()
    for n, c in read_file('community.dat'):
        g.add_node(n-1, cmty=c-1)
    for n1, n2 in read_file('network.dat'):
        g.add_edge(n1-1, n2-1)
    return g

def weighted_graph(**kwargs):
    """
    -N              [number of nodes]
    -k              [average degree]
    -maxk           [maximum degree]
    -mut            [mixing parameter for the topology]
    -muw            [mixing parameter for the weights]
    -beta           [exponent for the weight distribution]
    -t1             [minus exponent for the degree sequence]
    -t2             [minus exponent for the community size distribution]
    -minc           [minimum for the community sizes]
    -maxc           [maximum for the community sizes]
    -on             [number of overlapping nodes]
    -om             [number of memberships of the overlapping nodes]

    -N, -k, -maxk, -muw have to be specified. For the others, the program can use default values:
    t1=2, t2=1, on=0, om=0, beta=1.5, mut=muw, minc and maxc will be chosen close to the degree sequence extremes.
    If you set a parameter twice, the latter one will be taken.

    To have a random network use:
    -rand
    Using this option will set muw=0, mut=0, and minc=maxc=N, i.e.
    there will be one only community.
    Use option -sup (-inf) if you want to produce a benchmark whose
    distribution of the ratio of external degree/total degree is
    superiorly (inferiorly) bounded by the mixing parameter.

    Example1:
    ./benchmark -N 1000 -k 15 -maxk 50 -muw 0.1 -minc 20 -maxc 50
    Example2:
    ./benchmark -f flags.dat -t1 3
    """
    prog = join(progdir, 'weighted_networks/benchmark')
    args = [ prog ] + makeargs(kwargs)
    print "Arguments are: ", " ".join(args)

    retcode = subprocess.call(args)
    assert retcode == 0

    g = networkx.Graph()
    for x in read_file('community.dat'):
        n, cmtys = x[0], x[1:]
        cmtys = [c-1 for c in cmtys]
        g.add_node(n-1, cmtys=cmtys)
    for n1, n2, weight in read_file('network.dat'):
        g.add_edge(n1-1, n2-1, weight=weight)
    g.graph['statistics'] = open('statistics.dat').read()
    return g


def hierarchical_graph(**kwargs):
    """
    -N              [number of nodes]
    -k              [average degree]
    -maxk           [maximum degree]
    -t1             [minus exponent for the degree sequence]
    -t2             [minus exponent for the community size distribution]
    -minc           [minimum for the micro community sizes]
    -maxc           [maximum for the micro community sizes]
    -on             [number of overlapping nodes]
    -om             [number of memberships of the overlapping nodes]
    -minC           [minimum for the macro community size]
    -maxC           [maximum for the macro community size]
    -mu1            [mixing parameter for the macro communities (see Readme file)]
    -mu2            [mixing parameter for the micro communities (see Readme file)]

    Example2:
    ./hbenchmark -f flags.dat
    ./hbenchmark -N 10000 -k 20 -maxk 50 -mu2 0.3 -minc 20 -maxc 50 -minC 100 -maxC 1000 -mu1 0.1
    """
    prog = join(progdir, 'hierarchical_bench2_2/hbenchmark')
    args = [ prog ] + makeargs(kwargs)
    print "Arguments are: ", " ".join(args)

    retcode = subprocess.call(args)
    assert retcode == 0

    g = networkx.Graph()
    for n, c in read_file('community_first_level.dat'):
        g.add_node(n-1, microC=c-1)
    for n, c in read_file('community_second_level.dat'):
        g.add_node(n-1, macroC=c-1)
    for n1, n2 in read_file('network.dat'):
        g.add_edge(n1-1, n2-1)
    return g


def setCmtyAssignmentsOverlap(G, g):
    for node, data in g.node.iteritems():
        for c in data['cmtys']:
            G.cmtyListAdd(c, G._nodeIndex[node])
    return G

def setCmtyAssignments(G, g, keyname='cmty'):
    for node, data in g.node.iteritems():
        c = data[keyname]
        G.cmtyListAdd(c, G._nodeIndex[node])
    return G


if __name__ == "__main__":
    #A(1, 1, 100, 10, 1)

    g = hierarchical_graph(
        N=10000, k=10, maxk=100,
        mu1=.3, mu2=.4,
        minC=500, maxC=2000,
        minc=10, maxc=100,
        )
    print networkx.info(g)

    g = binary_graph(
        N=1000, k=15, maxk=50, mu=.1,
        )
    print networkx.info(g)



    g = benchmark()
    networkx.write_gml(g, "graph.gml")
