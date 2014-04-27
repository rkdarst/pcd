
import os
import random
import subprocess

import networkx

import pcd.util

#binary = './gb_code/programmi_last/a.out'
binary = '/home/darstr1/proj/growsf/gb_code/programmi_last/a.out'
# Usage: ./a.out T (maximal time) p beta kappa

with_seed = True

def growsf(N, p, beta, kappa):
    """Gnerate a growing graph with neighbor attachment.

    N: int
        number of nodes
    p: float
        1-(neighbor attachment probability)
    beta:
        fitness parameter.  0=all fitness same
    kappa:
        another fitness parameter.

    The fitness is calculated via:
        fitness = exp(-beta * x**(1./(kappa+1.)) )
        with x being a uniform random variable [0,1)
    """
    seed = random.randint(0, 2**32-1)
    with pcd.util.tmpdir_context(prefix='tmp-growsfgb-', chdir=True):
        #print os.getcwd()
        args = [ binary ]
        args.extend(("%d %f %f %f"%(N, p, beta, kappa)).split())
        if with_seed:
            args.append("%d"%seed)

        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        assert proc.wait() == 0

        g = networkx.read_edgelist('clust_net.txt', nodetype=int)
        g.graph['N'] = N
        g.graph['p'] = float(p)
        g.graph['beta'] = float(beta)
        g.graph['kappa'] = float(kappa)
        return g

# static_original.c
def marsili(N, T, lambda_, xi):
    """Marsili model

    T: int
        Number of equilibration steps.
    lambda_: float
        Lambda value (frequency of deleting nodes)
    xi: float
        Probably of links via friends.
    """
    seed = random.randint(0, 2**32-1)
    binary = '/home/darstr1/proj/growsf/gb_code/marsili/static_original'
    with pcd.util.tmpdir_context(prefix='tmp-marsili-', chdir=True):
        #print os.getcwd()
        args = [ binary ]
        args.extend(("%d %d %f %f"%(N, T, lambda_, xi)).split())
        if with_seed:
            args.append("%d"%seed)

        proc = subprocess.Popen(args)#, stdout=subprocess.PIPE)
        assert proc.wait() == 0

        g = networkx.read_edgelist('data.edges', nodetype=int,
                                   create_using=networkx.Graph())
        g.graph['T'] = T
        g.graph['lambda'] = float(lambda_)
        g.graph['xi'] = float(xi)
        if with_seed:
            g.graph['seed'] = seed
    # Remove self-loops?
    for n in g.nodes_with_selfloops():
        g.remove_edge(n, n)
    return g




# duplication.c
def sole(T, delta, alpha):
    """Sole model

    T: int
        Size of network
    delta: float
        Probability of link deletion.  z_old*(1-delta) is average degree
        for a new node.  Misnamed nu in the help.
    alpha: float
        Probability of linking old and new node.
    """
    seed = random.randint(0, 2**32-1)
    binary = '/home/darstr1/proj/growsf/gb_code/sole/duplication'
    with pcd.util.tmpdir_context(prefix='tmp-sole-', chdir=True):
        #print os.getcwd()
        args = [ binary ]
        args.extend(("%d %f %f"%(T, delta, alpha)).split())
        if with_seed:
            args.append("%d"%seed)

        proc = subprocess.Popen(args)#, stdout=subprocess.PIPE)
        assert proc.wait() == 0

        g = networkx.read_edgelist('data_duplication.edges', nodetype=int,
                                   create_using=networkx.Graph())
        g.graph['T'] = T
        g.graph['delta'] = float(delta)
        g.graph['alpha'] = float(alpha)
        if with_seed:
            g.graph['seed'] = seed
        return g


