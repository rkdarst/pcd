
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

