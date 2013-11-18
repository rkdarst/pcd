
import os

import pcd.util
import subprocess

import networkx

#binary = './gb_code/programmi_last/a.out'
binary = '/home/darstr1/proj/growsf/gb_code/programmi_last/a.out'
# Usage: ./a.out T (maximal time) p beta kappa

def growsf(N, p, beta, kappa):
    with pcd.util.tmpdir_context(prefix='tmp-growsfgb-', chdir=True):
        #print os.getcwd()
        args = [ binary ]
        args.extend(("%d %f %f %f"%(N, p, beta, kappa)).split())

        p = subprocess.Popen(args)
        p.wait()

        g = networkx.read_edgelist('clust_net.txt', nodetype=int)
        return g

