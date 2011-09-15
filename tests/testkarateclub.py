# Richard Darst, September 2011

import pcd
import graphs

fast = globals().get('fast', False)

g = graphs.karate_club()
for a,b,d in g.edges(data=True):
    print a, b, d
G = pcd.Graph.fromNetworkX(g, defaultweight=1)
G.minimize(gamma=.1)
if not fast:
    G.viz()

lnumber = 50
if fast:
    lnumber = 10

MR = pcd.MultiResolution(.05, 50, number=lnumber)
MR.do([G]*10, trials=10, threads=2)
MR.write('tmp-karate.txt')
MR.plot("karate_.png")
