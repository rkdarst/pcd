# Richard Darst, September 2011

import models
import graphs

g = graphs.karate_club()
for a,b,d in g.edges(data=True):
    print a, b, d
G = models.Graph.fromNetworkX(g, defaultweight=1)
G.minimize(gamma=.1)
G.viz()

MR = models.MultiResolution(.05, 50, number=50)
MR.do([G]*10, trials=10, threads=2)
MR.write('tmp-karate.txt')
MR.plot("karate_.png")
