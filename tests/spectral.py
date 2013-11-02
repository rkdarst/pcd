
# compare with newman2013spectral

import pcd.support.spectral
import pcd.sbmdegree
import pcd.support.powerlaw


for i in range(1):
    g = pcd.sbmdegree.graph(2, 500, pcd.support.powerlaw.PowerLaw(power=-2.5, xmin=10,xmax=1e6),
                            lambda_=.999)
    sf = pcd.support.spectral.Flow(g, num_ev=5)
    sb = pcd.support.spectral.NonBacktrackingWalk(g, num_ev=5)

    # sf.cmty_classify(-1).values()[-10:]
    asf = sf.cmty_classify(-2)
    asb = sb.cmty_classify(-2)
    #pylab.plot(a.values())
    #pyplot.ylim(-.5, .5)

    print sf.q

from fitz import interactnow
