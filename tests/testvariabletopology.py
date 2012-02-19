# Richard Darst, February 2011

import pcd.graphs

basename = 'tests-output/vartop/'

G = pcd.graphs.polopa_tribes_G()

MR = pcd.MultiResolution()
MR.run([G]*4, gammas=dict(low=.003, high=5, density=30))
MR.plot(fname=basename+"tribes_mr.png")


print G.imatrix
G.imatrix[:] -= G.imatrix.max()
print G.imatrix

G.make_sparse(default=1, cutoff_op=lambda x,y:x!=y)
G.enableVT()
print G.rmatrix
print G.srmatrix


G.trials(gamma=.1, trials=10)
G2 = G.copy()
#G2.make_sparse('auto')
#from fitz import interactnow


MR = pcd.MultiResolution(
    #overlap=False,
    #savefigargs=dict(fname="tests-outputs/vartop/tribes_gamma%(gamma)09.4f"),
    #plotargs=dict(fname=basename+"tribes_mr.png"),
    )
MR.run([G]*4, gammas=dict(low=.3, high=5, density=30))
MR.plot(fname=basename+"tribes_VT_mr.png")


##G.viz('out.png')
#G.trials(gamma=1, trials=10)
#
#g = G._graph #make_networkx()
#gnew = g
#import pcd.layout
##gnew = pcd.layout.layoutgraph(g)
#G.colors_to_networkx(gnew)
#pcd.layout.invertWeights(gnew)
#pcd.layout.delCoords(gnew)
##from fitz import interactnow
#
#gnew = pcd.layout.layoutgraph(gnew)
#open('tmp.dot','w').write(pcd.layout.toDotStr(gnew))
#
#pcd.layout.writePng(gnew, basename+'tribes_out.png')
#
#pcd.layout.getCoordMap(G, gnew)



G = pcd.graphs.dolphins_G()
MR = pcd.MultiResolution()
MR.run([G]*4, gammas=dict(low=.01, high=10, density=30))
MR.plot(fname=basename+"dolphins_mr.png")


g = pcd.graphs.dolphins(weightFriend=-1)
G = pcd.Graph.fromNetworkX(g, defaultweight=0)

G.make_sparse(default='auto')
G.enableVT()

MR = pcd.MultiResolution(
    #savefigargs=dict(fname="tests-outputs/vartop/tribes_gamma%(gamma)09.4f"),
    #plotargs=dict(fname=basename+"tribes_mr.png"),
    )
MR.run([G]*4, gammas=dict(low=.01, high=1, density=30))
MR.plot(fname=basename+"dolphins_VT_mr.png")
