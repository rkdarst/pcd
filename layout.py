# Richard Darst, February 2012

import subprocess
import cStringIO as StringIO

import networkx
import pydot

def toDotStr(g):
    return networkx.to_pydot(g).to_string()
def fromDotStr(s):
    return networkx.from_pydot(pydot.graph_from_dot_data(s))

def coordsFromFile(fname):
    #coords = networkx.read_dot('/home/richard/dolphins-coords-messy.gv')
    #graph = pcd.graphs.dolphins(weightFriend=-1)
    pass

def layoutgraph(g):
    """Given a Graph object, """

    #g = G.make_networkx()
    #G.colors_to_networkx(g)
    #networkx.write_gml(g.copy(), dirname+'dolphins_gamma%s.gml'%gamma)

    pd = networkx.to_pydot(g)
    #pd.to_string()

    prog = 'neato'
    #prog = 'fdp'
    p = subprocess.Popen('%s -Goverlap=scalexy -Gsplines=true -Gstart=5'%prog,
                         stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                         shell=True)
    out = p.communicate(pd.to_string())
    #from fitz import interactnow

    gnew = networkx.from_pydot(pydot.graph_from_dot_data(out[0]))
    return gnew

def writePng(g, fname):
    p = subprocess.Popen(['neato','-Tpng',
                          '-Goverlap=scalexy',
                          '-Gsplines=true',
                          '-Gstart=5',
                          '-Epenwidth=2',
                          '-Npenwidth=2',
                          '-Nfontsize=20',
                          '-o%s'%fname],
              stdin=subprocess.PIPE)
    p.stdin.write(toDotStr(g))

def invertWeights(g, func=lambda x: -x):
    for a,b,dat in g.edges(data=True):
        if 'weight' in dat:
            dat['weight'] = func(float(dat['weight'].strip('"\'')))
def invertWeights(g, func=lambda x: -x):
    for a,b,dat in g.edges(data=True):
        if 'weight' in dat:
            try:
                weight = float(dat['weight'])
            except ValueError:
                weight = float(dat['weight'].strip('"\''))
            print weight,
            if weight < 0:
                newweight = 1
            else:
                newweight = 0
            print newweight
            dat['weight'] = newweight

def delCoords(g):
    for a,dat in g.nodes(data=True):
        dat.pop('pos', None)
    for a,b,dat in g.edges(data=True):
        dat.pop('pos', None)


def getCoordMap(G, g):
    """Takes coords from g"""
    coords = { }
    for n in range(G.N):
        c = g.node[G._nodeLabel[n]]['pos']
        c = tuple(float(x) for x in c.strip('"\'').split(','))
        coords[n] = c
