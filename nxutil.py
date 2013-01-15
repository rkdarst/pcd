# Richard Darst, July 2012

"""Module to handle networkx graphs and their interfaces.

Networkx graphs have two ways of representing the community:

- node attribute 'cmty' which indicates the respective community of
  the graph.
- node attribute 'cmtys' which is an iterable listing muultiple
  communities for the node.
"""

import collections
import numpy

def _iterCmtys(d):
    """Iterate over the communities in a node' data dict.

    This abstracts out whether or not there is one community or
    multiple communities for the node."""
    if 'cmty' in d:
        #raise
        if 'cmtys' in d:
            raise ValueError('Multiple "cmty" and "cmtys" both exist.')
        yield d['cmty']
    else:
        for c in d['cmtys']:
            yield c

def setCmtyAssignments(G, g):
    """Loads communities from a networkx graph g into a pcd graph G.

    Make sure you clear the G  first.  G.cmtyListClear()"""
    G.oneToOne = 0
    for node, data in g.nodes_iter(data=True):
        #print node, data['cmtys']
        for c in _iterCmtys(data):
            G.cmtyListAddOverlap(c, G._nodeIndex[node])
    return G

def cmtyInit(g):
    """Clear/initialize all communities in a networkx graph."""
    for n, d in g.nodes_iter(data=True):
        d.pop('cmty', 0)
        d['cmtys'] = set()
def cmtyClear(g):
    for n, d in g.nodes_iter(data=True):
        d.pop('cmty', 0)
        d.pop('cmtys', 0)
def cmtyAddFromList(g, cmtyID, nodes):
    """For a networkx graph, add all 'nodes' to 'cmtyID'."""
    print cmtyID, nodes
    #for n in g.nbunch_iter(nodes):
    for n in nodes:
        g.node[n]['cmtys'].add(cmtyID)


def graphproperties(g):
    import collections
    cmty_d = collections.defaultdict(set)
    for node, data in g.nodes_iter(data=True):
        for c in _iterCmtys(data):
            cmty_d[c].add(node)
    print "number of communities:", len(cmty_d)
    print "mean cmty size:", numpy.mean([len(v) for v in cmty_d.values()])
    print "std cmty size:", numpy.std([len(v) for v in cmty_d.values()])
    U = set(g.nodes_iter())
    U_c = set()
    #print cmty_d
    #print [ v for v in cmty_d.itervalues() ]
    [ U_c.update(v) for v in cmty_d.itervalues() ]
    print "total nodes:", len(U)
    print "total nodes in communities:", len(U_c)


def communities(g):
    raise Exception
    """Given a graph, return a mapping c -> nodes_in_c."""
    cmtys = collections.defaultdict(set)
    for node, data in g.nodes_iter(data=True):
        for c in _iterCmtys(data):
            cmtys[c].add(node)
    return cmtys
cmtynodes = communities
def nodecmtys(g):
    raise Exception
    """Given a graph, return mapping n -> its communities"""
    cmtys = collections.defaultdict(set)
    for node, data in g.nodes_iter(data=True):
        for c in _iterCmtys(data):
            cmtys[node].add(c)
    return cmtys

def ed_defined(g, nodecmtys):
    raise Exception
    cmtynodes = collections.defaultdict(set)
    for n, cs in nodecmtys.iteritems():
        for c in cs:
            cmtynodes[c].add(n)
    #cmtysizes = dict((c, len(ns)) for c,ns in communities(g).iteritems())
    cmtysizes = dict((c, len(ns)) for c,ns in cmtynodes.iteritems())
    well_defined = 0
    illnodes = set()
    for n0, n0data in g.nodes_iter(data=True):
        #c0 = g.node[n0]['cmty']
        print "zzz", n0, n0data
        #for c0 in _iterCmtys(n0data):
        for c0 in nodecmtys[n0]:
            print "zzz cmty", c0
            degree = 0
            degree_int = 0
            degree_ext_by_cmty = collections.defaultdict(int)
            #from fitz import interact ; interact.interact('b> ')
            for n1 in g.neighbors_iter(n0):
                print "xx", n0, n1
                n1data = g.node[n1]
                #c1 = g.node[n1]['cmty']
                #n1cmtys = set(_iterCmtys(n1data))
                n1cmtys = set(nodecmtys[n1])
                print n0, n0data, n1, n1data
                degree += 1
                #if c0 == c1:
                #    continue
                #degree_ext_by_cmty[c1] += 1
                for c1 in n1cmtys:
                    if c0 == c1:
                        degree_int += 1
                    else:
                        degree_ext_by_cmty[c1] += 1
            int_density = float(degree_int)/(cmtysizes[c0]-1)
            print 'y', int_density, degree_ext_by_cmty
            if len(degree_ext_by_cmty) == 0:
                max_ext_density = 0
            else:
                max_ext_density = \
                            max((float(d) / cmtysizes[c])
                                for (c,d) in degree_ext_by_cmty.iteritems())
        if int_density > max_ext_density:
            g.node[n0]['ill'] = False
            well_defined += 1
        else:
            g.node[n0]['ill'] = True
            illnodes.add(n0)
    return dict(
        illnodes=illnodes,
        fracwell=float(well_defined)/len(g),
        )

def cmty_mapping(g0, g1):
    raise Exception
    """For each community in g0, find the community in g1 which most
    overlaps it.  Return a mapping g0 communities -> g1 communities.

    cmty_mapping[c0] = the c1 which contains the most intersected nodes

    Or: every community it g0 -> every best commmunity in g1.
    """
    cmtys0 = communities(g0)
    cmtys1 = communities (g1)
    cmty0_to_cmty1 = dict()
    for c0 in cmtys0:
        maxoverlap = 0
        bestcmty = None
        c0nodes = cmtys0[c0]
        for c1 in cmtys1:
            c1nodes = cmtys1[c1]
            overlap = float(len(c0nodes & c1nodes))/len(c0nodes)
            if overlap > maxoverlap:
                maxoverlap = overlap
                bestcmty = c1
        cmty0_to_cmty1[c0] = c1
    return cmty0_to_cmty1

def frac_detected(gPlanted, gDetected):
    raise Exception
    n_detected = 0
    n_well_detected = 0
    n_well = sum(1 for (n,d) in gPlanted.nodes_iter(data=True) if not d['ill'])
    bestDetected = cmtyMapping(gPlanted, gDetected)
    cmtysPlanted = communities(gPlanted)
    cmtysDetected = communities(gDetected)
    # for each planted communitiy
    # select the best detected community
    for cPlanted, cDetected in baseDetected:
        # how many nodes were correctly detected?  Add to correct number.
        detected = cmtysPlanted[cPlanted] & cmtysDetected[cDetected]
        n_detected += len(detected)
        n_well_detected += len(set(1 for n in cPlanted
                                  if not gPlanted.node[n]['ill'])
                              &
                              cDetected)

    # How many nodes are well defined (count them).  How many are well
    # defined and correctly detected?  (increment).

    # Return values
    return (
        float(n_detected) / len(gPlanted),
        float(n_well_detected) / n_well,
        )
