# Richard Darst, December 2012

import collections
import numpy
import textwrap

import pcd.nxutil

class Communities(object):
    """Store and operate on a community assignment.

    cmtynodes: a dictionary mapping community_id -> set(cmty_nodes)
    nodes: a set of all possible nodes in the system
    name: if given, a string to be stored at self.name to use to identify
    """
    def __init__(self, cmtynodes, nodes, name=None):
        self._cmtynodes = cmtynodes
        self.nodes = nodes
        if name:
            self.name = name

    @classmethod
    def from_pcd(cls, G):
        """Convert a pcd.Graph into Communities object.

        Communities are loaded from pcd.Graph, using G._nodeLabel."""
        d = G.cmtyDict()
        cmtynodes = { }
        if hasattr(G, '_nodeLabel'):
            for c, nodes in d.iteritems():
                cmtynodes[c] = set(G._nodeLabel[n] for n in nodes)
            nodes = set(G._nodeIndex.iterkeys())
            assert len(G._nodeIndex) == len(G._nodeLabel)
        else:
            for c, nodes in d.iteritems():
                cmtynodes[c] = set(nodes)
            nodes = set(range(G.N))
        return cls(cmtynodes=cmtynodes, nodes=nodes)

    @classmethod
    def from_networkx(cls, g):
        """Create a new Communities object from a networkx.Graph object.

        The communities are stored in each node data dict (g.node[n])
        in one of these keys (the first takes priority):
          g.node[n]['cmty'] = c
          g.node[n]['cmtys'] = set((c1,c2,c3,...))
        """
        nodes = set(g.nodes_iter())
        cmtynodes = { }
        for node, d in g.nodes_iter(data=True):
            for c in pcd.nxutil._iterCmtys(d):
                if c not in cmtynodes: cmtynodes[c] = set()
                cmtynodes[c].add(node)
        return cls(cmtynodes=cmtynodes, nodes=nodes)

    # I/O

    def to_pcd(self):
        """Load this community structure into a pcd.Graph structure.

        Returns pcd.Graph with this community structure stored within
        it."""
        G = pcd.Graph(N=len(self.nodes), sparse=True)
        G._makeNodeMap(self.nodes)
        G = self.load_pcd(G, clear=True)
        return G
    def to_pcd_cmtystate(self):
        """Return the community structure as pcd.Graph state.

        This is suitable for usage in pcd.Graph.setcmtystate().  Thin
        wrapper around self.to_pcd().getcmtystate() (but could someday
        be optimized beyond that)."""
        G = self.to_pcd()
        return G.getcmtystate()

    def load_pcd(self, G, clear=True):
        """Loads communities from a networkx graph g into a pcd graph G.

        clear: if True, zero all present communities before loading.
        Otherwise, simply add all nodes to the new communities."""
        cmtynodes = self.cmtynodes()
        non_overlapping = self.is_non_overlapping()
        assert max(cmtynodes.iterkeys()) < G.N
        assert all(isinstance(c, int) for c in cmtynodes.iterkeys())
        if clear:
            G.cmtyListClear()

        G.oneToOne = 0
        if non_overlapping and self.is_cover():
            G.oneToOne = 1
        for c, nodes in cmtynodes.iteritems():
            for n in nodes:
                if not clear and G.cmtyContains(c, G._nodeIndex[n]):
                    # Ignore nodes which are already in the community
                    # (raises an error in pcd code)
                    continue
                G.cmtyListAddOverlap(c, G._nodeIndex[n])
                if non_overlapping:
                    G.cmty[G._nodeIndex[n]] = c
        return G

    # Status information
    @property
    def N(self):
        """Number of nodes"""
        return len(self.nodes)
    @property
    def q(self):
        """Number of communities"""
        return len(self._cmtynodes)

    def is_cover(self, nodes=None, strict=True):
        """Is every node in at least one community?

        nodes: if given, use this as the universe of nodes to be covered.
        strict: if true, raise an except if we span more than set of nodes.
        """
        spannednodes = set()
        if nodes is None:
            nodes = self.nodes
        for ns in self._cmtynodes.itervalues():
            spannednodes |= ns
        if strict:
            if len(spannednodes - nodes) > 0: # No spanned nodes left over.
                raise ValueError("We cover more than the universe of nodes, and strict mode is enabled.")
        return len(nodes - spannednodes) == 0 # Every node in spanned nodes.
    def is_partition(self, nodes=None):
        """Every node in exactly one community.

        nodes: if given, use this as the universe of nodes to be
        covered.

        Should be equivalent to is_non_overlapping and is_cover."""
        # This should be identical to is_cover and non_overlapping
        if nodes is None:
            nodes = self.nodes
        spannednodes = set()
        total_nodes = 0
        for ns in self._cmtynodes.itervalues():
            spannednodes |= ns
            total_nodes += len(ns)
        return (    total_nodes == len(nodes)   # no node in multiple cmtys
                and len(nodes - spannednodes) == 0 ) # every node covered.

    def is_non_overlapping(self):
        """No node is in more than one community."""
        total_nodes = 0
        spannednodes = set()
        for ns in self._cmtynodes.itervalues():
            spannednodes |= ns
            total_nodes += len(ns)
        return total_nodes == len(spannednodes)

    def stats(self):
        cmtynodes = self.cmtynodes()
        cmtynodes = dict((k,v) for (k,v) in cmtynodes.iteritems() if v)

        stats = [ ]
        stats.append("number of nodes: %d"%self.N)
        if hasattr(self, 'g'):
            stats.append("number of edges: %d"%self.g.number_of_edges())
        stats.append("number of communities: %d"%len(cmtynodes))
        stats.append("mean community size: %f"%numpy.mean(
            [len(c) for c in cmtynodes.values()]))
        stats.append("std community size: %f"%numpy.std(
            [len(c) for c in cmtynodes.values()]))
        cmtys_by_rev_size = [c for (c,ns) in sorted(cmtynodes.iteritems(),
                                                    reverse=True,
                                                    key=lambda (c,ns):len(ns))]
        stats.append("community sizes: %s"%' '.join(
            str(len(cmtynodes[c])) for c in cmtys_by_rev_size))
        stats.append("fraction of nodes in largest community: %f"%(
            max(len(c) for c in cmtynodes.values())/float(self.N)))
        return stats
    def list_communities(self, width=120):
        stats = [ ]
        stats.append("List of communities:")
        singleton_cmtys = [ ]
        cmtynodes = self.cmtynodes()
        cmtys_by_rev_size = [c for (c,ns) in sorted(cmtynodes.iteritems(),
                                                    reverse=True,
                                                    key=lambda (c,ns):len(ns))]
        for c in cmtys_by_rev_size:
            cmtysize = len(cmtynodes[c])
            if cmtysize == 1:
                singleton_cmtys.append(c)
                continue
            stats.append(textwrap.fill(
                (("%3s: "%c)+' '.join(str(n) for n in sorted(cmtynodes[c]))),
                width=width,
                initial_indent="", subsequent_indent="     ",
                ))
        stats.append(textwrap.fill(
            "Nodes in singleton communities: "+' '.join(
                str(cmtynodes[c].pop()) for c in singleton_cmtys),
            width=width,
            initial_indent="", subsequent_indent="     ",
            ))
        return stats
    def list_overlapping_nodes(self, width=120):
        stats = [ ]
        stats.append("How many communities does each node belong to?:")
        cmty_counts = collections.defaultdict(set)
        nodecmtys = self.nodecmtys()
        for node, cmtys in nodecmtys.iteritems():
            cmty_counts[len(cmtys)].add(node)
        for num_cmtys, nodes in cmty_counts.iteritems():
            stats.append(textwrap.fill(
                (("%4s: "%num_cmtys)+' '.join(str(n) for n in sorted(nodes))),
                width=width,
                initial_indent="", subsequent_indent="      ",
                ))
        return stats

    def nodeInfo(self, n, g):
        """Print detailed information about a node, including what
        communities it is connected to.

        Uses edges (not neighbors).
        """
        adjacency = [ ]
        adjacency.extend(neighbor for n,neighbor in g.edges(n))
        nodecmtys = self.nodecmtys()

        print "Node %s in communities %s"%(n, nodecmtys[n])
        #print "Edges connect to %s"%sorted(adjacency)
        cmtys = collections.defaultdict(list)
        #print "Connected communities for each adjacent node:"
        for n2 in sorted(adjacency):
            #print "  %10s: cmtys %s"%(n2, sorted(nodecmtys[n2]))
            for c in nodecmtys[n2]:
                cmtys[c].append(n2)
        print "Connections to each community:"
        for c, nodes in sorted(cmtys.iteritems()):
            print " %3s: %3d edges, %s"%(c, len(nodes), sorted(nodes))



    # Transformations

    def cmtynames(self):
        return set(self._cmtynodes.keys())
    def cmtynodes(self, copy=False):
        """Return mapping cmty -> node set.

        Returns a dictionary {c0:set(n00,n01,...), c1:set(n10,n11), ...}

        Warning: do not mutate the returned dictionary unless copy=True."""
        if not copy:
            return self._cmtynodes
        else:
            cmtynodes = { }
            for c, nodes in self._cmtynodes.iteritems():
                cmtynodes[c] = set(nodes)
            return cmtynodes
    def cmtysizes(self):
        return dict((c, len(ns)) for c,ns in self._cmtynodes.iteritems())

    def nodecmtys(self, copy=False):
        """Return mapping node -> cmty set.

        Returns a dictionary {n0:set(c00,c01,...), n1:set(c10,c11), ...}

        Warning: do not mutate the returned dictionary unless copy=True.
        """
        nodecmtys = { }
        for c, nodes in self._cmtynodes.iteritems():
            for n in nodes:
                nodecmtys.setdefault(n, set())
                nodecmtys[n].add(c)
        return nodecmtys
    def nodecmtys_onetoone(self):
        """Return mapping node -> community if there are no overlaps.

        If there are overlaps, raise an exception."""
        nodecmtys = { }
        for c, nodes in self._cmtynodes.iteritems():
            for n in nodes:
                if n in nodecmtys:
                    raise ValueError("Overlapping: node %s in cmtys %s and %s"%
                                     (n, nodecmtys[n], c))
                nodecmtys[n] = c
        return nodecmtys

    def cmty_mapping(self, otherCmtys, mode="overlap"):
        """For each community in self, find the community in g1 which most
        overlaps it.  Return a mapping g0 communities -> g1 communities.

        g0 is planted, g1 is detected.

        cmty_mapping[c0] = the c1 which contains the most intersected nodes

        Or: every community it g0 -> every best commmunity in g1.

        Comparison modes:
        overlap: which c1 has the most overlapping nodes with c0
        F1: which c0/c1 has the greatest F-score overlap.
        """
        cmtys0 = self.cmtynodes()
        cmtys0list = list(cmtys0.iterkeys())
        cmtys1 = otherCmtys.cmtynodes()
        cmtys1list = list(cmtys1.iterkeys())
        cmty0_to_cmty1 = dict()
        assigned_cmtys = set()


        overlaps = numpy.zeros(shape=(len(cmtys0), len(cmtys1)), dtype=int)
        #overlaps[detected, planted]

        for i0, c0 in enumerate(cmtys0list):
            bestScore = 0
            bestcmty = None
            c0nodes = cmtys0[c0]
            for i1, c1 in enumerate(cmtys1list):
                c1nodes = cmtys1[c1]
                overlap = len(c0nodes & c1nodes)
                overlaps[i0,i1] = overlap
                score = 0
                if mode == "overlap":
                    # Using overlap / planted community size (note:
                    # constant divisor, so the division doesn't do
                    # anything really.)
                    if c1 in assigned_cmtys:
                        continue
                    if len(c0nodes) > 0:
                        score = float(overlap)/len(c0nodes)
                elif mode == "F1":
                    # Using f-score
                    precision = float(overlap)/len(c1nodes)
                    recall    = float(overlap)/len(c0nodes)
                    F1 = 2. * precision * recall / (precision + recall)
                    score = F1
                elif mode == "newman":
                    pass
                    #overlaps[c0,c1] = overlap
                else:
                    raise ValueError("Unknown mode: %s"%mode)

                if bestScore is None or score > bestScore:
                    bestScore = score
                    bestcmty = c1
            cmty0_to_cmty1[c0] = c1
            assigned_cmtys.add(c1)
        if mode == "xxx":
            pass
        if mode == "newman":
            cmty0_to_cmty1 = { }
            # Newman algorithm in newman2004fast, footnote #19.
            plant_to_detect = collections.defaultdict(set)
            best_overlap_of_detected = overlaps.argmax(1)
            for i_planted, i_detected in enumerate(best_overlap_of_detected):
                c_planted = cmtys0list[i_planted]
                c_detected= cmtys1list[i_detected]
                plant_to_detect[c_detected].add(c_planted)
            for c_detect, c_planteds in plant_to_detect.iteritems():
                if len(c_planteds) == 1:
                    #print "mapping: c %s to c %s"%(c_planteds, c_detect)
                    cmty0_to_cmty1[c_planteds.pop()] = c_detect
                else:
                    #print "mapping: NOT done: %s all map to %s"%(
                    #    sorted(c_planteds), c_detect)
                    pass
            #from fitz import interactnow
        return cmty0_to_cmty1


    def frac_detected(self, detected, cmtyMapping, n_nodes, limit_nodes=None,
                      full=False):
        """Fraction of nodes properly identified.

        self: planted
        detected: detected communities, another Communities object.
        """
        if isinstance(cmtyMapping, tuple):
            cmtyMapping = self.cmty_mapping(cmtyMapping[0], **cmtyMapping[1])
        n_detected = 0
        n_total = 0
        cmtysPlanted = self.cmtynodes()
        cmtysDetected = detected.cmtynodes()
        allDetected = set()

        # for each planted communitiy select the best detected community
        for cPlanted, cDetected in cmtyMapping.iteritems():
            # How many nodes are planted?
            plantedNodes = cmtysPlanted[cPlanted]
            detectedNodes = cmtysDetected[cDetected]

            if limit_nodes is not None:
                plantedNodes = plantedNodes & limit_nodes  # explicit copy
            # how many nodes were correctly detected?  Add to correct number.
            #detected = cmtysPlanted[cPlanted] & cmtysDetected[cDetected]
            detected = plantedNodes & detectedNodes
            if limit_nodes is not None:
                detected &= limit_nodes
            # If we use the line below, we miscalculate the total
            # number of nodes we are trying to detect.
            n_total += len(plantedNodes)
            n_detected += len(detected)
            allDetected |= detected
            #print sorted(plantedNodes-detectedNodes)
            #print sorted(detectedNodes-plantedNodes)
            #print cPlanted, cDetected, len(plantedNodes), len(detected)

        # This is done because plantedNodes only includes planted
        # cmtys where there was any match.  We need to know the size of the entire universe.
        n_total = n_nodes
        if n_total == 0:
            #print limit_nodes
            return float('nan')
            #return 0
        if full:
            return float(n_detected) / n_total, allDetected
        return float(n_detected) / n_total


    def illdefined(self, g,
                   mode='internaledgedensity',
                   insist_cover=True,
                   insist_non_overlapping=True):

        nodecmtys = self.nodecmtys()
        cmtynodes = self.cmtynodes()
        cmtysizes = self.cmtysizes()

        n_well_defined = 0
        n_ill_defined = 0
        illnodes = set()
        for n0, n0data in g.nodes_iter(data=True):
            #c0 = g.node[n0]['cmty']
            #print "zzz", n0, n0data
            #for c0 in _iterCmtys(n0data):
            if insist_cover and len(nodecmtys[n0]) < 1:
                raise ValueError("Node $r is in no communities"%n0)
            if insist_non_overlapping and len(nodecmtys[n0]) > 1:
                raise ValueError("Node $r is in multiple communities"%n0)
            if len(nodecmtys[n0]) > 1:
                raise NotImplementedError("We don't handle this case yet.")
            for c0 in nodecmtys[n0]:
                #print "zzz cmty", c0
                degree = 0
                degree_int = 0
                degree_ext = 0
                degree_ext_by_cmty = collections.defaultdict(int)
                #from fitz import interact ; interact.interact('b> ')
                # WARNING: fix for DiGraph or MultiGraph
                for n1 in g.neighbors_iter(n0):
                    n1cmtys = set(nodecmtys[n1])
                    #print "xx", n0, n1, n1cmtys
                    if insist_cover and len(n1cmtys) < 1:
                        raise ValueError("Node $r is in no communities"%n1)
                    if insist_non_overlapping and len(n1cmtys) > 1:
                        raise ValueError("Node $r is in multiple communities"
                                                                           %n1)
                    degree += 1
                    for c1 in n1cmtys:
                        if c0 == c1:
                            degree_int += 1
                        else:
                            degree_ext += 1
                            degree_ext_by_cmty[c1] += 1
                int_density = float(degree_int)/(cmtysizes[c0]-1)

                density_ext_by_cmty = dict(
                      (_c, float(_deg)/cmtysizes[_c])
                      for _c,_deg in degree_ext_by_cmty.iteritems()
                      if _deg != 0
                    )
                if len(density_ext_by_cmty) == 0:
                    max_ext_density = 0
                else:
                    max_ext_density = max(density_ext_by_cmty.itervalues())


                #print 'y', int_density, max_ext_density

            # At what level of indention should this be?
            if mode == "degree":
                if degree_int > max(degree_ext_by_cmty.itervalues()):
                    n_well_defined += 1
                else:
                    n_ill_defined +=1
                    illnodes.add(n0)
            elif mode == "internaledgedensity":
                #if int_density > max_ext_density:
                if int_density >= max_ext_density:
                    #g.node[n0]['ill'] = False
                    n_well_defined += 1
                else:
                    #g.node[n0]['ill'] = True
                    n_ill_defined +=1
                    illnodes.add(n0)
            else:
                raise ValueError("Unrecognized mode: %s"%mode)
        assert n_well_defined + n_ill_defined == len(g)
        assert len(self.nodes) == len(g)
        assert n_ill_defined == len(illnodes)


        return dict(
            illnodes=illnodes,
            fracwell=float(n_well_defined)/len(g),
            )
    def illnodes(self, g, *args, **kwargs):
        return self.illdefined(g, *args, **kwargs)['illnodes']



    def VI(self, other):
        import pcd.util
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return pcd.util.VI(G1, G2)
    def ovIn(self, other):
        import pcd.util
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return pcd.util.N(G1, G2)
    def ovInw(self, other):
        import pcd.util
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return pcd.util.N(G1, G2, weighted=True)
    def In(self, other):
        import pcd.util
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return pcd.util.In(G1, G2)
    def F1(self, other, weighted=False):
        import pcd.F1
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return .5 * (  pcd.F1.F1(G1, G2)[0]
                     + pcd.F1.F1(G2, G1)[0])
    def F1w(self, other, weighted=False):
        """Like F1, but weighted by community size"""
        import pcd.F1
        G1 = self.to_pcd()
        G2 = other.to_pcd()
        return .5 * (  pcd.F1.F1(G1, G2, weighted=True)[0]
                     + pcd.F1.F1(G2, G1, weighted=True)[0])

