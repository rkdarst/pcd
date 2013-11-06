# Richard Darst, October 2012

import collections
import cPickle as pickle
import itertools
import numpy
import os
import os.path
import random
import re
import subprocess
import time

import networkx

import pcd.nxutil
import pcd.cmty
import pcd.io
import pcd.support.algorithms

import contextlib
@contextlib.contextmanager
def chdir(dirname):
    olddir = os.getcwd()
    os.chdir(dirname)
    yield
    os.chdir(olddir)

def pathsplitall(f):
    """Fully split all componetns of a path"""
    f = os.path.normpath(f)
    a, b = os.path.split(f)
    #print a, b
    if a and b:        return pathsplitall(a) + (b, )
    if a == '/':       return (a, )
    if a and not b:    return pathsplitall(a)
    return (b, )


def read_edgelist(fname, constructor=networkx.DiGraph):
    return networkx.read_edgelist(fname, create_using=constructor())
def read_edgelistreverse(fname, constructor=networkx.DiGraph):
    g = constructor(fname=fname)
    split_re = re.compile('\s+')
    for line  in open(fname):
        line = line.strip()
        if not line: continue
        nodes = split_re.split(line)
        #from fitz import interactnow
        if len(nodes) != 2: raise ValueError('More than two nodes on a line: "%s"'%line)
        node_source = nodes[1]
        node_dest   = nodes[0]
        assert not g.has_edge(node_source, node_dest)
        #assert node_source != node_dest
        g.add_edge(node_source, node_dest)
    return g

def read_pajek(fname, constructor=networkx.DiGraph):
    g = networkx.read_pajek(fname)
    print g.__class__
    # Test for multi-plicitiy
    for node in g.nodes_iter():
        neighbors = g.neighbors(node)
        assert len(neighbors) == len(set(neighbors)), "Not a simple graph..."
    return constructor(g)

#def stats(g):
#    stats = [ ]
#    stats.append("number of nodes: %d"%len(g))
#    stats.append("number of edges: %d"%g.number_of_edges())
#    cmtynodes = pcd.nxutil.cmtynodes(g)
#    stats.append("number of communities: %d"%len(cmtynodes))
#    stats.append("mean community size: %f"%numpy.mean([len(c) for c in cmtynodes.values()]))
#    stats.append("std community size: %f"%numpy.std([len(c) for c in cmtynodes.values()]))
#    cmtys_by_rev_size = [c for (c,ns) in sorted(cmtynodes.iteritems(),
#                                                reverse=True,
#                                                key=lambda (c,ns): len(ns))]
#    stats.append("community sizes: %s"%' '.join(str(len(cmtynodes[c])) for c in cmtys_by_rev_size))
#    stats.append("fraction of nodes in largest community: %f"%(max(len(c) for c in cmtynodes.values())/float(len(g))))
#
#    stats.append("List of communities:")
#    singleton_cmtys = [ ]
#    for c in cmtys_by_rev_size:
#        cmtysize = len(cmtynodes[c])
#        if cmtysize == 1:
#            singleton_cmtys.append(c)
#            continue
#        import textwrap
#        stats.append(textwrap.fill(
#            (("%3s: "%c)+' '.join(str(n) for n in sorted(cmtynodes[c]))),
#            width=120,
#            initial_indent="", subsequent_indent="     ",
#            ))
#    stats.append(textwrap.fill(
#        "Nodes in singleton communities: "+' '.join(
#            list(cmtynodes[c])[0] for c in singleton_cmtys),
#        subsequent_indent="     ",
#        ))
#    return stats

import pcd.support.algorithms
class OslomSeed30_dir(pcd.support.algorithms.Oslom_dir):
    randseed = 30
class OslomSeed73_dir(pcd.support.algorithms.Oslom_dir):
    randseed = 73
class OslomSeed30_undir(pcd.support.algorithms.Oslom_dir):
    randseed = 30
class OslomSeed73_undir(pcd.support.algorithms.Oslom_dir):
    randseed = 73

class BulkAnalyzer(object):
    output = None
    common_output = None
    output_dir = None
    loader = 'edgelist'
    filename_chooser = None

    hook_result = [ ]

    def __init__(self, fname, **kwargs):
        self.fname = fname
        for k, v in kwargs.iteritems():
            if hasattr(self, k):
                setattr(self, k, v)


        # Decide where to write our files.
        if self.filename_chooser:
            # Filename chooser must set these things:
            # - output = basename of output
            self.output = self.filename_chooser(self=self)

        if self.output:
            # if output is given, use that filename.
            self.outname = self.output
        elif self.output_dir:
            # if output_dir is given, use that as a basename for files.
            self.outname = os.path.join(self.output_dir, *pathsplitall(fname)[1:])
            if not os.access(os.path.dirname(self.outname), os.F_OK):
                os.makedirs(os.path.dirname(self.outname))
        else:
            # Otherwise, write to same place with extenision .CD.txt
            self.outname = fname + '.CD.txt'
        print self.outname
        self.results = [ ]
        if isinstance(self.common_output, str):
            self.common_output = open(self.common_output, 'w')

    def load(self):
        if isinstance(self.loader, type(read_edgelist)):
            # if it's a function...
            load_function = self.loader
        else:
            load_function = globals()['read_'+self.loader]
        self.g_dir = load_function(self.fname, constructor=networkx.DiGraph)
        self.g_undir = networkx.Graph(self.g_dir)
        assert len(self.g_dir) == len(self.g_undir), "Success for directed/undirected, remove this line."


    methods = [
        pcd.support.algorithms.Infomap,
        pcd.support.algorithms.InfomapSingle,
        pcd.support.algorithms.Louvain,
        pcd.support.algorithms.Oslom,
        pcd.support.algorithms.Oslom_dir,
        pcd.support.algorithms.COPRA,
        pcd.support.algorithms.COPRA_dir,
        pcd.support.algorithms.ModularitySA,
        pcd.support.algorithms.APM,
        ]
    def runMethods(self):
        for method in self.methods:
            if getattr(method, '_is_directed', False):
                g = self.g_dir.copy()
            else:
                g = self.g_undir.copy()
            self.runMethod(method, g)
    def runMethods_oslom(self):
        self.runMethod(OslomSeed30_dir,     self.g_dir.copy())
        self.runMethod(OslomSeed73_dir,     self.g_dir.copy())
        self.runMethod(OslomSeed30_undir,   self.g_undir.copy())
        self.runMethod(OslomSeed73_undir,   self.g_undir.copy())

    def run(self):
        self.load()
        self.network_output = open(self.outname+'.stats', 'w') # truncate file
        self.results = [ ]

        self.runMethods()
        #self.methods_oslom()

        self.makeGrid('ovIn')
        self.makeGrid('ovInw')
        self.makeGrid('F1')
        self.makeGrid('F1w')
        self.makeGrid('VI')
        self.makeGrid('In')
        pickle.dump(self.resultsDict,
                    open(self.outname+'.results.pickle', 'w'))
        if self.common_output:
            print >> self.common_output, '\n'

    def runMethod(self, Method, g, kwargs={}):
        name = Method.name()
        r = Method(g, basename=self.outname, **kwargs)
        pickle.dump(r.results, open(self.outname+'.'+name+'.pickle','w'))
        open(self.outname+'.result.'+name+'.stats', 'w').close() # truncate file
        self.results += r.results
        self.resultsDict = { }
        for result in r.results:
            #assert len(result.nodecmtys()) == len(result.nodes), "ca1"
            if hasattr(result,'shortlabel'):
                label = result.label
                shortlabel = result.label
                fullname = name+'.'+shortlabel
            elif hasattr(result, 'label'):
                label = result.label
                shortlabel = result.label
                fullname = name+'.'+label
            else:
                label = '-'
                shortlabel = name
                fullname = name
            result.label = fullname

            s =  [ ]
            s += [ "Basename: "+self.outname ]
            s += [ "Analyzer: "+name ]
            s += [ "Result: "+fullname ]
            #s += [ "Filename: "]
            s += result.stats()
            s += result.list_communities()
            s += result.list_overlapping_nodes()
            s += ["", ""]
            s = '\n'.join(s)

            print s
            print >> open(self.outname+'.result.'+name+'.stats', 'a'), s
            if self.network_output:
                print >> self.network_output, s
            if self.common_output:
                print >> self.common_output, s

            #fname = (os.path.join(r.dirname, 'result.'+result.name+'.txt'))
            self.resultsDict[fullname] = result

            # Write these communities out
            fname = self.outname + '.result.'+fullname+'.txt'
            result.write_clusters(fname, headers=['Result-Name: %s'%fullname])

            # Draw pictures
            if 1 or self.draw_nodes:
                import pcd.draw
                if getattr(self, 'pos', None) is None:
                    self.pos = pcd.draw.layout(g)
                fname = self.outname + '.result.'+fullname+'.png'
                pcd.draw.draw(g, fname, result, pos=self.pos, figsize=(50,50))
            # Hooks for results
            for hook in self.hook_result:
                hook(self=self, g=g, cmtys=result,
                     basename=self.outname+'.result.'+fullname,
                     resultname=self.outname+'.result.'+name)


    def makeGrid(self, measure, **kwargs):
        results = self.results
        nResults = len(self.results)
        matrix = numpy.zeros(shape=(nResults, nResults),
                             dtype=float)
        for i in range(nResults):
            cmtys1 = results[i]
            for j in range(i+1, nResults):
                cmtys2 = results[j]
                if cmtys1.cmtysizes_sum()==0 or cmtys2.cmtysizes_sum()==0:
                    # The measures break down for n=0 for all
                    # communities, ignore these cases.
                    matrix[i, j] = float('nan')
                    continue
                if measure == 'ovIn':
                    # use the external code for this measure.  FIXME:
                    # find difference between my code and external
                    # code.
                    matrix[i,j] = cmtys1.ovIn_LF(cmtys2)
                else:
                    matrix[i,j] = getattr(cmtys1, measure)(cmtys2)
        s = [ ]
        s.append("Comparison of different community detections, using "+measure)
        #print matrix
        s.append("")
        s[-1] += "      "
        for j in range(nResults):
            s[-1] += "%4d  "%j
        s.append("")
        for i in range(nResults):
            s[-1] += "%4d  "%i
            s[-1] += ''.join("%6.2f"%matrix[i,j] if matrix[i,j] else '   -  '
                          for j in range(nResults))
            s.append("")
        for i, r in enumerate(self.results):
            s.append("%4d: %s"%(i, r.label))
        s.append('')
        s.append('')
        s = '\n'.join(s)

        print s
        if self.network_output:
            print >> self.network_output, s
        if self.common_output:
            print >> self.common_output, s


    def stats(self, name):
        raise
        s = '\n'.join(stats(self.g))
        print >> open(self.outname+'.'+name+'.stats', 'w'), s
        print s

        if self.common_output:
            print >> self.common_output, self.outname
            print >> self.common_output, "Analyzer:", name
            print >> self.common_output, s
            self.common_output.flush()

        #import matplotlib
        #matplotlib.use('Agg')

        #import
        #cm = matplotlib.colors.Colormap('jet')()
        #import matplotlib.cm
        #cm = matplotlib.cm.get_cmap('jet')
        #
        #colors = dict((node, tuple(cm(self.g.node[node]['cmty']))) for node in self.g.nodes())
        #print colors
        #networkx.draw_random(self.g,
        #              #node_color=colors
        #              )
        #import matplotlib.pyplot as plt
        #plt.savefig(self.outname+'.infomap.png')

def hook_result_pajek(g, cmtys, basename, **kwargs):
    #g = g.copy()
    #cmtys.load_networkx_custom(g, type_=str)
    #networkx.write_pajek(g, basename+'.net')
    pcd.io.write_pajek(basename+'.net', g, cmtys)
    #raise

def output_subdir(self, base='output/', cut_dirs=1):
    fname = self.fname
    fname_parts = pathsplitall(fname)[cut_dirs:]
    output = os.path.join(*([base] + list(fname_parts) + [fname_parts[-1]]))
    if not os.access(os.path.dirname(output), os.F_OK):
        os.makedirs(os.path.dirname(output))
    return output

if __name__ == "__main__":
    from fitz import cliargs
    args, options = cliargs.get_cliargs()
    #global_output_dir = "output-oslomseed"
    #global_output_dir = "output-test"
    #global_output_dir = options['output']


    common_output = None
    if 'allresults' in options:
        all_results = options['allresults']
        if all_results is True:
            all_results = "allresults.%05d.txt"

        if '%' in all_results:
            # If a % is in the filename, substiute it with an incrementing
            # counter.
            for i in itertools.count():
                fname = all_results%i
                if not os.access(fname, os.F_OK): break
        common_output = open(fname, 'a')
        print >> common_output, ("==== "+time.ctime()+' ').ljust(80, '=')

    BulkAnalyzer.methods = [
        pcd.support.algorithms.Infomap,
        pcd.support.algorithms.InfomapSingle,
        pcd.support.algorithms.Louvain,
        pcd.support.algorithms.Oslom,
        pcd.support.algorithms.Oslom_dir,
        pcd.support.algorithms.COPRA,
        pcd.support.algorithms.COPRA_dir,
        pcd.support.algorithms.ModularitySA,
        #pcd.support.algorithms.APM,
        ]

    #BulkAnalyzer(args[1], output=options['output'], output_dir=options.get('outputdir', None),
    #             common_output=common_output).run()

    infname = args[1]
    outfname = infname+'.d/'+infname
    if not os.access(os.path.dirname(outfname), os.F_OK):
        os.mkdir(os.path.dirname(outfname))
    BulkAnalyzer(infname, output=infname+'.d/'+infname, output_dir=options.get('outputdir', None),
                 common_output=common_output).run()


    #for d in os.listdir('data/'):
    #    if os.access(os.D_OK)
