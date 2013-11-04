# Richard Darst, October 2012

import collections
import itertools
import numpy
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time

import networkx

import pcd.nxutil
import pcd.cmty
from pcd.util import chdir_context

global_code_path = ["/home/richard/research/cd/cd-code/",
                    "/proj/networks/darst/cd-code/",
                    ]
if os.uname()[1] == 'amor':
    pass

def _get_file(name):
    """File opener, using different search paths.

    This is used to find the external programs called by this module.

    - relative to os.getcwd(),
    - relative to the directory of this file (algorithms.py)
    - relative to the parent directory of this file (algorithms.py)
    - relative to the contents of the global variable in this module
      `global_code_path`
    - relative to the colon-separated contents of the environment variable
      PCD_EXTERNAL_CODE_PATH
    """
    search_path = [
        '',  # search globally first.
        os.path.dirname(__file__),
        os.path.dirname(os.path.dirname(__file__)),
        ]
    search_path.extend(global_code_path)
    if 'PCD_EXTERNAL_CODE_PATH' in os.environ:
        search_path.extend(os.environ[PCD_EXTERNAL_CODE_PATH].split(':'))
    for dirname in search_path:
        path = os.path.join(dirname, name)
        if os.access(path, os.F_OK):
            return path
    raise ValueError("Can not find program: %s"%name)


class GraphFile(object):
    """Class representing already-exsting edge list, to use as input
    to the algorithm."""
    def __init__(self, fname, vmap_inv=None, index=0):
        self.fname = fname
        self.vmap_inv = None
        self.index = 0
    def setup(self, meth):
        if not os.access(meth.graphfile, os.F_OK):
            os.symlink(self.fname, meth.graphfile)

def pathsplitall(f):
    """Fully split all componetns of a path"""
    f = os.path.normpath(f)
    a, b = os.path.split(f)
    #print a, b
    if a and b:        return pathsplitall(a) + (b, )
    if a == '/':       return (a, )
    if a and not b:    return pathsplitall(a)
    return (b, )


def read_graph(fname, constructor=networkx.DiGraph):
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

def _identity(x):
    return x
class NullMap(object):
    """Class to use as a placeholder when node names do NOT need to be mapped."""
    def __init__(self, map=_identity, unmap=int):
        self.map  = map
        self.unmap = unmap
    def __getitem__(self, x):
        return self.map(x)
    def __call__(self, iterable):
        """Unmap entire set of items - return generator, so wrap in object creator."""
        map = self.map
        return (map[x] for x in iterable)
    def inverse(self):
        import copy
        inv = copy.copy(self)
        inv.map, inv.unmap = self.unmap, self.map
        return inv

class CDMethod(object):
    """Base class to present standard interface to CD methods.

    The CD interface is defined as such.

    A method is initialized with arguments g, a networkix graph
    """
    delete_tmpdir = None         # Delete the temporary directory after running?
    _nodemapZeroIndexed = False   # Should edgelist be zero-indexed?  default start at 1.
    _interact_results = False     # If True, after self.run() is called, interact in directory.
    seed = None                   # Seed to override.
    verbosity = 2
    map_nodes_to_int = True       # Do we need to remap all nodes to ints?
    _run_method = True
    _load_communities = True

    def __init__(self, g, dir=None, basename=None, **kwargs):
        """
        Arguments:

        - g: a networkx graph object.  This will be automatically
        written to disk in the proper location, and with the proper
        type

        - dir: chdir to this directory before running, default
          tmp-GRAPH-CDMETHOD.RANDOM.  This directory is saved after
          the running, unless delete_tmpdir is True.

        - basename: root of all input and output filenames within the
          temporary directory.  Default 'graph'.

        - **kwargs: any other keyword arguments are saved on self, as
          attributes.  The attributes MUST be already defined in the
          class, in order to prevent mistakes.

        For unique runs, you should use a unique tmpdir.  Some methods
        do not use the basename as a prefix to all output, so if run
        in the same directory then they will override each other's
        files.

        g will not be modified.

        How methods are run:

        - self.gen_filenames() is called to figure out the temporary
          directory name, and graph file name.

        - Within the temporary directory,

          - self._input_format determines which method is called to
            write the file names: getattr(self,
            'write_'+self.`self._input_format`) and it is called with
            the arguments (g, graphfile)

          - self.run() is called.  This is the central entry point to
            all CD methods, and is the only thing that needs
            overriding in subclasses.  This method simply uses
            self.graphfile to run a CD method, and parses the results
            and stores them in self.cmtys (for single-best results),
            and self.results for a list of all hierarchical
            levels/partitions/layers returned.

        """
        self.g = g
        self.dir = dir
        self.basename = basename
        self.gen_filenames()
        # **kwargs sets options (attributes on the class)
        for k, v in kwargs.iteritems():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise ValueError("Unknown option %s"%k)

        # Make mapping of nodes to integers, if g is a GraphFile
        # (already existing file), it won't make a map.
        self.make_map()
        with chdir_context(self.dir):
            # Get the respective writing method and run it.
            if isinstance(g, GraphFile):
                self.g.setup(self)
            else:
                # if self.graphfile is a symlink, it might be pointing
                # to some other important file that would then be
                # overwritten.
                if os.path.islink(self.graphfile):
                    raise RuntimeError("%s is a symlink, aborting."%self.graphfile)
                getattr(self, 'write_'+self._input_format)(self.g, self.graphfile)
            if self._run_method:
                self.run()
            if self._load_communities:
                self.read_cmtys()
            if self._interact_results:
                import code; code.interact(banner=("In dir=%s\n"
                                   "File is %s")%(self.dir, self.graphfile),
                                   local=locals())
        if self.delete_tmpdir:
            shutil.rmtree(self.dir)

    # Utility methods
    @classmethod
    def name(cls):
        """Return CD method name, as a string"""
        return getattr(cls, '_name', cls.__name__)
    @property
    def randseed(self):
        """Random (integer) seed for CD methods.

        We always pass the method a seed generated from python
        random.randint, when possible.  This is a better seed than the
        time, it uses os.urandom when possible."""
        if self.seed is not None:
            return self.seed
        # Seed is from os.urandom.
        return random.randint(0, 2**31-1)
    def make_map(self):
        """Create a mapping of nodes -> integers.

        This is conrolled by several things:"""
        # If we want no mapping: have a placeholder that does nothing
        # when mapping, converts to integers when unmapping.
        if not self.map_nodes_to_int or isinstance(self.g, GraphFile):
            self.vmap = NullMap()
            self.vmap_inv = self.vmap.inverse()
            return
        # Do we want indexes to start at zero or one?
        if self._nodemapZeroIndexed:
            start_at = 0
        else:
            start_at = 1
        # Are we already in the proper format?  Test if the set of all
        # nodes is equal to the set of ints from start_at, start_at+N.
        # There is still one problem here: we don't know test if x's
        # are all ints, or strings representing ints.  I ignore this
        # problem for now.  FIXME.
        try:
            node_set = set(int(x) for x in self.g.nodes_iter())
            int_set = set(range(start_at, len(self.g)+start_at))
        except ValueError:
            node_set = 1
            int_set = 0
        if node_set == int_set:
            # This will convert resultant values to Python ints.
            self.vmap = NullMap()
            self.vmap_inv = self.vmap.inverse()
            return
        # Otherwise, map all nodes to integers:
        self.vmap = vmap = { }
        self.vmap_inv = vmap_inv = { }
        for i, n in enumerate(self.g.nodes_iter()):
            i += start_at
            vmap[n] = i
            vmap_inv[i] = n
        #print vmap, vmap_inv
    def gen_filenames(self):
        """Generate filenames for running:

        This method sets self.dir and self.graphfile to be directories
        to run in.  Here is our algorithm:

        If neither `dir` nor `basename` are given, use a
        basename of 'graph' and create a temporary directory of
        'tmp-graph.CDMETHODNAME.RANDOM.

        If `dir` is given, run in that directory with a basename of
        `graph` within that directory.

        If basename is given, create a temporary directory of name
        BASENAME.CDMETHODNAME.d and run within that directory.

        How are these names used?

        - Before the method runs, we chdir into `self.dir`
          (automatically).  All filenames are relative to this
          directory.
        - Within that directory, `self.graphfile` is the name of the
          input file.  This is basename + self.graphfileExtension
          (default '.net').
        """
        #if not self.dir:
        #    self.dir = tempfile.mkdtemp()
        if not self.dir and not self.basename:
            # need both basename and dir.
            self.basename = 'graph'
            self.dir = tempfile.mkdtemp(
                prefix=''.join(('tmp-', self.basename, '.', self.name(),'-')),
                dir='.', )
            if self.delete_tmpdir is None:
                self.delete_tmpdir = True
        if not self.basename:
            # have dir, need basename.
            self.basename = 'graph'
        if not self.dir:
            # have basename, need dir.
            self.dir = self.basename+'-'+self.name()+'.d'
            self.basename = os.path.basename(self.basename)
            #self.basename = os.path.basename(self.basename)
        self._absdir = os.path.abspath(self.dir)
        self.graphfile = self.basename+getattr(self, 'graphfileExtension', '.net')
        if not os.access(self.dir, os.F_OK): os.mkdir(self.dir)
    def call_process(self, args, binary_name=None):
        """Call a process, saving the standard output to disk.

        This handles running a program, while also doing other
        important tasks like saving the times, command line arguments,
        return value, and stdout to disk.
        """
        if self.verbosity >= 2:
            print " ".join(args)
        # Stdout is saved to BINARY_NAME.stdout, BINARY_NAME defaults
        # to args[0]
        if binary_name is None:
            binary_name = os.path.basename(args[0])
        binary_name = binary_name.replace('/', '_')
        # Open a stream for the standard output
        self._stdout_fname = binary_name+'.stdout'
        disk_stdout = open(binary_name+'.stdout', 'w')
        # Write the command lines we are running
        disk_stdout.write(repr(list(args))+'\n')
        disk_stdout.write("self attributes: %s\n"%(pcd.util.listAttributes(self),))
        disk_stdout.write("Process times: %s\n"%(os.times(),))
        disk_stdout.write("Beginning run at %s\n\n\n"%time.ctime())
        disk_stdout.flush()

        if self.verbosity >= 2:
            # In order to split the subprocess output between two
            # files (stdout and disk), we have to resort to deep unix
            # trickery.  There is simply no other way to "tee" it
            # without another thread.  (edit: perhaps you could do a
            # loop reading pipe.stdout...)
            from pcd.support.fdop import FDSplitter
            stdout = splitter = FDSplitter([disk_stdout, sys.stdout])
        else:
            stdout = disk_stdout

        ret = None
        try:
            ret = subprocess.call(args, stdout=stdout, stderr=subprocess.STDOUT)
        finally:
            #if self.verbosity >= 2:
            #    print open(binary_name+'.stdout').read()
            if self.verbosity >= 2:
                splitter.close()
                splitter.wait()

            disk_stdout.write('\n\n\n')
            disk_stdout.write("Ending run at %s\n"%time.ctime())
            disk_stdout.write("Process times: %s\n"%(os.times(),))
            disk_stdout.write("Return value: %s\n"%ret)
            disk_stdout.flush()
        if ret != 0:
            raise ValueError("CD method %s (%s) command failed: %s."%(
                self.name(), self.graphfile, ret))
        return ret

    # Format writer methods
    def write_pajek(self, g, fname):
        f = open(fname, 'w')
        print >> f, '*Vertices %d'%len(g)
        vmap = self.vmap
        for i, n in enumerate(g.nodes_iter()):
            #if ' ' in repr(n): raise ValueError("Node '%s' has a space in its name"%n)
            print >> f, i, '%s'%vmap[n]#, "0.00", "0.00", "ellipse"
            #print >> f, i, repr(n)#, "0.00", "0.00", "ellipse"
        print >> f, "*Edges"
        for u, v, d in g.edges_iter(data=True):
            if 'weight' in d:
                print >> f, vmap[u], vmap[v], float(d['weight'])
            else:
                print >> f, vmap[u], vmap[v] #, "1.0"
        f.flush()
        f.close()

    def write_edgelist(self, g, fname):
        f = open(fname, 'w')
        vmap = self.vmap
        for u, v, d in g.edges_iter(data=True):
            #if getattr(self, 'weighted', True):
            if 'weight' in d:
                # weighted
                print >> f, vmap[u], vmap[v], float(d.get('weight', 1.0))
            else:
                # unweighted
                print >> f, vmap[u], vmap[v]
    def write_gml(self, g, fname):
        # networkx writes nodes in order of g.nodes()
        #self.node_order = g.nodes()
        #for i, n in enumerate(g.nodes_iter()):
        #    g.node[n]['id'] =
        #assert [ int(x) for x in g.nodes() ]
        g = networkx.relabel_nodes(g, mapping=self.vmap)  # Make a copy.
        networkx.write_gml(g, fname)
    def write_null(self, g, fname):
        pass


    # To create your own CD methods, override 
    def run(self):
        """Run the actual CD method"""
        raise NotImplementedError()
    def read_cmtys(self):
        """This method reads and saves results.

        This method should have two side-effects:

        - Read the 'best' community structure and save it in
          self.cmtys.

        - Save all possible results to a list self.results.

        This method should return self.cmtys"""
        return self.results
    def read_cmtys_iter(self):
        return self.results


class NullCD(CDMethod):
    """Null community detection: always returns all nodes in one community."""
    _input_format = 'null'
    def run(self):
        cmtys = pcd.cmty.Communities(cmtynodes={0:set(self.g.nodes())})
        cmtys.label = '(q=1)'
        self.results = [ cmtys ]
        self.cmtys = cmtys

class _InfomapOld_(CDMethod):
    """Old non-hierarchical infomap code."""
    trials = 50
    _input_format = 'pajek'
    def run(self):
        args = (#'/usr/bin/gdb', '--args',
                _get_file(self.binary),
                str(self.randseed),
                self.graphfile,
                str(self.trials))
        self.call_process(args)

    def read_cmtys(self):
        #if isinstance(self, Infomap_undir):
        #    out_filename = self.basename.split('.')[0]+'.tree'  # Infomap splits incorrectly
        #else:
        #    out_filename = self.basename+'.tree'
        fname = self.basename+'.tree'

        linere = re.compile(r'(\d+):(\d+) ([\d.e-]+) "([^"]+)"')
        cmtynodes = collections.defaultdict(set)

        for line in open(fname):
            if line[0] == '#': continue
            if not line.strip(): continue
            m = linere.match(line)
            cmty = int(m.group(1))-1
            nodeCmtyID = int(m.group(2))
            #weight = float(m.group(3))
            nodeID = int(m.group(4))
            #cmtynodes[cmty].add(nodename)
            cmtynodes[cmty].add(self.vmap_inv[nodeID])

        self.cmtys = pcd.cmty.Communities(dict(cmtynodes))
        self.results = [ self.cmtys ]
        return self.results

    def run(self):
        self.run_infomap()


class _InfomapOld(_InfomapOld_):
    binary = 'infomap/infomap_undir/infomap'

class _InfomapOld_dir(_InfomapOld_):
    binary = 'infomap/infomap_dir/infomap'
    _is_directed = True

class Infomap(CDMethod):
    """New infomap code.

    This runs Infomap code from http://mapequation.org/code.html.

    hierarchical: If true, use hierarchical infomap code.  If false,
    return just one partition of the network.  Note that this
    corresponds to the '--two-level' option to the Infomap binary.

    directed: If true, have infomap treat the network as directed.

    full_lower_levels: If true(default), each layer is a complete partition
    of the system.  If false, each lower level includes only the
    communities that are different from the higher levels.

    trials: number of trials

    args: list of strings to be added to the arguments.

    initial: initial configuration.
    """
    binary = 'infomap/Infomap-0.11.5/Infomap'
    _input_format = 'edgelist'
    _nodemapZeroIndexed = True
    directed = False
    hierarchical = True
    full_lower_levels = True
    initial = None
    trials = 10
    include_self_links = False
    args = [ ]
    def run(self):
        args = [_get_file(self.binary),
                '--num-trials=%d'%self.trials,
                '--input-format=link-list', '--zero-based-numbering',
                '--seed=%d'%self.randseed,
                self.graphfile,
                '.',
                ]
        if self.directed:
            args.append('--directed')
        if not self.hierarchical:
            # This argument says two-level, but the lowest level is
            # always just singleton nodes.
            args.append('--two-level')
        if self.include_self_links:
            args.append('--include-self-links')
        if self.initial:
            self.initial.write_clusters("initial.txt", mapping=self.vmap)
            args.append('--cluster-data=initial.txt')
        args.extend(self.args)
        self.call_process(args)

    def read_cmtys(self):
        fname = self.basename+'.tree'
        self.results = self.read_infomap_cmtys(
            fname=fname, vmap_inv=self.vmap_inv,
            full_lower_levels=self.full_lower_levels)
        self.cmtys = self.results[0]
        return self.results

    def read_infomap_cmtys(self, fname, full_lower_levels=True, vmap_inv=None):
        """General function for reading Infomap tree files.

        This is a static method, so can easily be used in other code
        separate from this framework.

        WARNING: Infomap writes an output with node indexes one-based,
        even if the input is zero-based.  We subtract one from each
        output node ID to correct for that."""
        # Correct for Infomap changing zero-based to one-based
        # indexes.  If we don't have self._nodemapZeroIndexed, then
        # the code below must be corrected (search "zero" below).
        assert self._nodemapZeroIndexed

        linere = re.compile(r'([\d:]+) ([\d.e-]+) "([^"]+)"')
        #cmtynodes = collections.defaultdict(set)

        #raise NotImplemented()

        levels = [ ]
        all_nodes = [ ]
        for line in open(fname):
            if line[0] == '#': continue
            if not line.strip(): continue
            m = linere.match(line)
            # cmtys is a node level0:level1:level2:...:index
            # level0 is topmost level, the last level is an individual
            # node index within that level.
            cmtys = m.group(1).split(':')
            # index is last item.  A hierarchical community ID can be
            # constructed by using everything except that last index.
            #community_id = ':'.join(cmtys[:-1])
            #cmty = int(m.group(1))-1
            #nodeCmtyID = int(m.group(2))
            #weight = float(m.group(3))
            nodeID = int(m.group(3))
            # Correct for Infomap changing zero-based to one-based
            # indexes.  The assertion makes sure that future infomap
            # behavior does not change, and begin outputting
            # zero-based indexes.
            nodeID -= 1
            assert 0 <= nodeID, "Infomap indexing behavior may have changed, see source."
            #
            node = vmap_inv[nodeID]
            for level in range(len(cmtys)-1):
                # for each level and community, but excluding the
                # unique index on the end...
                if len(levels) <= level:
                    levels.append(dict())
                levels[level][node] = ':'.join(cmtys[:level+1])  # cmtys is the full hierarchy:

        # Now that we have a list of cmtynodes dicts set up, so we go
        # through and generate the "real" communities.
        if full_lower_levels:
            for l0 in range(len(levels)):
                for node in levels[0]:
                    # search level.
                    for l1 in range(l0, -1, -1):
                        if node in levels[l1]:
                            levels[l0][node] = levels[l1][node]
                            break
        # Reverse, to put most granular level first.
        levels.reverse()
        results = [ ]
        for l in range(len(levels)):
            cmtys = pcd.cmty.Communities.from_nodecmtys(levels[l])
            cmtys.label = 'level%02d'%l
            results.append(cmtys)
        # Test for zero-based indexing fixing.
        #if isinstance(self.g, networkx.Graph):
        #    assert set(cmtys.nodes) == set(self.g.nodes_iter())

        return results

class Infomap_dir(Infomap):
    """Subclass of Infomap, options set for directed analysis."""
    _is_directed = True
    directed = True

class InfomapSingle(Infomap):
    """Subclass of Infomap, options set for one-layer."""
    hierarchical = False

class InfomapSingle_dir(Infomap):
    """Subclass of Infomap, options set for one-layer, and directed."""
    hierarchical = False
    _is_directed = True
    directed = True



class _Oslom(CDMethod):
    #binary = 'oslom/OSLOM2/oslom_undir'
    _input_format = 'edgelist'
    trials = 10

    no_singletons = False
    initial = None
    weighted = None
    strip_singletons = False

    #def _use_weigthed(self, g):
    #    if self.weigthed is not None:
    #        return self.weighted

    def run(self):
        self.outdir = self.graphfile + '_oslo_files/'
        # Delete all pre-existing output files - if we don't, there
        # might be higher hierarchical levels that aren't deleted from
        # a previous run, giving us results mixed from this run and
        # previous runs.
        if os.access(self.outdir, os.F_OK):
            files = os.listdir(self.outdir)
            for f in files:
                os.unlink(os.path.join(self.outdir, f))

        args = (_get_file(self.binary),
                "-seed", str(self.randseed),
                "-w" if self.weighted else '-uw', #unweighted or weighted
                "-f", self.graphfile,
                '-r', str(self.trials),
                )
        if self.no_singletons:
            # In the new code, this is called "-singlet" and has the
            # opposite effect, so negate the conditional.
            args = args + ('-all', )
        if self.initial:
            args = args + ('-hint', "initial_state.txt",)
            f = open("initial_state.txt", 'w')
            for cmty, nodes in self.initial.cmtynodes():
                f.write(' '.join(str(self.vmap[n]) for n in nodes))
                f.write('\n')
            f.close()

        self.call_process(args)

    def read_cmtys(self):
        self.results = [ ]
        for level in itertools.count():
            fname = self.outdir+'tp%s'%(level if level else '')
            if not os.access(fname, os.F_OK):
                if self.verbosity >= 3:
                    print "*** level %s"%level
                break

            cmtys = self.read_oslom_cmtys(fname=fname)
            cmtys.label = "level%02d"%level

            self.results.append(cmtys)
            if not self.strip_singletons:
                assert len(self.g) <= self.results[-1].cmtysizes_sum()
        self.cmtys = self.results[0]
        return self.results
    def read_oslom_cmtys(self, fname):
        """Read Oslom communities from a single hierachical level"""
        f = open(fname)
        cmtynodes = { }
        for i, line in enumerate(f):
            if i%2 == 0:
                m = re.match(
                    r'#module ([\d]+) size: (\d+) bs: ([\d0.e-]+)',
                    line)
                cmty = int(m.group(1))
            else:
                nodes = [ int(x) for x in line.split() ]
                nodes = set(self.vmap_inv[n] for n in nodes)
                #print nodes
                if self.strip_singletons and len(nodes) == 1:
                    continue
                cmtynodes[cmty] = set(nodes)
        # If we detected all singletons, we have no community
        # structure.  Return one community covering the whole graph,
        # since we tend to have problems with null-graphs.
        if len(cmtynodes) == 0:
            cmtynodes[0] = set(self.g.nodes())
        # If we can strip singletons, then we can end up not spanning
        # the whole network.  In that case we must specify nodes= to
        # the communities object, or else it won't be able to
        # recalculate the universe correctly when asked.
        nodes = None
        if self.strip_singletons:
            nodes = set(self.g.nodes())
        return pcd.cmty.Communities(cmtynodes,
                                    nodes=nodes,
                                    )

class Oslom(_Oslom):
    binary = 'oslom/OSLOM2/oslom_undir'
class OslomWeighted(Oslom):
    weighted = True
class Oslom_dir(_Oslom):
    binary = 'oslom/OSLOM2/oslom_dir'
    _is_directed = True


class _Copra(CDMethod):
    # http://www.cs.bris.ac.uk/~steve/networks/software/copra-guide.pdf
    _input_format = 'edgelist'
    trials = 10
    binary = 'copra.jar'
    weighted = False
    max_overlap = 1   # default in COPRA code is 1.
    max_overlap_range = None

    def run(self):
        args = ["/usr/bin/java",
                "-cp", _get_file(self.binary),
                "COPRA",
                os.path.abspath(self.graphfile),
                "-repeat", str(self.trials),
                ]
        #print " ".join(args)
        #os.spawnv(os.P_WAIT, "/usr/bin/java", args)
        if self.weighted:
            args += ['-w', ]
        if self.max_overlap:
            args += ['-v', str(self.max_overlap)]
        if self.max_overlap_range:
            raise NotImplementedError("You can specify this option, but I don't properly parse the output yet.")
            args += ['-vs', str(self.max_overlap_range[0]), str(self.max_overlap_range[1])]
        self.call_process(args)

    def read_cmtys(self, fname):
        cmtys = self.read_cmtys('clusters-'+os.path.basename(self.graphfile))
        cmtynodes = collections.defaultdict(set)
        for cmty, line in enumerate(open(fname)):
            nodes = [ int(x) for x in line.split() ]
            cmtynodes[cmty] = set(self.vmap_inv[i] for i in nodes)
            #print line
            #print cmtynodes
        self.cmtys = pcd.cmty.Communities(cmtynodes)
        self.results = [ self.cmtys ]
        return self.results

class Copra(_Copra):
    pass
class Copra_dir(_Copra):
    _is_directed=True
    pass
COPRA = Copra
COPRA_dir = Copra_dir
class CopraWeighted(_Copra):
    weigthed = True

class APM(CDMethod):
    _input_format = 'null'
    options = { }  # All options for the pcd.auto.Auto object.
    preferred_order = 'In', 'VI'
    initial = None
    map_nodes_to_int = False

    def run(self):
        import pcd.multiresolution
        try:
            self.run_APM()
        except pcd.multiresolution.NoResultsError:
            self.cmtys = pcd.cmty.Communities(cmtynodes={0:set(self.g.nodes())})
            self.cmtys.label = "NullResults"
            self.results = [ self.cmtys ]
    def run_APM(self):
        import pcd.auto

        initial = self.initial
        if initial is None:
            initial = 'random'
        if isinstance(initial, pcd.cmty.Communities):
            initial = initial.to_pcd().getcmtystate()
            raise NotImplementedError("to_pcd doesn't use the same map as the otherG.  Fix this.")

        # Contruct options dictionary.
        defaultOptions = dict(overlap=True,
                              threads=6,
                              gammas=dict(low=.01, high=10, density=20),
                              initial=initial,
                              quiet=(self.verbosity <= 0))
        defaultOptions.update(self.options)
        options = defaultOptions
        options.update(dict(basename='graph',
                            plot=False))

        # At this point we are chdir'ed into the directory.
        a = pcd.auto.Auto(options=options)
        a.run_g(self.g)
        self.auto = a

        self.results = [ ]
        self.results_d = { }
        for overlap in (True, False):
            if overlap  and  not options['overlap']:
                continue
            # Sort in order of "preferred measures", lowest measures =
            # most preferred, so that first element is the default
            # return value.
            best_cmtys = a.best_communities(overlap=overlap)
            #from fitz import interactnow
            preferred_order = dict((v,k) for k,v in enumerate(self.preferred_order))
            keyfunc = lambda x: preferred_order.get(x.varname, len(preferred_order))
            best_cmtys.sort(key=keyfunc)
            #for cc in best_cmtys: ####
            #    print sorted(cc.cmtysizes().values(), reverse=True)  ####
            self.results[:] = best_cmtys
            results = self.results

            if overlap == options['overlap']:
                self.cmtys = results[0]

            #results = \
            #      list(v for k,v in sorted(best_cmtys.iteritems()))
            # Reorder things so that the "N" (renamed to "ovIn") measure
            # is first.
            #for i in range(len(results)):
            #    result = results[i]
            #    if result.shortname in ('N', 'ovIn'):
            #        results.insert(0, results.pop(i))
            #        break
            # Make names more descriptive, include "In"
            if overlap: overlapName = '(overlap) '   ; ovNameShort = 'ov '
            else:       overlapName = '(no overlap) '; ovNameShort = ''
            for r in results:
                # If this is the "preferred result", then save it in
                # self.cmtys.
                #if r.shortname == self.preferred_result and overlap == options['overlap']:
                #    self.cmtys = r
                self.results_d[r.name] = r
                r.label = overlapName+r.name
                r.shortlabel = overlapName+r.shortname
            self.results += results


        # save results as self.cmtys and self.results
class APM_dir(APM):
    _is_directed = True
    pass


class Louvain(CDMethod):
    """Louvain community detection method.

    weighted: if True, runs in weigthed mode (-w option)

    stop_dQ: float, optional: if given stop when modularity change is
    this much.

    trials: int, default 5: do this many trial runs.  This is
    implemented in python, not in Louvain itself.

    which_partition: All partitions are stored in self.results.  One
    partition is stored in self.cmtys.  which_partition is either an
    index (0, 1, ... -1) of the results list.  0 is the most granular
    (smallest communty size), and -1 is the least granular (largest
    community size).  which_partition can also be the string 'modmax',
    in which case the partition of maximum modularity is taken.  If we
    have multiple trials, return the set of partition with the maximum
    modularity of this level.


    At the end, the class has these attributes as information:

    num_levels: number of levels in the results returned.

    avg_num_levels: number of levels detected, average over all trials.

    modularity: modularity of the partition returned (self.cmtys)

    results_modularity: modularity of each partition in self.results.

    """
    _input_format = 'edgelist'
    binary_convert = 'louvain/Community_latest/convert'
    binary_community = 'louvain/Community_latest/community'
    weighted = False
    stop_dQ = None
    _nodemapZeroIndexed = True
    trials = 10
    # which_partition: which of the louvain partations should we
    # return?  If this is 'modmax', return the partition of maximum
    # modularity out of all trials, and all levels in all trials.  If
    # this is an integer, return the ith (0=first, most granular,
    # -1=last, largest cmtys) level of detection, BUT, the one (out of
    # all the self.trials attempts) of maximum Modularity.
    #which_partition = 'modmax'
    which_partition = 0
    initial = None

    def run(self):
        """Run a number of Louvain trials, and use the best of them."""
        best_cmtys = None
        best_Q = -1e6
        self.avg_num_levels = []
        for i in range(self.trials):
            self.run_louvain()
            results = self.read_cmtys_and_return()
            assert len(results) > 0
            self.avg_num_levels.append(len(results))
            # Pick which partition to use
            if self.which_partition == 'modmax':
                best_of_results = max(results, key=lambda x: x.modularity)
            else:
                best_of_results = results[self.which_partition]
            #print results
            if best_of_results.modularity > best_Q or (
                  best_of_results.modularity is None and not hasattr(self, 'cmtys')
                ):
                self.results = results
                self.results_modularity = [r.modularity for r in results]
                self.cmtys = best_of_results #self.results[0]
                self.modularity = best_Q = results[0].modularity
            #print best_Q
        self.num_levels = len(self.results)
        self.avg_num_levels = numpy.mean(self.avg_num_levels)

    def run_louvain(self):
        self.graphfilebin = self.graphfile+'.bin'
        #Remove old files.
        for fname in (self.graphfilebin, self.graphfilebin+'.weights'):
            if os.access(fname, os.F_OK):
                os.unlink(fname)
        args = (_get_file(self.binary_convert),
                '-i', self.graphfile, '-o', self.graphfilebin,
                )
        if self.weighted:
            args = args + ('-w', self.graphfilebin+'.weights', )
        self.call_process(args)

        args = (_get_file(self.binary_community),
                self.graphfilebin,
                '-l', '-1',
                '-v', )
        if self.weighted:
            args = args + ('-w', self.graphfilebin+'.weights', )
        if self.stop_dQ:
            args = args + ('-q', str(self.stop_dQ),)
        if self.initial:
            raise NotImplementedError("Check this before using.")
            # FIXME: write it
            args = args + ('-p', 'initial.txt')

        self.call_process(args)

    def read_cmtys(self):
        self.results = self.read_cmtys_and_return()
        self.cmtys = self.results[0]
        return self.results
    def read_cmtys_and_return(self):
        stdout = open(self._stdout_fname).read()

        results = [ ]
        level = None
        modularity = None
        #print self.stdout
        for line in stdout.split('\n'):
            #print 'line:', line
            # Find 'level' starting lines
            r = re.match(r'level ([0-9]+):', line)
            if r:
                level = int(r.group(1))
                # Initialization
                cmtynodes = collections.defaultdict(set)
            # If not in a 'level' block, then break.
            if level is None: continue
            # End condition
            if re.match(r'  end computation: ', line):
                cmtys = pcd.cmty.Communities(dict(cmtynodes))
                cmtys.modularity = modularity
                cmtys.label = 'level%02d'%level
                #level += 1
                results.append(cmtys)
                level = None
                del cmtynodes, modularity, cmtys
                modularity = None
            # Modularity:
            r = re.match(r'  modularity increased from [0-9.+-]+ to ([0-9.+-]+)', line)
            if r:
                modularity = float(r.group(1))
            # Actual node parsing
            r = re.match(r'([\d]+) ([\d]+)', line)
            if r:
                node, c = int(r.group(1)), int(r.group(2))

                def find_real_node(node, level):
                    #print 'x', repr(node), repr(level)
                    if level == 0:
                        return set((self.vmap_inv[node],) )
                    assert node in results[level-1]._cmtynodes
                    lower_nodes = results[level-1]._cmtynodes[node]
                    return lower_nodes
                nodes = find_real_node(node, level)
                #print 'z', nodes
                cmtynodes[c].update(nodes)
                #print cmtynodes[c]

        #from fitz import interactnow
        # Check that all have the same total size:
        if len(results) > 0:
            assert all(results[0].N == r.N  for r in results[1:])
        return results

class LouvainWeighted(Louvain):
    weighted = True

#class ModularitySA(CDMethod):
#    _input_format = 'edgelist'
#    binary = '/home/richard/research/cd/code-dl/good_modularity_SA/modularity_sampling_v1.0.0/anneal.py'
#    weighted = False
#    quiet = True
#    trials = 1
#    graphfileExtension = '.pairs'
#    _nodemapZeroIndexed = True
#
#    def run_sa(self):
#        args = ('python',
#                self.binary,
#                #'-o', # store local optima only
#                self.graphfile,
#                )
#        p = subprocess.Popen(args, stdout=subprocess.PIPE,
#                             stderr=subprocess.STDOUT)
#
#        stdout, stderr = p.communicate()
#        if not self.quiet: print stdout
#        if p.returncode != 0:
#            raise ValueError("SA didn't return propertly: %d"%p.returncode)
#
#        self.stdout = stdout
#        raise
#
#    def read_cmtys(self):
#        results = [ ]
#        level = None
#        for line in self.stdout.split('\n'):
#            #print 'line:', line
#            # Find 'level' starting lines
#            r = re.match(r'level ([0-9]+):', line)
#            if r:
#                level = int(r.group(1))
#                # Initialization
#                cmtynodes = collections.defaultdict(set)
#            # If not in a 'level' block, then break.
#            if level is None: continue
#            # End condition
#            if re.match(r'  end computation: ', line):
#                cmtys = pcd.cmty.Communities(cmtynodes,
#                                             nodes=set(self.g.nodes()),
#                                             name=self.__class__.__name__+'_level%d'%level)
#                cmtys.modularity = modularity
#                results.append(cmtys)
#                level = None
#                del cmtynodes, modularity, cmtys
#            # Modularity:
#            r = re.match(r'  modularity increased from [0-9.+-]+ to ([0-9.+-]+)', line)
#            if r:
#                modularity = float(r.group(1))
#            # Actual node parsing
#            r = re.match(r'([\d]+) ([\d]+)', line)
#            if r:
#                node, c = int(r.group(1)), int(r.group(2))
#
#                def find_real_node(node, level):
#                    #print 'x', repr(node), repr(level)
#                    if level == 0:
#                        #print 'b'
#                        return set((node,))
#                    assert node in results[level-1]._cmtynodes
#                    lower_nodes = results[level-1]._cmtynodes[node]
#                    #nodes = set()
#                    #for lower_node in lower_nodes:
#                    #    print 'a', lower_node
#                    #    nodes.update(find_real_node(lower_node, level-1))
#                    #print 'y', nodes, level
#                    #return nodes
#                    return lower_nodes
#                nodes = find_real_node(node, level)
#                #print 'z', nodes
#                cmtynodes[c].update(nodes)
#                #print cmtynodes[c]
#
#        #from fitz import interactnow
#        return results
#
#    def run(self):
#        self.run_sa()
#        results = self.read_cmtys()
#        self.results = results
#        self.cmtys = self.results[0]
class ModularitySA(CDMethod):
    _input_format = 'edgelist'
    #binary = 'lancichinetti_codes/clustering_programs_5_1/bin/modopt'
    binary = 'lancichinetti_modSA/modopt/modopt'
    _nodemapZeroIndexed = True
    trials = 5
    resolution = 1
    temp_init = 1e-6
    temp_step = .9999
    initial = None

    def run(self):
        binary = _get_file(self.binary)
        if self.initial is not None:
            # Write initial seed structure
            f = open('community.dat', 'w')
            nodecmtys = self.initial.nodecmtys_onetoone()
            cmty_map = self.initial.cmtyintmap()
            for node, cmty in nodecmtys.iteritems():
                print >> f, self.vmap[node], cmty_map[cmty]
        args = (binary,
                self.graphfile,
                str(self.randseed),  # seed
                str(self.resolution), # resolution parameter (lambda)
                str(self.trials),   # number of runs
                str(self.temp_init),    # initial temperature
                str(self.temp_step),   # temperature step
                )
        self.call_process(args)

    def read_cmtys(self):
        results = [ ]
        level = None
        outfile = self.basename+'.netpart'
        cmtynodes = collections.defaultdict(set)
        for i, line in enumerate(open(outfile)):
            if not line.strip(): continue
            cmtynodes[i] = set(self.vmap_inv[int(x)] for x in line.split())

        self.cmtys = pcd.cmty.Communities(cmtynodes)
        self.results = [ self.cmtys ]
        return self.results

class _PCD_single(CDMethod):
    _input_format = 'null'
    trials = 10
    gamma = 1
    const_q = False
    minargs_default = { }
    minargs = { }
    initial = None # initial state for each minimization attempt.  Default is 'random'.
                   # Same as Graph.trials(initial=)
    map_nodes_to_int = False
    def _initG(self, G):  # for subclassing
        pass
    def run(self):
        G = pcd.Graph.fromNetworkX(self.g)
        if self.verbosity <= 0:
            G.verbosity = -1
        if self.const_q:
            G.const_q = self.const_q
        kwargs = self.minargs_default.copy()
        kwargs.update(self.minargs)

        self._initG(G)
        initial = self.initial
        if isinstance(initial, pcd.cmty.Communities):
            initial = initial.to_pcd_cmtystate(G)
        if initial is None:
            initial = 'random'
        G.trials(gamma=self.gamma, trials=self.trials,
                 minimizer=self.minimizer, #kwargs are args to minimizer
                 initial=initial,
                 **kwargs
                 )
        self.G = G
        self.cmtys = pcd.cmty.Communities.from_pcd(G)
        self.results = [ self.cmtys ]
        if hasattr(self, '_finalize'):
            self._finalize(G)
class _PCD_SA(object):
    minimizer = 'anneal'
    # From Sa code:
    attempts = 1
    betafactor = 1.001
    Escale = None
    new_cmty_prob = .01
    p_binary = .05
    p_collective = .05
    mode = None
    # Constant q mode:
    const_q_SA = False
    constqEcoupling = 1
    # Minimum n mode:
    min_n_SA = False
    minnEcoupling = 1
    maxrounds = None
    _args_names = ('betafactor', 'attempts', 'Escale',
                   'mode',
                   'new_cmty_prob', 'p_binary', 'p_collective',
                   'const_q_SA', 'constqEcoupling',
                   'min_n_SA', 'minnEcoupling', 'maxrounds')
    @property
    def minargs_default(self):
        a = super(_PCD_SA, self).minargs_default.copy()
        a.update((name, getattr(self, name)) for name in self._args_names)
        return a
# PCD APM methods
class APMGreedy(_PCD_single):
    minimizer = 'greedy'
    trials = 10
class APMAnneal(_PCD_SA, _PCD_single):
    minimizer = 'anneal'
    trials = 10
# PCD modularity methods
class _PCDmod(object):
    gamma = 1  # required for modularity
    def _initG(self, G):
        G.enableModularity(mode=None)
    def _finalize(self, G):
        self.modularity = G.modularity()
class PCDmodGreedy(_PCDmod, _PCD_single):
    minimizer = 'greedy'
    trials = 10
class PCDmodSA(_PCD_SA, _PCDmod, _PCD_single):
    pass



class BeliefPropogationq(CDMethod):
    _input_format = 'null'
    q = None # must be specified (so far)
    initial = None   # Initial values for communities
    epsilon = None   # Initial c_out/c_in ratio.  c_in+c_out = avg degree of graph.
    p_in = None      # Seed learning c-matrix with this value on diagonal
    p_out = None     # See learning c-matrix with this value off-diagonal
    sizes = None     # Initial communiy sizes
    def run(self):
        #import sys
        #sys.path.append('/home/richard/cavitymethod')
        #import cavity
        #del sys.path[-1]
        import socket
        if 'becs' in socket.getfqdn():
            cavity = __import__('pcd.support.cavity_becs',
                                globals(), locals(), ['cavity'], -1).cavity
        else:
            from pcd.support import cavity
        if self.q is None:
            raise ValueError("q is a required argument (and must be integer).")
        for i in range(10):
            c = cavity.Cavity(self.g, self.q,
                              **dict((name, getattr(self, name)) for name in
                                     ('epsilon', 'p_in', 'p_out', 'initial', 'sizes')))
            c.verbosity = self.verbosity
            #from fitz import interactnow
            try:
                c.run_learn()
            except ValueError:
                print "NaN error occured.  Restarting.  Iteration %s"%i
                continue
            break
        self.cavity = c
        self.p = c.p
        self.n = c.n
        self.f = c.f()

        self.cmtys = pcd.cmty.Communities(c.cmtynodes(), set(self.g.nodes()))
        self.results = [ self.cmtys ]


class BeliefPropogationZ(CDMethod):
    graphfileextension = '.gml'
    _input_format = 'gml'
    full = False
    q = None  # number of communities - Must be specified.
    epsilon = None   # Initial c_out/c_in ratio.  c_in+c_out = avg degree of graph.
    p_in = None
    p_out = None
    binary = 'BP_z__sbm/sbm'
    def run(self):
        binary = _get_file(self.binary)
        args = [binary, 'learn', '-l', self.graphfile, '-q', str(self.q),
                '-w', self.graphfile+'.out',
                #'-W', self.graphfile+'.spm.out',
                '-M', self.graphfile+'.marginal.out',
                #'--wcab', self.graphfile+'.cab.out',
                #'--wna', self.graphfile+'.na.out',
                '--wall', self.graphfile,
                #'-M', self.graphfile+'.marginal.out',
                '-n', str(len(self.g)),
                '-p'+(','.join(("%f"%(1./self.q), )*self.q )),
                ]
        # -f: full_flag, default False.  if T, then non-sparse implementation.  full is O(n^2).
        # -c: set c_ab
        # -p: set p_a
        # -P<epsilon><avgdegree>
        if self.full: args.append('-f')
        #if self.verbosity<= 0:  args.extend(('-v', '-1'))
        args.extend(('-v', str(self.verbosity)))
        #if self.seed: args.extend(('-D', str(self.seed)))
        args.extend(('-D', str(self.randseed)))
        args.extend(self.get_params())

        self.call_process(args)
    def get_params(self):
        """Return command arguments related to initial parameters."""
        args = [ ]
        avg_degree = 2*self.g.number_of_edges()/float(len(self.g))
        if self.p_in is not None and self.p_out is not None:
            args.append('-P%s,%s'%(self.p_out/float(self.p_in), avg_degree))
        elif self.epsilon is not None:
            args.append('-P%s,%s'%(self.epsilon, avg_degree))
        return args

    def read_cmtys(self):
        # Read communities.
        # free energies and metadata
        vmap_inv = self.vmap_inv
        f = open(self.graphfile+'.marginal.out')
        line = f.readline()
        m = re.match('f=([-\d.]+|-?nan)\s+L=([-\d.]+|-?nan)\s+argmax_ovl=([-\d.]+|-?nan)', line)
        fe = float(m.group(1))
        L = float(m.group(1))
        argmax_ovl = float(m.group(1))
        # Read the communities.  The code implies that these should be
        # written out in the same order as in the gml file.  The gml
        # file should be written in the order of g.nodes(), which is
        # saved in self.node_order.  Let's hope that these invariants
        # are not changed.
        line = f.readline()
        assert line.startswith('argmax_configuration:')
        line = f.readline()
        cmtylist = [ int(c) for c in line.split() ]
        cmtynodes = dict()
        for node, c in zip(self.node_order, cmtylist):
            if c not in cmtynodes: cmtynodes[c] = set()
            #cmtynodes[c].add(node)
            cmtynodes[c].add(vmap_inv[node])
        cmtys = pcd.cmty.Communities(cmtynodes)
        # Marginals
        while 'marginals:' not in f.readline():
            pass
        marginals = [ ]
        for i, node in enumerate(self.node_order):
            marg = f.readline().split()
            marg = [ float(x) for x in marg ]
            marginals.append(marg)
        # c_ab matrix - scaled affinity matrix, p_{a,b}*N
        f = open(self.graphfile+'.cab')
        cab = f.readlines()
        cab = [ line.split() for line in cab if line.strip() ]
        cab = [ [float(x) for x in line] for line in cab ]
        cab = numpy.asarray(cab)
        self.c = cab
        # n_a matrix - group occupancies
        f = open(self.graphfile+'.na')
        na = f.readlines()
        na = [ float(line.strip()) for line in na if line.strip() ]
        na = numpy.asarray(na)
        self.n = na

        # Save all our data
        self.f = fe
        self.marginals = marginals
        self.cmtys = cmtys
        self.results = [ cmtys ]

        #from fitz import interactnow
        return self.results

class _2WalkSpectrum(CDMethod):
    _input_format = 'null'
    q_hint = None  # calculate q_hint+1 eigenvalues.
    q = None  # Specify the number of eigenvalues.
    kwargs = { }
    map_nodes_to_int = False
    def run(self):
        if self.q:
            # Add one to self.q, so that we have enough information to
            # check if q was correct.
            num_ev = self.q+1
        elif self.q_hint:
            num_ev = self.q_hint+1
        else:
            num_ev = None
        S = self.SpectralMethod(self.g, num_ev=num_ev, **self.kwargs)
        self.cmtys = S.cmtys(q=self.q)
        self.results = [ self.cmtys ]

        self.q = S.q()
        if hasattr(S, 'c'):
            self.c = S.c()
        self.ev = S.ev_sorted
        self.S = S
        self.time_ev = S.time_ev
class NonBacktrackingWalkSpectrum(_2WalkSpectrum):
    from pcd.support.spectral import NonBacktrackingWalk as SpectralMethod


class FlowSpectrum(_2WalkSpectrum):
    from pcd.support.spectral import Flow as SpectralMethod



if __name__ == "__main__":
    if len(sys.argv) < 2 or '-h' in sys.argv or '--help' in sys.argv:
        print "%s MethodName infile outfile"%os.path.basename(sys.argv[0])
        print "available methods:", ', '.join(sorted(
            name for name, cda in globals().iteritems()
            if isinstance(cda, type)
            and issubclass(cda, CDMethod)
            and name[0]!='_'))
        sys.exit()
    #import optparse
    method = sys.argv[1]
    input = sys.argv[2]
    output = sys.argv[3]
    options = sys.argv[4]
    from pcd.util import leval
    options = leval(options)
    print options

    import pcd.io
    g = pcd.io.read_any(input)

    options['basename'] = os.path.basename(input)

    method = locals()[method]
    r = method(g, **options)
    #initial = r.cmtys
    #r = method(g, initial=initial, **options)
    from fitz import interactnow


