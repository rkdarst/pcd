# Richard Darst, October 2012

import collections
import itertools
import math
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

tmpdir_base = '.'

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
        search_path.extend(os.environ['PCD_EXTERNAL_CODE_PATH'].split(':'))
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
    delete_tmpdir = None          # Delete the temporary directory after running?
    _nodemapZeroIndexed = False   # Should edgelist be zero-indexed?  default start at 1.
    _interact_results = False     # If True, after self.run() is called, interact in directory.
    seed = None                   # Seed to override.
    verbosity = 2
    _map_nodes_to_int = True      # Do we need to remap all nodes to ints?
    _run_method = True
    _load_results = True
    _recovery_mode = False        # recovery mode - do nothing, user should
                                  # call .run(), .read_cmtys() etc.
    weighted = False              # Write a weighted edgelist/turn on weights?

    def __init__(self, g, dir=None, basename=None, cache=None, **kwargs):
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
        if g is None:
            # This is used to force the method to not actually do
            # anything, but return an object from which the .name()
            # attribute can be accessed.  This may be used when the
            # class is wrapped in a partial object.
            return
        if cache:
            if isinstance(cache, str):
                if os.path.exists(cache):
                    data = pickle.load(open(cache, 'rb'))
                    self.results = data['results']
                    self.cmtys = data['cmtys']
                    return
            if len(cache) > 1 and hasattr(cache[0], '__getitem__'):
                key = cache[1]
                cache_d = cache[0]
                if key in cache_d:
                    self.results = cache_d[key]['results']
                    self.cmtys   = cache_d[key]['cmtys']
                    return
        #
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
                getattr(self, 'write_'+self._input_format)(self.g, self.graphfile, weighted=self.weighted)
            if self._run_method        and not self._recovery_mode:
                self.run()
            if self._load_results  and not self._recovery_mode:
                self.read_cmtys()
            if self._interact_results:
                import code; code.interact(banner=("In dir=%s\n"
                                   "File is %s")%(self.dir, self.graphfile),
                                   local=locals())
        if self.delete_tmpdir:
            shutil.rmtree(self.dir)
        # Write data to cache, if needed
        if cache:
            data = dict(cmtys=self.cmtys, results=self.results)
            if isinstance(cache, str):
                if os.path.exists(cache):
                    pickle.dump(data, open(cache, 'wb'))
            if len(cache) > 1 and hasattr(cache[0], '__getitem__'):
                key = cache[1]
                cache_d = cache[0]
                cache_d[key] = data

    # Utility methods
    @classmethod
    def name(cls):
        """Return CD method name, as a string"""
        return getattr(cls, '_name', cls.__name__)
    @property
    def _randseed(self):
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
        if not self._map_nodes_to_int or isinstance(self.g, GraphFile):
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
        if all(isinstance(x, int) for x in self.g.nodes_iter()):
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
                dir=tmpdir_base, )
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
        if hasattr(self, '_filter_graphfilename'):
            self.graphfile = self._filter_graphfilename(self.graphfile)
        if not os.access(self.dir, os.F_OK): os.mkdir(self.dir)
    def call_process(self, args, binary_name=None, stdout_suffix="",
                     env=None):
        """Call a process, saving the standard output to disk.

        This handles running a program, while also doing other
        important tasks like saving the times, command line arguments,
        return value, and stdout to disk.
        """
        # Environment update, if requested.
        env_new = None
        env_new_string = ''
        if env:
            env_new = dict(os.environ)    # make a copy
            env_new.update(env)
            for k,v in sorted(env.iteritems()):
                env_new_string += '%s=%s '%(k,v)


        if self.verbosity >= 2:
            print env_new_string + " ".join(args)
        # Stdout is saved to BINARY_NAME.stdout, BINARY_NAME defaults
        # to args[0]
        if binary_name is None:
            binary_name = os.path.basename(args[0])
        binary_name = binary_name.replace('/', '_')
        # Open a stream for the standard output
        self._stdout_fname = binary_name+'.stdout'+stdout_suffix
        disk_stdout = open(self._stdout_fname, 'w')
        # Write the command lines we are running
        disk_stdout.write("Command line: %s\n"%repr(list(args)))
        if env_new_string:
            disk_stdout.write("Environment variables: %s\n"%env_new_string)
        disk_stdout.write("Algorithm instance attributes: %s\n"%(
            pcd.util.listAttributes(self, exclude=['_randseed', 'vmap', 'vmap_inv'],
                                    exclude_deep=['initial', ]),))
        disk_stdout.write("Process times (os.times()): %s\n"%(os.times(),))
        disk_stdout.write("Beginning run at %s\n"%time.ctime())
        disk_stdout.write("==========\n\n\n")
        disk_stdout.flush()
        #time_start = time.time()

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
            ret = subprocess.call(args, stdout=stdout,
                                  stderr=subprocess.STDOUT, env=env_new)
        finally:
            #if self.verbosity >= 2:
            #    print open(binary_name+'.stdout').read()
            if self.verbosity >= 2:
                splitter.close()
                splitter.wait()

            #print "call_process time elapsed (s):", time.time() - time_start
            disk_stdout.write('\n\n\n==========\n')
            disk_stdout.write("Ending run at %s\n"%time.ctime())
            disk_stdout.write("Process times: %s\n"%(os.times(),))
            disk_stdout.write("Return value: %s\n"%ret)
            disk_stdout.flush()
        if ret != 0:
            raise ValueError("CD method %s (%s) command failed: %s."%(
                self.name(), self.graphfile, ret))
        return ret

    # Format writer methods
    def write_pajek(self, g, fname, weighted=False):
        f = open(fname, 'w')
        print >> f, '*Vertices %d'%len(g)
        vmap = self.vmap
        for i, n in enumerate(g.nodes_iter()):
            #if ' ' in repr(n): raise ValueError("Node '%s' has a space in its name"%n)
            print >> f, i, '%s'%vmap[n]#, "0.00", "0.00", "ellipse"
            #print >> f, i, repr(n)#, "0.00", "0.00", "ellipse"
        print >> f, "*Edges"
        for u, v, d in g.edges_iter(data=True):
            if weighted:
                print >> f, vmap[u], vmap[v], float(d.get('weight', 1.0))
            else:
                print >> f, vmap[u], vmap[v] #, "1.0"
        f.flush()
        f.close()

    def write_edgelist(self, g, fname, weighted=False):
        f = open(fname, 'w')
        vmap = self.vmap
        for u, v, d in g.edges_iter(data=True):
            #if getattr(self, 'weighted', True):
            if weighted:
                # weighted
                print >> f, vmap[u], vmap[v], float(d.get('weight', 1.0))
            else:
                # unweighted
                print >> f, vmap[u], vmap[v]
    def write_edgelistWeighted(self, g, fname, weighted=True):
        """Write an edgelist with forced weights."""
        self.write_edgelist(g, fname, weighted=True)
    def write_gml(self, g, fname):
        # networkx writes nodes in order of g.nodes()
        #self.node_order = g.nodes()
        #for i, n in enumerate(g.nodes_iter()):
        #    g.node[n]['id'] =
        #assert [ int(x) for x in g.nodes() ]
        g = networkx.relabel_nodes(g, mapping=self.vmap)  # Make a copy.
        networkx.write_gml(g, fname)
    def write_null(self, g, fname, weighted=False):
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
    def run_infomap(self):
        args = (#'/usr/bin/gdb', '--args',
                _get_file(self._binary),
                str(self._randseed),
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
        self.cmtys.label = "InfomapOld"
        self.results = [ self.cmtys ]
        return self.results

    def run(self):
        self.run_infomap()


class _InfomapOld(_InfomapOld_):
    _binary = 'infomap/infomap_undir/infomap'

class _InfomapOld_dir(_InfomapOld_):
    _binary = 'infomap/infomap_dir/infomap'
    _is_directed = True

class Infomap(CDMethod):
    """New infomap code.

    This runs Infomap code from http://mapequation.org/code.html.

    hierarchical: bool
        If true, use hierarchical infomap code.  If false, return just
        one partition of the network.  Note that this corresponds to
        the '--two-level' option to the Infomap binary.

    directed: bool
        If true, have Infomap treat the network as directed.

    full_lower_levels: bool
        If true(default), each level is a complete partition of the
        system.  If false, each lower level (smaller communitie)
        includes only the communities that are different from the
        higher levels.

    trials: int
        Number of minimization trials.

    _args: list of str
        List of strings to be added to the arguments.

    initial: Communities object
        Initial configuration for minimization.

    include_self_links: bool
        If true, Infomap will consider self-links in random walk dynamics.

    """
    threads = None
    _binary = 'infomap/Infomap-0.11.5/Infomap'
    _input_format = 'edgelist'
    _nodemapZeroIndexed = True
    directed = False
    hierarchical = True
    full_lower_levels = True
    initial = None
    trials = 10
    include_self_links = False
    _args = [ ]
    which_level = 0
    def run(self):
        args = [_get_file(self._binary),
                '--num-trials=%d'%self.trials,
                '--input-format=link-list', '--zero-based-numbering',
                '--seed=%d'%self._randseed,
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
        env = { }
        if self.threads is not None:
            env.update({'OMP_NUM_THREADS': str(int(self.threads))})
        args.extend(self._args)
        self.call_process(args, env=env)

    def read_cmtys(self):
        fname = self.basename+'.tree'
        self.results = self.read_infomap_cmtys(
            fname=fname, vmap_inv=self.vmap_inv,
            full_lower_levels=self.full_lower_levels)
        self.cmtys = self.results[self.which_level]
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

class InfomapHighest(Infomap):
    which_level = -1


class _Oslom(CDMethod):
    """Oslom.

    Code available at http://oslom.org/

    Lisence:

    References:

    merge_singletons: bool
        If true, return the partitions which have singletons (homeless
        nodes) merged into their best community.

    strip_singletons: bool
        If there are singletons in the result, do not return those in
        the detected structure (and thus the result will not be a
        cover).  This is a custom option, not from the core oslom
        code.  If the merge_singletons option is true, this option
        will do nothing.

    initial: Communities object
        Initial configuration.  Used with -hint option.  Cleans up and
        only returns statistically significant communities.  Should be
        used with trials=0 and trials_hier=0 to prevent a new search for
        communities, and these are not automatically set.

    weighted: bool
        Use -w option, which specifies floating point weights.  If
        False, integer weights are transformed into multiple edges and
        non-integer weights are an error.

    trials: int, optional
        Number of trials to run for first hierarchical level, '-r'
        option.  For fastest execution, set this option and
        trials_hier to 1.

    trials_hier: int, optional
        Number of trials for each higher level.  These levels have
        many fewer nodes, so the default is much higher than for
        `trials`.

    fast: int, optional
        If fast=1, specifies 'fast mode', trials=1 and trials_hier=1.
        This is the fastest possible run time.  If fast is an integer
        greater than one, then set both 'trials' and `trials_hier` to
        this many trials.  Thus, fast=2 can be used for 'fast, but not
        quiet as fast as fast=True'.  'fast' is incompatible with
        specifying 'trials' and 'trials_hier'.

    p_value: float, optional
        Significance level of modules, default to upstream default of
        0.1.

    coverage_parameter: float, optional
        Submodule merging criteria.  Default in Oslom code .5, range 0
        to 1.  Larger leads to bigger clusters.
    """
    #_binary = 'oslom/OSLOM2/oslom_undir'
    _input_format = 'edgelist'
    trials = 10
    trials_hier = 50
    fast = None

    merge_singletons = False
    initial = None
    weighted = None
    strip_singletons = False
    p_value = None
    coverage_parameter = None
    which_level = 0

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

        args = (_get_file(self._binary),
                "-seed", str(self._randseed),
                "-w" if self.weighted else '-uw', #unweighted or weighted
                "-f", self.graphfile,
                )
        # Specification of how many trials to run.
        # fast mode.  trials and trials_hier must not be specified.
        if self.fast:
            if ('trials' in self.__dict__
                or 'trials_hier' in self.__dict__):
                raise ValueError("`fast` incompatible with `trials` or `trials_hier`")
            if self.fast == 1:
                args = args + ('-fast', )
            else:
                args = args + ('-r', '%d'%self.fast, '-hr', '%d'%self.fast, )
        # Use trials and trials_hier
        else:
            args = args + (
                '-r', '%d'%self.trials,
                '-hr', '%d'%self.trials_hier,
                )
        #if not self.merge_singletons:
        if True:
            # The old code defaulted to singletons for homeless nodes,
            # and had an option "-all" to specify that they should be
            # merged.  In the new code, the situation is reversed, and
            # there is an option called "-singlet" which causes
            # homeless nodes to be included.  However, no matter what,
            # if you use -singlet, there is tpN_without_singletons
            # written.  So, we always run with -singlet, and we select
            # which files to read at the end.
            args = args + ('-singlet', )
        if self.initial:
            args = args + ('-hint', "initial_state.txt",)
            f = open("initial_state.txt", 'w')
            if self.trials != 0 or self.trials_hier != 0:
                raise ValueError("Oslom initial state requires trials=0 and trials_hier=0")
            for cmty, nodes in self.initial.iteritems():
                f.write(' '.join(str(self.vmap[n]) for n in nodes))
                f.write('\n')
            f.close()
        if self.p_value:
            args = args + ('-t', '%f'%self.p_value)
        if self.coverage_parameter:
            args = args + ('-cp', '%f'%self.coverage_parameter)

        self.call_process(args)

    def read_cmtys(self):
        self.outdir = self.graphfile + '_oslo_files/'
        self.results = [ ]
        for level in itertools.count():
            if self.merge_singletons:
                # Files with singletons merged
                fname = self.outdir+'tp%s_without_singletons'%(level if level else '')
            else:
                # Files with singletons in homeless communities
                fname = self.outdir+'tp%s'%(level if level else '')
            if not os.access(fname, os.F_OK):
                if self.verbosity >= 3:
                    print "*** Breaking at: level %s"%level
                break

            cmtys = self.read_oslom_cmtys(fname=fname)
            cmtys.label = "level%02d"%level

            self.results.append(cmtys)
            if not self.strip_singletons and isinstance(self.g, networkx.Graph):
                assert len(self.g) <= self.results[-1].cmtysizes_sum(), \
                       "Too few nodes detected - Graph disconnected?"
        self.cmtys = self.results[self.which_level]
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
            if not isinstance(self.g, networkx.Graph): raise NotImplemented
            cmtynodes[0] = set(self.g.nodes())
        # If we can strip singletons, then we can end up not spanning
        # the whole network.  In that case we must specify nodes= to
        # the communities object, or else it won't be able to
        # recalculate the universe correctly when asked.
        nodes = None
        if self.strip_singletons:
            if not isinstance(self.g, networkx.Graph): raise NotImplemented
            nodes = set(self.g.nodes())
        return pcd.cmty.Communities(cmtynodes,
                                    nodes=nodes,
                                    )

class Oslom(_Oslom):
    _binary = 'oslom/OSLOM2/oslom_undir'
class OslomWeighted(Oslom):
    weighted = True
class Oslom_dir(_Oslom):
    _binary = 'oslom/OSLOM2/oslom_dir'
    _is_directed = True
class OslomMerge(Oslom):
    """Oslom with merge_singletons=True.

    Singletons are merged into the communities they best connect to."""
    merge_singletons = True
class OslomHighest(Oslom):
    which_level = -1



class _Copra(CDMethod):
    """

    trials: int
        Number of trial runs.

    max_overlap: int
        Maximum number of communities per node.  With max_overlap=1,
        Copra returns a strict partition.

    weighted: bool
        Treat graph as weighted.

    max_overlap_range: list of int
        Not handled yet.

    """
    # http://www.cs.bris.ac.uk/~steve/networks/software/copra-guide.pdf
    _input_format = 'edgelist'
    trials = 10
    _binary = 'copra.jar'
    weighted = False
    max_overlap = 1   # default in COPRA code is 1.
    max_overlap_range = None

    def run(self):
        args = ["/usr/bin/java",
                "-cp", _get_file(self._binary),
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

    def read_cmtys(self):
        fname = 'clusters-'+os.path.basename(self.graphfile)
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
    _map_nodes_to_int = False

    def run(self):
        import pcd.old.multiresolution as multiresolution
        try:
            self.run_APM()
        except multiresolution.NoResultsError:
            self.cmtys = pcd.cmty.Communities(cmtynodes={0:set(self.g.nodes())})
            self.cmtys.label = "NullResults"
            self.results = [ self.cmtys ]
    def run_APM(self):
        import pcd.old.auto as auto

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
                              verbosity=self.verbosity,
                              )
        defaultOptions.update(self.options)
        options = defaultOptions
        options.update(dict(basename='graph',
                            plot=False))

        # At this point we are chdir'ed into the directory.
        a = auto.Auto(options=options)
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

    weighted: bool
        Runs in weigted mode (-w option).

    stop_dQ: float, optional
        If specified, stop when modularity change is this much.

    trials: int
        Do this many trial runs.  This is implemented by multiple runs
        of the Louvain binary, since it does not implement this
        itself.

    which_partition: int
        Louvain returns a hierarchy of partitions.  If which_partition
        is 0, maximize modularity of the the lowest level (smallest
        communities).  If which_partition is 1, maximize modularity of
        the level with second-smallest communities.  If
        which_partition is -1, maximize the modularity of the
        partition with largest communities (this is the level with the
        greatest modularity.)

        All partitions are stored in self.results.  One partition is
        stored in self.cmtys.  which_partition is either an index (0,
        1, ... -1) of the results list.  0 is the most granular
        (smallest communty size), and -1 is the least granular
        (largest community size).  which_partition can also be the
        string 'modmax', in which case the partition of maximum
        modularity is taken.  If we have multiple trials, return the
        set of partition with the maximum modularity of this level.


    At the end, the class has these attributes as information:

    num_levels: number of levels in the results returned.

    avg_num_levels: number of levels detected, average over all trials.

    modularity: modularity of the partition returned (self.cmtys)

    results_modularity: modularity of each partition in self.results.

    """
    _input_format = 'edgelist'
    _binary_convert = 'louvain/Community_latest/convert'
    _binary_community = 'louvain/Community_latest/community'
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
        # Do set up, creating initial binary file
        self.graphfilebin = self.graphfile+'.bin'
        #Remove old files.
        for fname in (self.graphfilebin, self.graphfilebin+'.weights'):
            if os.access(fname, os.F_OK):
                os.unlink(fname)
        args = (_get_file(self._binary_convert),
                '-i', self.graphfile, '-o', self.graphfilebin,
                )
        if self.weighted:
            args = args + ('-w', self.graphfilebin+'.weights', )
        self.call_process(args)
        # run louvain self.trials times.
        for i in range(self.trials):
            self.run_louvain(trial=i)

    def run_louvain(self, trial):
        """Run one louvian process"""
        args = (_get_file(self._binary_community),
                self.graphfilebin,
                '-l', '-1',
                '-v', )
        if self.weighted:
            args = args + ('-w', self.graphfilebin+'.weights', )
        if self.stop_dQ:
            args = args + ('-q', str(self.stop_dQ),)
        if self.initial:
            raise NotImplementedError("Check this before using.")
            f = open('initial.txt', 'w')
            for cname, cnodes in self.initial.iteritems():
                for node in cnodes:
                    print >> f, self.vmap[node], cname
            args = args + ('-p', 'initial.txt')
        self.call_process(args, stdout_suffix='.%03d'%trial)

    def read_cmtys(self):
        """Read all of the communities"""
        best_cmtys = None
        best_Q = -1e6
        self.avg_num_levels = []
        for i in range(self.trials):
            results = self.read_cmtys_and_return(trial=i)
            #print [x.modularity for x in results], [x.q for x in results]
            #print [ x.modularity-x.Q(self.g) for x in results ]
            assert len(results) > 0
            assert not any(x.modularity is None for x in results)
            self.avg_num_levels.append(len(results))
            # Pick which partition to use
            if self.which_partition == 'modmax':
                best_of_results = max(results, key=lambda x: x.modularity)
            else:
                best_of_results = results[self.which_partition]
            #print results
            if best_of_results.modularity > best_Q:
                self.results = results
                self.results_modularity = [r.modularity for r in results]
                self.cmtys = best_of_results #self.results[0]
                self.modularity = best_Q = best_of_results.modularity
            #print best_of_results.modularity, best_Q
        # Finalize
        self.num_levels = len(self.results)
        self.avg_num_levels = numpy.mean(self.avg_num_levels)

    def read_cmtys_and_return(self, trial):
        """Read one output file (of a specific trial) and parse."""
        fname = os.path.basename(self._binary_community)+'.stdout'+'.%03d'%trial
        stdout = open(fname).read()

        results = [ ]
        level = None
        modularity = None
        #print self.stdout
        #f_tmp = open('tmp.txt', 'w')
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
                if modularity is None:
                    raise ValueError("No modularity found for level")
                cmtys.modularity = modularity
                cmtys.label = 'level%02d'%level
                #level += 1
                results.append(cmtys)
                level = None
                del cmtynodes, modularity, cmtys
                modularity = None
            # Modularity:
            r = re.match(r'  modularity increased from [e0-9.+-]+ to ([e0-9.+-]+)', line)
            if r:
                modularity = float(r.group(1))
            # Actual node parsing
            r = re.match(r'([\d]+) ([\d]+)', line)
            if r:
                #print >> f_tmp, line
                node, c = int(r.group(1)), int(r.group(2))
                def find_real_node(node, level):
                    #print 'x', repr(node), repr(level)
                    if level == 0:
                        return set((self.vmap_inv[node],) )
                    assert node in results[level-1]._cmtynodes
                    lower_nodes = results[level-1]._cmtynodes[node]
                    # If you adjust this function for hierarchical
                    # community names, fix the duplicate-last-layer
                    # detection above.  In particular, .cmtysizes()
                    # will no longer be equal for the last two
                    # communities.
                    return lower_nodes
                nodes = find_real_node(node, level)
                #print 'z', nodes
                cmtynodes[c].update(nodes)
                #print cmtynodes[c]
        #f_tmp.flush()

        #from fitz import interactnow
        # Check that all have the same total size:
        #if len(results) > 0:
        #    assert all(results[0].N == r.N  for r in results[1:])
        # Is the last level duplicated?
        if results[-1].q == results[-2].q \
           and results[-1].cmtysizes() == results[-2].cmtysizes():
            # Our two communities should be equal.
            # A final test for equality (not used):
            #assert self.results[-1]._cmtynodes = self.results[-2]._cmtynodes
            del results[-1]
        else:
            # From my current understanding, the only time the two
            # final levels aren't identical is when stop_dQ is
            # nonzero.  So, insert a check for that, and if it turns
            # out to be not true, figure out what went wrong.
            assert self.stop_dQ

        ## Test the Louvain results against the hierarchy binary.
        #for layer, cmtynodes in enumerate(results):
        #    #print layer
        #    args = (_get_file('louvain/Community_latest/hierarchy'),
        #            #fname.rsplit('.', 1)[0],
        #            'tmp.txt',
        #            '-l', str(layer+1), )
        #    #print args
        #    self.call_process(args)
        #    _new_cmtynodes = dict()
        #    for line in open('hierarchy.stdout'):
        #        #print line
        #        try:
        #            line = line.split()
        #            node, cmty = int(line[0]), int(line[1])
        #            #print 'n,c:', node, cmty
        #            node = self.vmap_inv[node]
        #        except (ValueError, IndexError):
        #            continue
        #        _new_cmtynodes.setdefault(cmty, set()).add(node)
        #    cmtynodes = results[layer]._cmtynodes
        #    _new_cmtynodes = dict(_new_cmtynodes)
        #    assert cmtynodes == _new_cmtynodes
        #    if cmtynodes != _new_cmtynodes:
        #        print len(cmtynodes), len(_new_cmtynodes)
        #        #print sorted(cmtynodes.keys()), sorted(_new_cmtynodes.keys())
        #        print type(cmtynodes), type(_new_cmtynodes)
        #        for cname in _new_cmtynodes:
        #            if cmtynodes[cname] != _new_cmtynodes[cname]:
        #                print type(cmtynodes[cname]), type(_new_cmtynodes[cname])
        #                print cmtynodes[cname] in _new_cmtynodes.values(), \
        #                      _new_cmtynodes[cname] in cmtynodes.values()
        #                print cmtynodes[cname], _new_cmtynodes[cname]
        #                #print cname
        #                from fitz import interactnow
        #                #raise
        #        raise

        return results

class LouvainWeighted(Louvain):
    weighted = True

class LouvainLowest(Louvain):
    #weighted = True
    def read_cmtys(self):
        super(LouvainLowest, self).read_cmtys()
        self.results = [ self.cmtys ]

class LouvainModMax(Louvain):
    which_partition = 'modmax'
    #weighted = True
    def read_cmtys(self):
        super(LouvainModMax, self).read_cmtys()
        self.results = [ self.cmtys ]

# This does not currently work:
#class _LouvainOriginal(Louvain):
#    _binary_convert = 'louvain_original/Community_BGLL_CPP/convert'
#    _binary_community = 'louvain_original/Community_BGLL_CPP/community'
#    which_partition = 'modmax'

#class ModularitySA(CDMethod):
#    _input_format = 'edgelist'
#    _binary = '/home/richard/research/cd/code-dl/good_modularity_SA/modularity_sampling_v1.0.0/anneal.py'
#    weighted = False
#    quiet = True
#    trials = 1
#    graphfileExtension = '.pairs'
#    _nodemapZeroIndexed = True
#
#    def run_sa(self):
#        args = ('python',
#                self._binary,
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
    """

    trials: int
        Number of trials.

    resolution: float
        Resolution parameter.

    temp_init: float
        Initial simulated annealing temperature.

    temp_step: float
        Simulated annealing temperature step.

    initial: Communities object
        Initial starting configuration.

    """

    _input_format = 'edgelist'
    #_binary = 'lancichinetti_codes/clustering_programs_5_1/bin/modopt'
    _binary = 'lancichinetti_modSA/modopt/modopt'
    _nodemapZeroIndexed = True
    trials = 5
    resolution = 1
    temp_init = 1e-6
    temp_step = .9999
    initial = None

    def run(self):
        binary = _get_file(self._binary)
        if os.path.exists('community.dat'):
            # This file will be read regardless of if it is requested
            # or not.  Move it to another name if it alreday exists, I
            # think this is the simplest least unsafe way.
            os.rename('community.dat', 'community.dat.old')
        if self.initial is not None:
            # Write initial seed structure
            f = open('community.dat', 'w')
            nodecmtys = self.initial.nodecmtys_onetoone()
            cmty_map = self.initial.cmtyintmap()
            for node, cmty in nodecmtys.iteritems():
                print >> f, self.vmap[node], cmty_map[cmty]
        args = (binary,
                self.graphfile,
                str(self._randseed),  # seed
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
    _map_nodes_to_int = False
    def _initG(self, G):  # for subclassing
        pass
    def run(self):
        from pcd.old.models import Graph
        G = Graph.fromNetworkX(self.g)
        G.verbosity = self.verbosity
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
    """
    q: int
        Initial number of communities (Required so far)

    p_in: float, optional
        .
    """

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
    _binary = 'BP_z__sbm/sbm'
    def run(self):
        binary = _get_file(self._binary)
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
        args.extend(('-D', str(self._randseed)))
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
        self.node_order = [self.vmap[i] for i in range(len(self.g))]
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
    _map_nodes_to_int = False
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


class LinkCommunities(CDMethod):
    """Link communities

    Link communities in complex networks, Yong-Yeol Ahn, James
    P. Bagrow, Sune Lehmann

    http://barabasilab.neu.edu/projects/linkcommunities/
    Link communities reveal multiscale complexity in networks, Nature (doi:10.1038/nature09182)


    threshold: list of float
        Threshold is edge density cutoff for communities.  May be
        either a list or single number.  The first value is the
        'default' communities returned.
    """
    _input_format = 'edgelist'
    graphfileExtension = '.pairs'
    _nodemapZeroIndexed = True
    _binary_calc = 'linkcommunities/calcJaccards'
    _binary_cluster = 'linkcommunities/clusterJaccards'
    threshold = .5, .25, .75
    _threshold_doc = """\
    Threshold is edge density cutoff for communities.  May be either a
    list or single number.  The first value is the 'default'
    communities returned.
    """
    _threshold_type = 'list(float)'
    def run(self):
        """Run link community calculation.

        This only runs the first step and computes the .jaccs files.
        read_communities processes that file with the thresholds and
        calculates the actual communities.
        """
        jaccs_file = self.graphfile+'.jaccs'
        args = [_get_file(self._binary_calc),
                self.graphfile, jaccs_file,
                ]
        self.call_process(args)

    def _get_threshold(self, thr):
        """Compute communities at a threshold.

        Return filename.  If filename already exists, do not
        re-compute it."""
        jaccs_file = self.graphfile+'.jaccs'
        fname = self.graphfile+'.thr=%f'%thr
        # From the docs:
        # record all the clusters at THRESHOLD in net.clusters
        out_clusters = fname+'.clusters'
        # and the sizes of each cluster (number of edges and number of
        # induced nodes) to net.mc_nc.
        out_sizes = fname+'.mc_nc'

        if os.path.exists(out_clusters):
            return out_clusters
        # Run it
        args = [ _get_file(self._binary_cluster),
                 self.graphfile, jaccs_file, out_clusters, out_sizes,
                 str(thr),
                 ]
        self.call_process(args, binary_name=self._binary_calc+'.thr=%f'%thr)
        return out_clusters

    def read_cmtys(self):
        """Read and set all communities"""
        jaccs_file = self.graphfile+'.jaccs'
        self.results = [ ]
        # Ensure threshold is going to be an iterable
        threshold = self.threshold
        if isinstance(threshold, (int, float)):
            threshold = (threshold, )
        # For each threshold, do our calculation.
        for thr in threshold:
            # Calculate clusters at a threshold.  Do not re-calculate
            # it if the output files already exist on disk.
            fname = self._get_threshold(thr)
            # Read results in, and label it.
            cmtys = self.read_file(fname)
            cmtys.label = "Threshold=%f"%thr
            self.results.append(cmtys)
        self.cmtys = self.results[0]
        return self.results

    def read_file(self, fname):
        """Parse a single link community file and return Communities object."""
        cmtynodes = { }
        cname = 0
        for line in open(fname):
            if line.startswith('#'): continue
            line = line.split()
            if not line: continue
            if cname not in cmtynodes:
                cmtynodes[cname] = set()
            nodes = cmtynodes[cname]
            for link in line:
                a, b = link.split(',')
                a = int(a)
                b = int(b)
                nodes.add(self.vmap_inv[a])
                nodes.add(self.vmap_inv[b])
            cname += 1
        cmtys = pcd.cmty.Communities(cmtynodes)
        return cmtys


class Conclude(CDMethod):
    """CONCLUDE

    http://www.emilio.ferrara.name/conclude/ [Downloaded 2013-11-25]
    http://arxiv.org/abs/1303.1738

    CONCLUDE (COmplex Network CLUster DEtection) is a fast community
    detection algorithm.

    The strategy consists of three steps:
    - (re)weight edges by using a particular random walker;
    - calculate the distance between each pair of connected nodes;
    - partition the network into communities so to optimize the
      weighted network modularity.

    CONCLUDE is computationally efficient since its cost is near
    linear with respect to the number of edges in the network.

    The adoption of this community detection method has been proved
    worthy in different contexts, such as for studying online social
    networks and biological networks.

    References:

    Mixing local and global information for community detection in
    large networks.  P De Meo, E Ferrara, G Fiumara and A Provetti.
    Journal of Computer and System Sciences, 80(1):72-87 (2014).

    CONCLUDE is an improved version of the algorithm presented in the
    following paper, you might cite too:

    Generalized Louvain method for community detection in large
    networks.  P De Meo, E Ferrara, G Fiumara, and A Provetti.  ISDA
    '11: Proceedings of the 11th International Conference on
    Intelligent Systems Design and Applications, 2011.
    """
    _input_format = 'edgelist'
    _nodemapZeroIndexed = True
    _binary = 'conclude/CONCLUDE.jar'
    _mem_limit = None
    _mem_limit_doc = """\
    int (GiB), Set memory limit for the process.  Not needed necessarily.
    """
    def run(self):
        output = self.graphfile+'.out'
        args = ['java',
                #'-Xmx4G',  # suggested from the docs
                '-jar', _get_file(self._binary),
                self.graphfile, output,
                " ",  # delimiter (default \t)
                # 1=only compute weights, don't cluster.
                ]
        if self._mem_limit:
            # -Xmx4G
            args.insert(1, '-Xmx%dG'%self._mem_limit)
        self.call_process(args)

    def read_cmtys(self):
        output = self.graphfile+'.out'
        modularity = None
        cmtynodes = { }
        for line in open(output):
            if line.startswith('#'): continue
            if line.startswith('Q = '):
                modularity = float(line[4:])
                continue
            line = line.split()
            cmtynodes[len(cmtynodes)] = set(self.vmap_inv[int(x)] for x in line)
        self.cmtys = pcd.cmty.Communities(cmtynodes)
        self.cmtys.label = "Conclude"
        self.results = [ self.cmtys ]


class _SNAPmethod(CDMethod):
    _input_format = 'edgelist'
    graphfileExtension = '.txt'
    _nodemapZeroIndexed = True
    def read_cmtys(self):
        cmtys = pcd.cmty.Communities.from_clustersfile(
            'cmtyvv.txt',
            converter=lambda x: self.vmap_inv[int(x)],
            )
        #if hasattr(self, 'g'):   # This is NOT true.
        #    assert cmtys.N == len(self.g)
        cmtys.label = self.name()
        self.cmtys = cmtys
        self.results = [ self.cmtys ]

class SnapBigClam(_SNAPmethod):
    """BigClam via SNAP library.

    https://snap.stanford.edu/snap/index.html
    http://dl.acm.org/citation.cfm?id=2433471

    Formulates community detection problems into non-negative matrix
    factorization and discovers community membership factors of nodes.

    From the abstract:

    Network communities represent basic structures for understanding
    the organization of real-world networks. A community (also
    referred to as a module or a cluster) is typically thought of as a
    group of nodes with more connections amongst its members than
    between its members and the remainder of the network. Communities
    in networks also overlap as nodes belong to multiple clusters at
    once. Due to the difficulties in evaluating the detected
    communities and the lack of scalable algorithms, the task of
    overlapping community detection in large networks largely remains
    an open problem.

    In this paper we present BIGCLAM (Cluster Affiliation Model for
    Big Networks), an overlapping community detection method that
    scales to large networks of millions of nodes and edges. We build
    on a novel observation that overlaps between communities are
    densely connected. This is in sharp contrast with present
    community detection methods which implicitly assume that overlaps
    between communities are sparsely connected and thus cannot
    properly extract overlapping communities in networks. In this
    paper, we develop a model-based community detection algorithm that
    can detect densely overlapping, hierarchically nested as well as
    non-overlapping communities in massive networks. We evaluate our
    algorithm on 6 large social, collaboration and information
    networks with ground-truth community information. Experiments show
    state of the art performance both in terms of the quality of
    detected communities as well as in speed and scalability of our
    algorithm.
    """
    #_binary = 'snap/Snap/examples/bigclam/bigclam'
    _binary = 'snap/bigclam/bigclam'
    # Note: this binary had to be modified to support edgelists
    # separated by spaces.
    _threads = None
    q = -1   # Auto-detect q
    _q_doc = """Number of communities to detect (-1 autodetect, but within the given range)"""
    q_max = None  # default 100
    q_min = None  # default 5
    q_trials = None
    _q_trials_doc = """How many trials for the number of communities"""
    alpha = None
    _alpha_doc = "Alpha for backtracking line search"
    beta = None
    _beta_doc = "Beta for backtracking line search"
    def _filter_graphfilename(self, name):
        # This method special-cases reading of files contains
        # '.ungraph', so stript that from any input filenames.
        return name.replace('.ungraph', '')
    def run(self):
        #output = self.graphfile+'.out..'
        args = [_get_file(self._binary),
                '-i:%s'%self.graphfile,
                #'-o:%s'%output,
                #'-l:%s'%label_file,
                ]
        if self._threads:          args.append('-nt:%s'%str(self._threads))
        if self.q is not None:     args.append('-c:%d'%self.q)
        if self.q_min:             args.append('-mc:%d'%self.q_max)
        if self.q_max:             args.append('-xc:%d'%self.q_min)
        if self.q_trials:          args.append('-nc:%d'%self.q_trials)
        if self.alpha is not None: args.append('-sa:%f'%self.alpha)
        if self.beta  is not None: args.append('-sb:%f'%self.beta)
        self.call_process(args)

class SnapAGMfit(_SNAPmethod):
    """Fitting to an AGM, via SNAP library.

    https://snap.stanford.edu/snap/index.html
    Reference: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6413734&tag=1 (???)

    """
    _binary = 'snap/agmfit/agmfitmain'
    # Note: this binary had to be modified to support edgelists
    # separated by spaces.
    q = -1   # Auto-detect q
    _q_doc = """Number of communities to detect (-1 autodetect, but within the given range)"""
    epsilon = None
    _epsilon_doc = """Edge probability between the nodes that do not share any community (default (0.0): set it to be 1 / N^2"""
    def run(self):
        #output = self.graphfile+'.out..'
        #label_file = self.graphfile+'.labels'
        args = [_get_file(self._binary),
                '-i:%s'%self.graphfile,
                #'-o:%s'%output,
                '-l:',
                '-s:%d'%int(self._randseed),
                ]
        if self.q is not None:     args.append('-c:%d'%self.q)
        if self.epsilon is not None: args.append('-e:%f'%self.epsilon)
        self.call_process(args)


class _SNAPcommunity(_SNAPmethod):
    """Base class for the 'community' examples in SNAP.

    GN and CNM algorithms derived from this."""
    _binary = 'snap/Snap/examples/community/community'
    def run(self):
        args = [_get_file(self._binary),
                '-i:%s'%self.graphfile,
                #'-o:%s'%output,
                '-a:%d'%self._algorithm,
                ]
        self.call_process(args)
    def read_cmtys(self):
        f = open('communities.txt')
        cmtynodes = { }
        for line in f:
            if line.startswith('# Communities:'):
                # It tells us how many communities were found, so we
                # can pre-create all the community sets.
                Ncmty = int(line.split(':')[1])
                for cid in range(Ncmty):
                    cmtynodes[cid] = set()
            if line.startswith('#'): continue
            if not line.strip(): continue
            a, b = line.split()
            nid = int(a)
            cid = int(b)
            cmtynodes[cid].add(self.vmap_inv[nid])
        cmtys = pcd.cmty.Communities(cmtynodes)
        cmtys.label = self._algorithm_name
        self.cmtys = cmtys
        self.results = [ self.cmtys ]
class SnapGN(_SNAPcommunity):
    """Girvan-Newman algorithm in SNAP library.

    http://dx.doi.org/10.1073/pnas.122653799"""
    _algorithm = 1
    _algorithm_name = 'Girvan-Newman'
class SnapCNM(_SNAPcommunity):
    """Clauset-Newman-Moore algorithm in SNAP library.

    http://arxiv.org/abs/cond-mat/0408187"""
    _algorithm = 2
    _algorithm_name = 'Clauset-Newman-Moore'


class CliquePerc(CDMethod):
    """Clique percolation

    Source code from https://github.com/aaronmcdaid/MaximalCliques

    There is one parameter, k.  By default it is scanned from 3 to the
    highest value possible, but

    Parameters:

    min_k: int
        minimum clique size

    max_k: int, optional
        maximum clique size, or None or zero, in which case the method
        will self-terminate.

    k: int, optional
        As a shortcut, if k is set, this min_k and max_k are both set
        to this value.  Only this clique size is used for percolation
        then.

    Returns one layer per clique size.  The default layer is the
    lowest, or the min_k.
    """
    _binary = 'cliquepercolation/MaximalCliques/cp5'
    _input_format = 'edgelist'
    _nodemapZeroIndexed = True
    k = None
    min_k = 3
    max_k = None
    def run(self):
        args = [_get_file(self._binary),
                self.graphfile,
                '--comments',
                ".",   # output directory
                ]
        if self.k:
            args.extend(['-k', str(self.k), '-K', str(self.k)])
        else:
            if self.min_k != 3:  args.extend(['-k', str(self.min_k)])
            if self.max_k:  args.extend(['-K', str(self.max_k)])
        self.call_process(args)

    def read_cmtys(self):
        self.results = results = [ ]
        i_values = [ ]
        for i in itertools.count(self.min_k):
            fname = 'comm%d'%i
            if not os.path.exists(fname):
                break
            unmap = lambda x: self.vmap_inv[int(x)]
            cmtys = pcd.cmty.Communities.from_clustersfile(fname, converter=unmap)
            results.append(cmtys)
            i_values.append(i)
        # Set the labels.  k%d, with enough leading zeros so they sort properly.
        max_i = i
        label_format = "k%%0%dd"%math.ceil(math.log(max_i)/math.log(10))
        for i, cmtys in zip(i_values, results):
            cmtys.label = label_format%i
        # Set default
        self.cmtys = results[0]
        return results

class SeqCliquePerc(CDMethod):
    """Sequential clique percolation

    Kumpula, Jussi M., et al. 'Sequential algorithm for fast clique
    percolation.' Physical Review E 78.2 (2008): 026109.

    Source code from https://git.becs.aalto.fi/complex-networks/scp

    k: int
        .

    """
    # This method has different limits: self-loops must be removed,
    # and graph must not have multiedges.
    _binary = 'sequentialcliquepercolation/scp/k_clique'
    _input_format = 'edgelistWeighted'
    _nodemapZeroIndexed = True
    k = 4
    def run(self):
        args = [_get_file(self._binary),
                self.graphfile,
                '-o=output%03d.txt'%self.k,
                '-v',
                ]
        if self.k:
            args.extend(['-k=%d'%self.k])
        self.call_process(args)

    def read_cmtys(self):
        self.results = results = [ ]
        k_values = [ ]
        for k in [self.k]:
            fname = 'output%03d.txt'%k
            unmap = lambda x: self.vmap_inv[int(x)]
            cmtys = pcd.cmty.Communities.from_clustersfile(fname, converter=unmap)
            results.append(cmtys)
            k_values.append(k)
        # Set the labels.  k%d, with enough leading zeros so they sort properly.
        max_k = k
        label_format = "k%%0%dd"%math.ceil(math.log(max_k)/math.log(10))
        for k, cmtys in zip(k_values, results):
            cmtys.label = label_format%k
        # Set default
        self.cmtys = results[0]
        return results


class Ganxis(CDMethod):
    """GANXiS (formerly SLPA) method

    Source code from https://sites.google.com/site/communitydetectionslpa/ and
    https://sites.google.com/site/communitydetectionslpa/GANXiS_v3.0.2____.zip

    J. Xie and B. Szymanski, 'Towards Linear Time Overlapping
    Community Detection in Social Networks', PAKDD
    2012.  http://arxiv.org/abs/1202.2465/

    Parameters:

    directed: bool, default False
        True indicates that all links should be symmetric.  If you use
        an undirected network, this must be true to set the '-Sym 1'
        option, otherwise you might have un-recipriocated edges.

    weighted: bool, default True
        If true, use weights if present, and default to one.  If
        False, strip all weights so that all weights are equal to one.

    overlaps: bool
        .

    threshold: float
        .

    maxiter: int
        Maximum number of iterations.  Default in code 100.

    trials:  int
        Number of trials

    loopfactor: float
        .

    There are various other currently undocumented parameters.
    """
    _binary = 'GANXiS/GANXiS_v3.0.2/GANXiSw.jar'
    _input_format = 'edgelist'
    _nodemapZeroIndexed = True
    directed = False
    weighted = True
    overlaps = True
    threshold = None
    maxiter = None # default 100
    trials = 1
    loopfactor = 1.0
    # Other options: maxC, minC, ev (speaker and listerner rule of embedding)

    def run(self):
        args = ['java', '-jar', _get_file(self._binary),
                '-i', self.graphfile,
                '-seed', str(self._randseed),
                #'-d', ".",   # output directory
                #'-L', '1',  # Use only largest connected component.
                ]
        if not self.directed:
            args.extend(('-Sym', '1'))
        if not self.weighted:
            args.extend(('-W', '0'))
        if not self.overlaps:
            args.extend(('-ov', '0'))
        if self.threshold:
            args.extend(('-r', str(self.threshold)))
        if self.trials:
            args.extend(('-run', str(self.trials)))
        if self.maxiter:
            args.extend(('-t', str(self.maxiter)))
        if self.loopfactor != 1.0:
            args.extend(('-loopfactor', str(self.loopfactor)))

        self.call_process(args)

    def read_cmtys(self):
        results = [ ]
        sort_keys = [ ]
        for fname in sorted(os.listdir('output/')):
            m = re.search(r'_run[0-9]+_r([0-9.]+)_v([0-9]+)_T([0-9]+).icpm$', fname)
            thresh = float(m.group(1))
            unmap = lambda x: self.vmap_inv[int(x)]
            cmtys = pcd.cmty.Communities.from_clustersfile('output/'+fname, converter=unmap)
            cmtys.label = 'r%3.3f'%thresh
            sort_keys.append(thresh)
            results.append(cmtys)
        # Set the labels.  k%d, with enough leading zeros so they sort properly.
        sort_keys, results = zip(*sorted(zip(sort_keys, results)))
        # Set values
        self.cmtys = results[0]
        self.results = results
        return results

class Demon(CDMethod):
    """DEMON: a local-first discovery method for overlapping communities

    http://www.michelecoscia.com/?page_id=42
    http://dl.acm.org/citation.cfm?id=2339630

    This uses a parameter of epsilon=0.25 and does not use any
    weights.
    """
    _binary = 'demon/launch.py'
    _input_format = 'edgelist'
    _nodemapZeroIndexed = True

    def run(self):
        args = ['python', _get_file(self._binary),
                self.graphfile,
                ]
        self.call_process(args)

    def read_cmtys(self):
        fname = 'communities'
        cmtynodes = { }
        for line in open(fname):
            cid, cnodes = line.split('\t')
            cnodes = set(self.vmap_inv[int(x)] for x in cnodes.split(','))
            cmtynodes[int(cid)] = cnodes
        cmtys = pcd.cmty.Communities(cmtynodes)
        cmtys.label = 'Demon'
        self.results = [ cmtys ]
        self.cmtys = cmtys
        return self.results

class GreedyCliqueExp(CDMethod):
    """Greedy Clique Expansion

    Lee, Conrad, et al. 'Detecting highly overlapping community
    structure by greedy clique expansion.' arXiv preprint
    arXiv:1002.1827 (2010).
    Conrad Lee, Fergal Reid, Aaron McDaid, Neil Hurley
    http://arxiv.org/abs/1002.1827

    https://sites.google.com/site/greedycliqueexpansion/
    https://sites.google.com/site/greedycliqueexpansion/benchmarking---datasets-and-tools/GCEimpl_r2011-11-06.tgz

    Parameters (taken from the command help):

    minsize: int, default 4
        The minimum size of cliques to use as seeds. Recommend 4 as
        default, unless particularly small communities are required
        (in which case use 3).

    eta: float, default 0.6
        The minimum value for one seed to overlap with another seed
        before it is considered sufficiently overlapping to be
        discarded (eta). 1 is complete overlap. However smaller values
        may be used to prune seeds more aggressively. A value of 0.6
        is recommended.

    alpha: float, default 1.0
        The alpha value to use in the fitness function greedily
        expanding the seeds. 1.0 is recommended default. Values
        between .8 and 1.5 may be useful. As the density of edges
        increases, alpha may need to be increased to stop communities
        expanding to engulf the whole graph. If this occurs, a warning
        message advising that a higher value of alpha be used, will be
        printed.

    phi: float, default 0.75
        The proportion of nodes (phi) within a core clique that must
        have already been covered by other cliques, for the clique to
        be 'sufficiently covered' in the Clique Coveage Heuristic

    """
    _binary = 'greedycliqueexpansion/GCECommunityFinder/build/GCECommunityFinder'
    _input_format = 'edgelist'
    _nodemapZeroIndexed = True
    # Other options: maxC, minC, ev (speaker and listerner rule of embedding)
    minsize = 4
    eta = 0.6
    alpha = 1.0
    phi = 0.75

    def run(self):
        args = [_get_file(self._binary),
                self.graphfile,
                # Other parameters, all must be specified:
                # clique minsize 4
                # eta, overlap for joining ()  0.6
                # alpha, used in expanding seeds. 1.0
                # phi, propoortion of nodes already covered .75
                str(self.minsize),
                str(self.eta),
                str(self.alpha),
                str(self.phi),
                ]
        self.call_process(args)

    def read_cmtys(self):
        results = [ ]
        sort_keys = [ ]
        fname = 'GCECommunityFinder.stdout'
        found_cmtys = False
        cname = 0
        cmtynodes = { }
        for line in open(fname):
            if line == 'Finished\n':
                found_cmtys = True
                continue
            if not found_cmtys:
                continue
            if found_cmtys and len(line) == 1:
                # Done reading cmtys
                break
            cnodes = set(self.vmap_inv[int(x)] for x in line.split())
            cmtynodes[cname] = cnodes
            cname += 1
            continue

        cmtys = pcd.cmty.Communities(cmtynodes)
        cmtys.label = 'GCE'
        self.results = [ cmtys ]
        self.cmtys = self.results[0]
        return self.results


class _IGraphAlgorithm(CDMethod):
    """Community detection algorithms from the igraph library.

    Igraph is a Python library for handling complex networks.  It
    contains a variety of community detection methods under the GPLv2
    license.

    Documentation on the community detection methods in igraph is
    available starting here:
    http://igraph.org/python/doc/igraph.Graph-class.html#community_fastgreedy


    References:

    Csardi, Gabor, and Tamas Nepusz. 'The igraph software package for
    complex network research.' InterJournal, Complex Systems 1695.5
    (2006).
    http://www.necsi.edu/events/iccs6/papers/c1602a3c126ba822d0bc4293371c.pdf
    """
    _input_format = 'null'
    _processor = 'VertexClustering'
    _argmap = { }
    _nodemapZeroIndexed = True
    _kwargs = { }
    def run(self):
        if self.weighted and not self._weightable:
            NotImplementedError("This method can not use weights.")
        vmap = self.vmap

        # Create our keyword argument dictiory.  Sources:
        # self._kwargs, then using self._argmap to get attribute
        # options.
        kwargs = { }
        kwargs.update(self._kwargs)
        for attrname, optname in self._argmap.iteritems():
            if attrname == 'weights':
                continue
            kwargs[optname] = getattr(self, attrname)

        # Create actual igraph object.
        import igraph
        ig = igraph.Graph(edges=[(vmap[a], vmap[b]) for a,b in self.g.edges_iter()])
        try:
            method = getattr(ig, self._method)
        except AttributeError:
            raise RuntimeError("Your version of igraph does not have the Graph.%s method"%
                               self._method)
        resultobject = method(**kwargs)

        # do we do weights?
        if self.weighted and self._weightable:
            print "Using weighted: weighted"
            ig.es['weight'] = 1.0  # default weight
            for a, b, data in self.g.edges_iter(data=True):
                ig[vmap[a], vmap[b]] = data.get('weight', 1.0)
            weight_argument = self._argmap.get('weights', 'weights')
            kwargs[weight_argument] = 'weight'
            #from fitz import interactnow
            print ig.es['weight']

        processor = getattr(self, '_process_'+self._processor)
        results = processor(resultobject)

        self.results = results
        self.cmtys = results[0]

    def _process_VertexClustering(self, clu):
        vmap_inv = self.vmap_inv
        cmtys = { }
        for idx in range(len(clu)):
            cmtys[idx] = set()
        for i, c in enumerate(clu.membership):
            cmtys[c].add(vmap_inv[i])

        cmtys = pcd.cmty.Communities(cmtys)
        return [ cmtys ]

    def read_cmtys(self):
        return self.results

class _IGraphAlgorithmDendogram(_IGraphAlgorithm):
    dendogram_cut = None
    _processor = 'dendogram'
    def _process_dendogram(self, dend):
        if self.dendogram_cut is None:
            # Return either the upstream hint or the max modularity.
            iclusters = dend.as_clustering()
            results = self._process_VertexClustering(iclusters)
            return results
        if self.dendogram_cut == 'modmax':
            cutat = dend.optimal_count
            iclusters = dend.as_clustering(n=cutat)
            results = self._process_VertexClustering(iclusters)
            return results
        if isinstance(self.dendogram_cut, (list, tuple)):
            raise NotImplementedError("We don't support other dendogram cuts yet.")


class IgraphCNM(_IGraphAlgorithmDendogram):
    """Fast modularity method by Clauset-Newman-Moore


    This algorithm merges individual nodes into communities in a way
    that greedily maximizes the modularity score of the graph. It can
    be proven that if no merge can increase the current modularity
    score, the algorithm can be stopped since no further increase can
    be achieved.

    This algorithm is said to run almost in linear time on sparse
    graphs.

    This method returns the communities found by cutting the dendogram
    at the maximum modularity level.


    Reference: Reference: A. Clauset, M.E.J. Newman and C. Moore: Finding
    community structure in very large networks. Phys. Rev. E 70, 066111
    (2004).  http://arxiv.org/abs/cond-mat/0408187
    """
    _method = 'community_fastgreedy'
    _weightable = True
class IgraphInfomap(_IGraphAlgorithm):
    """Infomap method as included in the igraph library.


    References:

    Rosvall, Martin, and Carl T. Bergstrom. Maps of random walks on
    complex networks reveal community structure. Proceedings of the
    National Academy of Sciences 105.4 (2008): 1118-1123.
    http://www.pnas.org/content/105/4/1118.full

    M. Rosvall, D. Axelsson, and C. T. Bergstrom. 'The map
    equation.' The European Physical Journal-Special Topics 178.1
    (2009): 13-23.
    http://arxiv.org/abs/0906.1405
    """
    #vertex_weights  - not used
    _weightable = True
    trials = 10
    _method = 'community_infomap'
    _argmap = dict(trials='trials', weights='edge_weights')
class IgraphModEVNaive(_IGraphAlgorithm):
    """Naive leading Eigenvector of theb Modularity matrix.

    A naive implementation of Newman's eigenvector community structure
    detection. This function splits the network into two components
    according to the leading eigenvector of the modularity matrix and
    then recursively takes the given number of steps by splitting the
    communities as individual networks. This is not the correct way,
    however, see the reference for explanation. Consider using the
    correct community_leading_eigenvector method instead.


    q: int, optional
        Number of clusters to detect.  If not given, continue
        splitting until all leading eigenvector signs are the same.


    References:

    M.E.J Newman. 'Finding community structure in networks using
    the eigenvectors of matrices.' Phys. Rev. E 74.3 (2006):
    036104.
    http://arxiv.org/abs/physics/0605087

    """
    # This method has the ability to return a dendogram, but is not
    # enabled in pcd yet.
    _method = 'community_leading_eigenvector_naive'
    q = None
    _weightable = False
    _argmap = dict(q='clusters')
class IgraphModEV(_IGraphAlgorithm):
    """Leading Eigenvector of the Modularity matrix.

    Newman's leading eigenvector method for detecting community
    structure. This is the proper implementation of the recursive,
    divisive algorithm: each split is done by maximizing the
    modularity regarding the original network.


    q: int, optional
        Number of clusters to detect.  If not given, continue
        splitting until all leading eigenvector signs are the same.


    References:

    M.E.J Newman. 'Finding community structure in networks using
    the eigenvectors of matrices.' Phys. Rev. E 74.3 (2006):
    036104.
    http://arxiv.org/abs/physics/0605087
    """
    _method = 'community_leading_eigenvector'
    q = None
    _weightable = True
    _argmap = dict(q='clusters')
class IgraphLabelPropagation(_IGraphAlgorithm):
    """RAK label propagation

    Finds the community structure of the graph according to the label
    propagation method of Raghavan et al. Initially, each vertex is
    assigned a different label. After that, each vertex chooses the
    dominant label in its neighbourhood in each iteration. Ties are
    broken randomly and the order in which the vertices are updated is
    randomized before every iteration. The algorithm ends when
    vertices reach a consensus. Note that since ties are broken
    randomly, there is no guarantee that the algorithm returns the
    same community structure after each run. In fact, they frequently
    differ. See the paper of Raghavan et al on how to come up with an
    aggregated community structure.


    initial: communities object
        Initial configuration upon starting (not implemented yet)

    fixed: iterable
        Nodes to hold fixed (not implemented yet)


    References:

    Raghavan, Usha Nandini, Reka Albert, and Soundar Kumara. 'Near
    linear time algorithm to detect community structures in
    large-scale networks.' Physical Review E 76.3 (2007): 036106.
    http://arxiv.org/abs/0709.2938
    """
    #initial = None
    #fixed = None
    _method = 'community_label_propagation'
    _weightable = True
class IgraphLouvain(_IGraphAlgorithm):
    """Louvain method (fast unfolding)

    This is a bottom-up algorithm: initially every vertex belongs to a
    separate community, and vertices are moved between communities
    iteratively in a way that maximizes the vertices' local
    contribution to the overall modularity score. When a consensus is
    reached (i.e. no single move would increase the modularity score),
    every community in the original graph is shrank to a single vertex
    (while keeping the total weight of the adjacent edges) and the
    process continues on the next level. The algorithm stops when it
    is not possible to increase the modularity any more after
    shrinking the communities to vertices.


    References:

    Vincent D. Blondel,Jean-Loup Guillaume, Renaud Lambiotte, and
    Etienne Lefebvre. 'Fast unfolding of communities in large
    networks.' J. Stat Mech.: Theory and Experiment 2008.10 (2008):
    P10008.
    http://arxiv.org/abs/0803.0476
    """
    # list of VertexClusterings
    _method = 'community_multilevel'
    _weightable = True
    _kwargs = dict(return_levels=True)
    _processor = 'VertexClusteringList'
    def _process_VertexClusteringList(self, resultobj):
        results = [ ]
        for clu in resultobj:
            r = self._process_VertexClustering(clu)
            results.extend(r)
        return results
class IgraphModularity(_IGraphAlgorithm):
    """Optimal modularity via linear programming.

    Calculates the optimal modularity score of the graph and the
    corresponding community structure.

    This function uses the GNU Linear Programming Kit to solve a large
    integer optimization problem in order to find the optimal
    modularity score and the corresponding community structure,
    therefore it is unlikely to work for graphs larger than a few
    (less than a hundred) vertices. Consider using one of the
    heuristic approaches instead if you have such a large graph.


    References:
    """
    # membership vector
    _method = 'community_optimal_modularity'
    _processor = 'modlist'
    _weightable = True
    def _process_modlist(self, obj):
        """Process results format of the membership list
        """
        vmap_inv = self.vmap_inv
        import igraph.clustering
        if isinstance(obj, igraph.clustering.VertexClustering):
            memvec = obj.membership
            mod = None
        else:
            memvec, mod = obj
        cmtys = { }
        for i, c in enumerate(memvec):
            if c not in cmtys: cmtys[c] = set()
            cmtys[c].add(vmap_inv[i])
        cmtys = pcd.cmty.Communities(cmtys)
        cmtys.modularity = mod
        return [ cmtys ]
class IgraphGN(_IGraphAlgorithmDendogram):
    """Edge betweenness dendogram.

    Community structure based on the betweenness of the edges in the
    network.

    The idea is that the betweenness of the edges connecting two
    communities is typically high, as many of the shortest paths
    between nodes in separate communities go through them. So we
    gradually remove the edge with the highest betweenness and
    recalculate the betweennesses after every removal. This way sooner
    or later the network falls of to separate components. The result
    of the clustering will be represented by a dendrogram.


    q: int, optional
        Number of communities to detect.  This, in practice, specifies
        the level at which to cut the dendogram.  If this is not
        given, return maximal modularity cut.  FIXME: make this work


    References:

    Girvan, Michelle, and Mark EJ Newman. 'Community structure in
    social and biological networks.' Proceedings of the National
    Academy of Sciences 99.12 (2002): 7821-7826.
    http://arxiv.org/abs/cond-mat/0112110v1
    """
    # dendogram
    _method = 'community_edge_betweenness'
    _kwargs = dict(directed=False)
    _argmap = dict(q='clusters')
    _weightable = True
    q = None
class IgraphSpinglass(_IGraphAlgorithm):
    """Potts model (spin glass) community detection

    Finds the community structure of the graph according to the
    spinglass community detection method of Reichardt & Bornholdt.


    q: int, default 25
        Number of communities to detect.  This is an upper limit, not
        a fixed number.

    synchronous: bool
        Syncronous updates?  If True, update all spins at the same
        time, otherwise one-by-one.

    start_temp: float
        Start temperature of simulated annealing

    stop_temp: float
        Stop temperature of simulated annealing.

    cool_fact: float
        Simulated annealing temperature step multiplier.

    update_rule: str
        Null model for graph.  'config' for configuration model,
        'simple' for random graph.

    gamma: float
        Relative trade-off between internal and external interactions.
        High gamma produces smaller communities.  See references for
        details.

    implementation: str
        Implementation to use.  'orig' is faster but can not handle
        negative weights, or 'neg' can handle negative weights.

    lambda_: float
        Like gamma for the 'neg' implementation.


    References:

    J. Reichardt and S. Bornholdt: Statistical mechanics of community
    detection. Phys. Rev. E. 74:016110 (2006).
    http://arxiv.org/abs/cond-mat/0603718

    V.A. Traag and J. Bruggeman: Community detection in networks with
    positive and negative links. Phys Rev E 80:036115 (2009).
    http://arxiv.org/abs/0811.2329
    """
    _weightable = True
    _method = 'community_spinglass'
    _argmap = dict(q='spins', synchronous='parupdate',
                   start_temp='start_temp', stop_temp='stop_temp', cool_fact='cool_fact',
                   update_rule='update_rule', gamma='gamma', implementation='implementation',
                   lambda_='lambda_',
                   )
    q = 25   # upstream default
    synchronous = False
    start_temp = 1
    stop_temp = 0.1
    cool_fact = .99
    update_rule = 'config' # 'config' (configuration model) or 'simple' (random graph)
    gamma = 1
    implementation = 'orig' # or 'neg' for negative weights
    lambda_ = 1
class IgraphWalktrap(_IGraphAlgorithmDendogram):
    """Walktrap

    Community detection algorithm of Latapy & Pons, based on random
    walks.

    The basic idea of the algorithm is that short random walks tend to
    stay in the same community. The result of the clustering will be
    represented as a dendrogram.


    steps: int, default 4
        Number of steps in the random walk to take.


    References:

    Pascal Pons and Matthieu Latapy: Computing communities in large
    networks using random walks.  Computer and Information
    Sciences-ISCIS 2005. Springer Berlin Heidelberg, 2005. 284-293.
    http://arxiv.org/abs/physics/0512106
    """
    # dendogram
    _method = 'community_walktrap'
    steps = 4
    _weightable = True
    _argmap = dict(steps='steps')





def get(name, search=()):
    """Get an algorithm, optionally searching different namespaces.

    This is a convenience function to get an algorithm of a certain
    name.  It will first search all namespaces given in the `search`
    argument, and then this module itself.  The purpose of this is to
    provide a convenience method of getting algorithms when you have
    made some custom subclasses in your own module.

    name: string
        Algorithm name to get
    search: dict or iterable:
        Dict of other namespace to search.  Can also be an iterable of
        dicts.
    """
    cda = None
    if isinstance(search, dict):
        search = [search]
    for ns in search:
        if name in ns:
            cda = ns[name]
    if cda is None:
        cda = globals()[name]
    return cda

if __name__ == "__main__":
    if len(sys.argv) < 4 or '-h' in sys.argv or '--help' in sys.argv:
        print "%s MethodName infile outfile [options-dict]"%os.path.basename(sys.argv[0])
        print
        print "available methods:", ', '.join(sorted(
            name for name, cda in globals().iteritems()
            if isinstance(cda, type)
            and issubclass(cda, CDMethod)
            and name[0]!='_'))
        if len(sys.argv) > 1 and sys.argv[1] not in ('-h', '--help'):
            method = sys.argv[1]
            method = locals()[method]
            print
            print "Available options for %s:"%method.name()
            print "\n".join(("%s=%s"%(x,getattr(method,x))
                                      if pcd.util._isTrivialType(getattr(method,x))
                                      else x)
                             for x in sorted(dir(method))
                             if (not x.startswith('_')
                                 and not hasattr(getattr(method, x), '__call__')))
        else:
            print "Method help: %s MethodName -h"%os.path.basename(sys.argv[0])
        sys.exit()
    #import optparse
    method = sys.argv[1]
    input = sys.argv[2]
    output = sys.argv[3]
    if len(sys.argv) > 4:
        options = sys.argv[4]
        from pcd.util import leval
        options = leval(options)
    else:
        options = {}

    import pcd.ioutil
    g = pcd.ioutil.read_any(input)
    #print networkx.get_edge_attributes(g, 'weight')
    #print g.adj['51']['68']
    print len(g), g.number_of_edges()
    print [len(_) for _ in networkx.connected_components(g)]
    #from fitz import interactnow

    options['basename'] = os.path.basename(input)

    method = locals()[method]
    r = method(g, **options)

    for i, cmtys in enumerate(r.results):
        cmtys.write_clusters(output+'.%d.txt'%i)


    #initial = r.cmtys
    #r = method(g, initial=initial, **options)
    #import pcd.ioutil
    #pcd.ioutil.write_pajek(input+'.cmtys.net', g, r.cmtys)

    #g2 = g.copy()
    #r.cmtys.load_networkx_custom(g2, attrname='cmty')
    #networkx.write_gml(g2, options['basename']+'.gml')

    #from fitz import interactnow


