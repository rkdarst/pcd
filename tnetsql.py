# Richard Darst, April 2015
"""Temporal network structures in sqlite3.

This module contains a structure suitable for holding temporal
networks.  All data is loaded into an sqlite3 database, which can be
either in-memory or on-disk.  Then, events can be efficiently queried
in a variety of ways.

The basic structure contains four fields: a, b, t, dt.  These
correspond to first event, second event, time, and time duration.  The
'dt' field is not implemented yet.  Generally, `a` and `b` are assumed
to be ints, and `t` is assumed to be a float.  However, this is not a
restriction of sqlite, and most code here is general enough to handle
any input formats, but this is not tested or an explicit goal (yet).

The network can also be queried from the command line using the
`sqlite3` command line tool, if saved to a file.  This allows useful
investigation and analytics without having to use Python.

For example usage, see test_tnetsql.py
"""

import re
import sqlite3

# Modified version of ast.literal_eval.  This turns strings into ints,
# floats, or strings depending on what data is actually contained in
# there.  Used to avoid restricting ourselves to ints or floats in
# certain cases.
from pcd.util import leval

def _t_where(tmin=None, tmax=None):
    """Return a WHERE clause filtering on time:

    Returns '1' if there are no parameters given, or 'tmin<=t' and/or
    't<tmax' properly joined by AND."""
    if not tmin and not tmax:
        return '1'
    s = [ ]
    if tmin:
        s.append('%s<=t'%tmin)
    if tmax:
        if s: s.append('and')
        s.append('t<%s'%tmax)
    return " ".join(s)


class TNet(object):
    """Temporal network object.

    This object is initiated with a filename of ':memory:' (default)
    or some file on disk.  If a disk filename is given, all operations
    will be synced to disk, and if an object is re-created with the
    same filename, all events will automatically be available with
    zero start-up cost.
    """
    _indices = [('t', ),
                ('t', 'a',),
                ('t', 'b',),
                ('a', 't',),
                ('b', 't',),
                ('a', 'b',),
                ]
    def __init__(self, fname=':memory:', make_indices=True):
        """Init.

        fname: see clas docstring

        make_indices: if true, create indices.  Could be set to false
        for delayed index generation (call _make_indices() yourself
        later).
        """
        self.conn = sqlite3.connect(fname)
        if make_indices:
            self._make_indices
        c = self.conn.cursor()
        c.execute('''CREATE TABLE if not exists edge
                     (a int not null, b int not null, t int)''')
        c.close()
    def _make_indices(self):
        """Create database indices, if they don't exist already.


        Each index is efficient for one form of lookup, however, each
        index slows down graph creation and takes more memory.
        Ideally, only the indices used will be created.  We default to
        creating all possible indices.
        """
        c = self.conn.cursor()
        # These statements left as an example of index creation.
        #c.execute('''CREATE INDEX if not exists index_edge_t ON edge (t)''')
        #c.execute('''CREATE INDEX if not exists index_edge_t_a ON edge (t, a)''')
        #c.execute('''CREATE INDEX if not exists index_edge_t_b ON edge (t, b)''')

        # This allows index creation to be parameterized by the
        # _indices attribute.  This could be changed in subclasses to
        # suit a particular problem, for example.
        for fields in self._indices:
            c.execute('''CREATE INDEX if not exists '''
                      '''index_edge_%s ON edge (%s)'''%(
                          ('_'.join(fields), ', '.join(fields)),
                          ))

        c.close()
    def dump(self):
        """Print entire network to stdout.

        This is equivalaent to dubping the 'edge' table."""
        c = self.conn.cursor()
        c.execute('''SELECT * FROM edge''')
        for row in c:
            print row
    def _execute(self, *args):
        """Make new cursor, execute SQL, return cursor"""
        c = self.conn.cursor()
        c.execute(*args)
        return c


    # Modifying network structure
    def add_edge(self, a, b, t):
        """Add one edge to graph, a->b at time t"""
        c = self.conn.cursor()
        c.execute("INSERT INTO edge VALUES (?, ?, ?)", (a, b, t))
        self.conn.commit()
    def add_edges(self, it):
        """Add edges to graph from iterator.

        Each item in the iterator should return a (a,b,t) triple which
        is added as one edge to graph, a->b at time t."""
        c = self.conn.cursor()
        c.executemany("INSERT INTO edge VALUES (?, ?, ?)", it)
        self.conn.commit()

    # basic information about network.
    def __len__(self):
        """Number of events in network"""
        c = self.conn.cursor()
        c.execute("SELECT count(*) from edge")
        return c.fetchone()[0]

    def __getitem__(self, name):
        """get all neighbors of node 'a'.

        Returns an iterator over (neighbor, time) pairs."""
        c = self.conn.cursor()
        c.execute('''select b, t from edge where a=?''', (name, ))
        for b, t in c:
            yield b, t
    def t_min(self):
        """Time of first event in network"""
        c = self.conn.cursor()
        c.execute("SELECT min(t) from edge")
        return c.fetchone()[0]
    def t_max(self):
        """Time of last event in network"""
        c = self.conn.cursor()
        c.execute("SELECT max(t) from edge")
        return c.fetchone()[0]


    #
    # Input and output methods
    #
    def from_eventlist(self, f):
        """Read networks from list of (a, b, t) in a file.

        Unknown what the source of this format is.
        """
        re_header = re.compile(r'^time\s+([0-9.-]+),\s+([0-9]+) contact')
        header = next(f)
        def _iter_events():
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue
                line = line.split()
                yield leval(line[0]), leval(line[1]), leval(line[2])
        self.add_edges(_iter_events())
    def from_ts(self, f):
        """Read networks from .ts format

        Unknown what the source of this format is.
        """
        re_header = re.compile(r'^time\s+([0-9.-]+),\s+([0-9]+) contact')
        header = next(f)
        def _iter_events():
            i = 0
            for line in f:
                #print line.strip()
                m = re_header.match(line)
                t = int(m.group(1))
                N = int(m.group(2))
                for _ in range(N):
                    nextline = next(f)
                    #print '*', nextline
                    a, b = nextline.strip().split()
                    yield a, b, t
                #    i += 1
                #if i > 500000: break
        self.add_edges(_iter_events())

    # Different iterators of network structure.
    def static(self):
        """Return entire network collapased to a static network.

        Returns: iterator over (a, b) edges, where edges are directed
        (a and b both exist) and there may be self-loops"""
        c = self.conn.cursor()
        c.execute('''SELECT DISTINCT a, b FROM edge''')
        return c
    def events_after(self, t):
        """Iterate over all events after a certain time.

        Returns: iterator over (a,b,t,w)
        """
        c = self.conn.cursor()
        c.execute('''SELECT a, b, t FROM edge WHERE t>? ORDER BY t ASC''',
                  (t, ))
        return c

    def iter_out_events_after(self, a, t):
        """Iterate out events after time

        Returns: iterator over (b, t)
        """
        c = self.conn.cursor()
        c.execute('''SELECT b, t FROM edge WHERE a=? and ?<t ORDER BY t ASC
                  ''', (a, t))
        return c
    def successors_after(self, a, t, tmax):
        """Iterate over (out) connected nodes after time.

        Nodes are not repeated after multiple edges, unlike the
        `events*` methods.

        Returns: iterator over (b, )
        """
        c = self.conn.cursor()
        c.execute('''SELECT DISTINCT b FROM edge WHERE a=? and ?<t and t<=?
                  ORDER BY t ASC
                  ''', (a, t, tmax))
        return c
    def successors_after_t(self, a, t, tmax):
        """Iterate over (out) connected nodes after time.

        Nodes are not repeated after multiple edges.

        Returns: iterator over (b, t)
        """
        c = self.conn.cursor()
        c.execute('''SELECT DISTINCT b FROM edge WHERE a=? and ?<t and t<=?
                  ORDER BY t ASC
                  ''', (n, t, tmax))
        return c
    # does an edge exist?
    def has_edge(self, a, b):
        """Is there an event from a to b at any time?

        This is treated as a directed test.
        """
        c = self.conn.cursor()
        c.execute('''SELECT rowid from edge WHERE a=? and b=? LIMIT 1
                  ''', (a, b))
        return c.fetchone() is not None
    # does an edge exist in the next dt time?
    def has_edge_interval(self, a, b, tmin, tmax):
        """Is there an event from a to b at tmin<=t<tmax"""
        c = self.conn.cursor()
        c.execute('''SELECT rowid FROM edge WHERE a=? and b=? and ?<t and t<=?
                     LIMIT 1
                  ''', (a, b, tmin, tmax))
        return c.fetchone() is not None

    def iter_neighs_after(self, t, n, reverse=True):
        raise NotImplementedError

    def iter_out_after(self, t, n, reverse=True):
        raise NotImplementedError


    # all neighbors of node


