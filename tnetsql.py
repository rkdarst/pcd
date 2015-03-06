import re
import sqlite3


def _t_where(tmin=None, tmax=None):
    if not tmin and not tmax:
        return '1'
    s = [ ]
    if tmin:
        s.append('%s<=t'%tmin)
    if tmax:
        if s: s.append('and')
        s.append('t<=%s'%tmax)
    return " ".join(s)

class TNet(object):
    def __init__(self, fname=':memory:', make_indices=True):
        self.conn = sqlite3.connect(fname)
        if make_indices:
            self._make_indices
        c = self.conn.cursor()
        c.execute('''CREATE TABLE if not exists edge
                     (a int not null, b int not null, t int)''')
        c.close()
    def _make_indices(self):
        c = self.conn.cursor()
        #c.execute('''CREATE INDEX if not exists index_edge_a ON edge (a)''')
        #c.execute('''CREATE INDEX if not exists index_edge_b ON edge (b)''')
        c.execute('''CREATE INDEX if not exists index_edge_t ON edge (t)''')
        c.execute('''CREATE INDEX if not exists index_edge_t_a ON edge (t, a)''')
        c.execute('''CREATE INDEX if not exists index_edge_t_b ON edge (t, b)''')
        #c.execute('''CREATE INDEX if not exists index_edge_a_t ON edge (a, t)''')
        #c.execute('''CREATE INDEX if not exists index_edge_b_t ON edge (b, t)''')
        c.close()
    def add_edge(self, a, b, t):
        c = self.conn.cursor()
        c.execute("INSERT INTO edge VALUES (?, ?, ?)", (a, b, t))
        self.conn.commit()
    def add_edges(self, it):
        c = self.conn.cursor()
        c.executemany("INSERT INTO edge VALUES (?, ?, ?)", it)
        self.conn.commit()
    def dump(self):
        c = self.conn.cursor()
        c.execute('''SELECT * FROM edge''')
        for row in c:
            print row
    def _execute(self, *args):
        """Make new cursor, execute SQL, return cursor"""
        c = self.conn.cursor()
        c.execute(*args)
        return c

    def __len__(self):
        c = self.conn.cursor()
        c.execute("SELECT count(*) from edge")
        return c.fetchone()[0]

    def __getitem__(self, name):
        c = self.conn.cursor()
        c.execute('''select b, t from edge where a=?''', (name, ))
        for b, t in c:
            yield b, t
    def t_min(self):
        c = self.conn.cursor()
        c.execute("SELECT min(t) from edge")
        return c.fetchone()[0]
    def t_max(self):
        c = self.conn.cursor()
        c.execute("SELECT max(t) from edge")
        return c.fetchone()[0]


    #
    # Input and output methods
    #
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

    @property
    def t(self):
        tnet = self
        class TimeSliceProxy(object):
            def __getitem__(self, slc):
                c = tnet.conn.cursor()
                if not isinstance(slc, slice):
                    c.execute('''select a, b, t from edge where t == ?''', (slc, ))
                    return c
                start = slc.start
                stop = slc.stop

                if stop is None:  # lower bound only
                    c.execute('''select a, b, t from edge where ? <= t''', (start, ))
                elif start is None:  # upper bound only
                    c.execute('''select a, b, t from edge where t < ?''', (stop, ))
                else:
                    c.execute('''select a, b, t from edge where ? <= t AND t < ?''', (start, stop, ))
                return c

        return TimeSliceProxy()


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

        Nodes are not repeated after multiple edges.

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
