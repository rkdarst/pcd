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
    def __init__(self, fname=':memory:'):
        self.conn = sqlite3.connect(fname)
        c = self.conn.cursor()

        c.execute('''CREATE TABLE if not exists edge
                     (a int not null, b int not null, t int, w int)''')
        c.execute('''CREATE INDEX if not exists index_edge_a ON edge (a)''')
        c.execute('''CREATE INDEX if not exists index_edge_b ON edge (b)''')
        c.execute('''CREATE INDEX if not exists index_edge_t ON edge (t)''')
        c.close()
    def add_edge(self, a, b, t, w=1.0):
        c = self.conn.cursor()
        c.execute("INSERT INTO edge VALUES (?, ?, ?, ?)", (a, b, t, w))
        self.conn.commit()
    def add_edges(self, it):
        c = self.conn.cursor()
        c.executemany("INSERT INTO edge VALUES (?, ?, ?, ?)", it)
        self.conn.commit()
    def dump(self):
        c = self.conn.cursor()
        c.execute('''SELECT * FROM edge''')
        for row in c:
            print row

    def __len__(self):
        c = self.conn.cursor()
        c.execute("SELECT count(*) from edge")
        return c.fetchone()[0]

    def __getitem__(self, name):
        c = self.conn.cursor()
        c.execute('''select b, t, w from edge where a=?''', (name, ))
        for b, t, w in c:
            yield b, t, w

    @property
    def t(self):
        tnet = self
        class TimeSliceProxy(object):
            def __getitem__(self, slc):
                c = tnet.conn.cursor()
                if not isinstance(slc, slice):
                    c.execute('''select a, b, t, w from edge where t == ?''', (slc, ))
                    return c
                start = slc.start
                stop = slc.stop

                if stop is None:  # lower bound only
                    c.execute('''select a, b, t, w from edge where ? <= t''', (start, ))
                elif start is None:  # upper bound only
                    c.execute('''select a, b, t, w from edge where t < ?''', (stop, ))
                else:
                    c.execute('''select a, b, t, w from edge where ? <= t AND t < ?''', (start, stop, ))
                return c

        return TimeSliceProxy()


    def events_after(self, t):
        """Iterate over all events after a certain time.

        Returns: iterator over (a,b,t,w)
        """
        c = self.conn.cursor()
        c.execute('''SELECT a, b, t, w FROM edge WHERE t>? ORDER BY t ASC''',
                  (t, ))
        return c

    def iter_out_events_after(self, a, t):
        """Iterate out events after time

        Returns: iterator over (b, t, w)
        """
        c = self.conn.cursor()
        c.execute('''SELECT b, t, w FROM edge WHERE a=? and ?<t ORDER BY t ASC
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
