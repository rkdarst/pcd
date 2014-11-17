"""Temporal communities
"""

from . import cmty


class TemporalCommunities(object):

    def __init__(self, tcmtys):
        self._tcmtys = tcmtys

    # Basic operatoins

    def __iter__(self):
        for t in sorted(self._tcmtys):
            yield t
    def __getitem__(self, t):
        return self._tcmtys[t]
    def iteritems(self):
        for t in self:
            yield t, self[t]

    # Modification
    def add_time(self, t, cmtys):
        assert t not in self._tcmtys
        self._tcmtys[t] = cmtys

    # Input
    @classmethod
    def from_tcommlist(cls, fname, time_type=float, node_type=str,
                       cmty_type=str):
        """Load temporal communities from a tcommlist.

        Tcommlist consists of successive lines of this format:
          time-value node-id cmty-id.
        """
        tcmtys = { }
        f = open(fname)
        for line in f:
            line = line.strip()
            if line.startswith('#'): continue
            t, n, c = line.split()
            t = time_type(t)
            n = time_type(n)
            c = time_type(c)

            if t not in tcmtys:
                cmtys = tcmtys[t] = { }
            else:
                cmtys = tcmtys[t]
            if c not in cmtys:
                nodeset = cmtys[c] = set()
            else:
                nodeset = cmtys[c]
            nodeset.add(n)

        self = cls(dict((t, cmty.Communities(cmtynodes))
                        for t, cmtynodes in tcmtys.iteritems()))
        return self
    @classmethod
    def from_tmatrix(cls, fname, node_list=None,
                     time_list=None):
        """Load temporal communities from a temporal matrix.

        Node IDs and times are not included.  Communitties must be
        covering and non-overlapping.  Each line consistes of the
        communities, of nodes range(0, N).  Community IDs start at 0
        and advance incrementally.
        """
        tcmtys = { }
        f = open(fname)
        time_counter = 0
        for line in f:
            line = line.strip()
            if line.startswith('#'): continue

            t = time
            time_counter += 1
            if t not in tcmtys:
                cmtys = tcmtys[t] = { }
            else:
                cmtys = tcmtys[t]

            node_cmtys = (cmty_type(x) for x in line.split())
            for n, c in enumerate(node_cmtys):
                if c not in cmtys:
                    nodeset = cmtys[c] = set()
                else:
                    nodeset = cmtys[c]
                nodeset.add(n)

        self = cls(dict((t, cmty.Communities(cmtynodes))
                        for t, cmtynodes in tcmtys.iteritems()))
        return self

    # Output
    def write_tcommlist(self, fname):
        f = open(fname, 'w')
        for t, cmtys in self.iteritems():
            for c, nodes in cmtys.iteritems():
                for n in nodes:
                    print >> f, t, n, c

