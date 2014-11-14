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

    # Output
    def write_tcommlist(self, fname):
        f = open(fname, 'w')
        for t, cmtys in self.iteritems():
            for c, nodes in cmtys.iteritems():
                for n in nodes:
                    print >> f, t, n, c

