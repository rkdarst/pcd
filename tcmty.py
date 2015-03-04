"""Temporal communities
"""

from . import cmty

from pcd.util import leval

class TemporalCommunities(object):

    def __init__(self, tcmtys):
        self._tcmtys = tcmtys

    # Basic operatoins

    def __iter__(self):
        for t in sorted(self._tcmtys):
            yield t
    def __getitem__(self, t):
        return self._tcmtys[t]
    def __len__(self):
        return self._tcmtys.__len__()
    def iteritems(self):
        for t in self:
            yield t, self[t]
    def __str__(self):
        return '<%s object with Nt=%s at 0x%x>'%(
            self.__class__.__name__,
            len(self),
            id(self))

    # Modification
    def add_time(self, t, cmtys):
        assert t not in self._tcmtys
        self._tcmtys[t] = cmtys

    # Input
    @classmethod
    def open_any(cls, fname, time_type=leval, node_type=leval,
                 cmty_type=leval):
        """Automatically open file, detecting type via contents.
        """
        lines = open(fname).read(1024).split()
        lines = [ x.strip() for x in lines ]
        lines = [ x for x in lines if len(x)>0 and not x.startswith('#') ]
        if len(lines[0].split()) >= 5:
            return cls.from_tmatrix(fname)
        else:
            return cls.from_tcommlist(fname,
                                      time_type=time_type,
                                      node_type=node_type,
                                      cmty_type=cmty_type)
    from_auto = open_any
    @classmethod
    def from_tcommlist(cls, fname, time_type=leval, node_type=leval,
                       cmty_type=leval):
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
            n = node_type(n)
            c = cmty_type(c)

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
                     time_list=None, cmty_type=leval):
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

            t = time_counter
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
    def write_tmatrix(self, fname):
        f = open(fname, 'w')
        for t, cmtys in self.iteritems():
            nodecmtys = cmtys.nodecmtys_onetoone()
            nodes = sorted(nodecmtys)
            print >> f, ' '.join([str(nodecmtys[n]) for n in nodes])



if __name__ == '__main__':
    import sys
    tcmtys = TemporalCommunities.open_any(sys.argv[1])
    print tcmtys
    from fitz import interactnow
