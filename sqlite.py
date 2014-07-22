"""
Utilities for interacting with sqlite3.
"""

import cPickle as pickle

class SQLiteDict(object):
    """On-disk key-value store with SQLite.

    This class operates similar to the 'shelve' module, but uses
    sqlite3 as a backend for better locking ability and possible
    performance.  Arbitrary python objects supported.  Keys must be
    hashable and both keys and values must be pickleable.

    Usage:

    cache = SQLiteDict('filename.sqlite')
    cache[key] = value

    if key in cache:
        ....

    Most dictionary operations are supported and the rest should be
    added.

    """
    def __init__(self, fname):
        import sqlite3
        self.conn = sqlite3.connect(fname)
        c = self.conn.cursor()
        c.execute('create table if not exists dict (hash int unique, key blob, value blob)')
        c.execute('create index if not exists dict_index on dict (hash)')
        self.conn.commit()
    def __len__(self):
        c = self.conn.cursor()
        c.execute('select count(*) from dict')
        row = c.fetchone()
        c.close()
        return row[0]
    def __setitem__(self, name, value):
        h = hash(name)
        # Make sure that if the hash is already in the database, it
        # refers to the same.
        key_p = buffer(pickle.dumps(name,  -1))
        val_p = buffer(pickle.dumps(value, -1))
        with self.conn:
            c = self.conn.cursor()
            c.execute('select hash,key,value from dict where hash==?', (h,))
            row = c.fetchone()
            if row is not None:
                assert name == pickle.loads(str(row[1])), "Hash-object mismatch: %r, %r"%(name, pickle.loads(str(row[1])))
            # Set value in database.
            c.execute('insert or replace into dict values (?,?,?)',
                      (h, key_p, val_p))
            c.close()
    def __getitem__(self, name):
        """Get item, or KeyError if not present"""
        h = hash(name)
        c = self.conn.cursor()
        c.execute('select hash,key,value from dict where hash==?', (h,))
        row = c.fetchone()
        if row is None:
            raise KeyError("%s"%name)
        c.close()
        assert name == pickle.loads(str(row[1])), "Hash-object mismatch: %r, %r"%(name, pickle.loads(str(row[1])))
        return pickle.loads(str(row[2]))
    def get(self, name, default=None):
        """Get item, or return default if missing."""
        try:
            return self[name]
        except KeyError:
            return default
    def __contains__(self, name):
        """True if self cotains key of 'name'"""
        h = hash(name)
        c = self.conn.cursor()
        c.execute('select hash,key,value from dict where hash==?', (h,))
        row = c.fetchone()
        if row is None:
            return False
        c.close()
        if row[0] != hash(pickle.loads(str(row[1]))):
            print 'h'%40, 'corruption, hash mismatch'
        assert name == pickle.loads(str(row[1])), "Hash-object mismatch: %d %r\n%d %r"%(h, name, hash(pickle.loads(str(row[1]))), pickle.loads(str(row[1])))
        return True
    has_key = __contains__
    def iterkeys(self):
        """Iterate over keys"""
        c = self.conn.cursor()
        c.execute('select key from dict')
        for row in c:
            yield pickle.loads(str(row[0]))
        c.close()
    __iter__ = iterkeys
    def itervalues(self):
        """Iterate over values"""
        c = self.conn.cursor()
        c.execute('select value from dict')
        for row in c:
            yield pickle.loads(str(row[0]))
        c.close()
    def iteritems(self):
        """Iterate over (key,value) tuples"""
        c = self.conn.cursor()
        c.execute('select key,value from dict')
        for row in c:
            yield pickle.loads(str(row[0])), pickle.loads(str(row[1]))
        c.close()
    def keys(self):
        """List of keys."""
        return list(self.iterkeys())
    def values(self):
        """List of values."""
        return list(self.itervalues())
    def items(self):
        """List of (key,value) pairs."""
        return list(self.iteritems())
    def __delitem__(self, key):
        h = hash(name)
        with self.conn:
            c = self.conn.cursor()
            # Verify that the object in the DB matches
            c.execute('select hash,key,value from dict where hash==?', (h,))
            row = c.fetchone()
            if row is None:
                raise KeyError("%s"%name)
            assert name == pickle.loads(str(row[1])), "Hash-object mismatch"
            # Do the actual deletion
            c.execute('delete from dict where hash==?', (h,))
            c.close()
        return pickle.loads(str(row[2]))
    def update(self, E):
        """Same as dict.update"""
        if hasattr(E, 'keys') or hasattr(E, 'iterkeys'):
            for k, v in E.iteritems():
                self[k] = v
        else:
            for (k,v) in E:
                self[k] = v

    def _find_corruption(self, remove=False):
        """Find keys with hash-object mismatches.

        All keys are stored by their hash values.  If the database is
        corrupted, then some of the hash values may be wrong.

        remove: bool, default False
            If true, remove keys with hash """
        is_corrupted = False
        corrupted_keys = [ ]

        c = self.conn.cursor()
        c.execute('select hash,key,value from dict')
        # Find corrupted keys
        for row in c:
            h = row[0]
            k = pickle.loads(str(row[1]))
            if h != hash(k):
                is_corrupted = True
                print h, hash(k), k
                if remove:
                    corrupted_keys.append((h, row[1]))
        # Remove the currupted values:
        if remove:
            with self.conn:
                c.executemany('delete from dict where hash==? and key==?',
                              corrupted_keys)
        c.close()
        return is_corrupted
    def _remove_nonunpickleable(self):
        c = self.conn.cursor()
        c.execute('select hash,key,value from dict')
        for row in c:
            k = None
            try:
                h = row[0]
                k = pickle.loads(str(row[1]))
                v = pickle.loads(str(row[2]))
            except pickle.UnpicklingError:
                print h, k if k is not None else row[1]
                corrupted_keys.append((h, row[1]))
        with self.conn:
            c.executemany('delete from dict where hash==? and key==?',
                          corrupted_keys)
        c.close()


if __name__ ==  '__main__':
    import sys

    d = SQLiteDict(sys.argv[1])
    cmd = sys.argv[2]
    arg = None
    if len(sys.argv) > 3:
        arg = sys.argv[3]
    if cmd == 'list':
        for key in d.iterkeys():
            if arg and not eval(arg):
                    continue
            print repr(key)

    if cmd == 'del':
        for key in d.iterkeys():
            if not eval(arg):
                continue
            print repr(key)
            #del d[key]
