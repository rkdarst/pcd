"""
Utilities for interacting with sqlite3
"""


import cPickle as pickle
class SQLiteDict(object):
    def __init__(self, fname):
        import sqlite3
        self.conn = sqlite3.connect(fname)
        c = self.conn.cursor()
        c.execute('create table if not exists dict (hash int unique, key blob, value blob)')
        c.execute('create index if not exists dict_index on dict (hash)')
        self.conn.commit()
    def __setitem__(self, name, value):
        h = hash(name)
        c = self.conn.cursor()
        #c.execute('select name from dict where hash==? ', h)
        #if len(c.fetchhall()) > 0:
        c.execute('insert or replace into dict values (?,?,?)',
                  (h, buffer(pickle.dumps(name, -1)),
                   buffer(pickle.dumps(value, -1))))
        self.conn.commit()
    def __getitem__(self, name):
        h = hash(name)
        c = self.conn.cursor()
        c.execute('select hash,key,value from dict where hash==?', (h,))
        row = c.fetchone()
        if row is None:
            raise KeyError("%s"%name)
        assert name == pickle.loads(str(row[1])), "Hash-object mismatch"
        return pickle.loads(str(row[2]))
    def __contains__(self, name):
        h = hash(name)
        c = self.conn.cursor()
        c.execute('select hash,key,value from dict where hash==?', (h,))
        row = c.fetchone()
        if row is None:
            return False
        assert name == pickle.loads(str(row[1])), "Hash-object mismatch"
        return True
