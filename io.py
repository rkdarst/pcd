
"""Functions for dealing with reading/writing graphs.

This function is designed to compliment, not replace, the networkx
reading/writing functions.

TODO: handle compressed files."""

import os
import re

import networkx

class UnknownTypeError(Exception):
    pass

def read_any(fname):
    """Attempt to read file in any format.

    - First, file name can be specified as a schema prefix, as in
      'gml:FILE_NAME'.

    - Then, it can be specified via filename extension.

    - Then, it can be specified via filename extension after .bz2 or
      .gz.  Note: We do not decompress it, but rely on networkx to do
      decompression.

    - Then, look at the first 512 bytes to see if the file is either
      in pajek or gml format.

    - If read is unsuccessful, raise UnknownTypeError.

    If read is successful, return a networkx graph object.

    A list of all readers in networkx:
    import networkx
    [ x for x in dir(networkx) if x.startswith('read_') ]
    """#%'    \n'.join(x for x in dir(networkx) if x.startswith('read_'))
    # Try a schema prefix based match:
    m =  re.match(r'^(\w+?):(.*)$', fname)
    if m:
        schema, fname = m.group(1), m.group(2)
        # Look for networkx.read_schema.  If it exists, use that to
        # read and return.
        reader = _get_reader(schema)
        if reader is not None:
            return reader(fname)
    # Try known file suffix based reading:
    base, ext = os.path.splitext(fname)
    #if ext == '.gml':
    #    return networkx.read_gml(fname)
    # look for any reader in networkx with this extension.
    reader = _get_reader(ext[1:])
    if reader is not None:
        return reader(fname)
    # If file is compressed, look at true extension, see if it has a
    # reader.
    if ext in ('.gz', '.bz2'):
        newext = os.path.splitext(base)
        reader = _get_reader(newext[1:])
        if reader is not None:
            return reader(fname)
    # Look inside for the schema:
    if _test_pajek(fname):    return networkx.read_pajek(fname)
    if _test_gml(fname):      return networkx.read_gml(fname)
    if _test_edgelist(fname): return networkx.read_edgelist(fname)

    # So it is probably an edgelist.

    # Raise exception
    raise UnknownTypeError("Can't open %s"%fname)


def _get_reader(name):
    """Return a reader for the format.

    Current search path:
    - networkx.read_NAME
    """
    return getattr(networkx, 'read_'+name, None)

def _test_pajek(fname):
    """Does this look like a pajek file?

    Pajek files have a '*Vertices' string (case insensitive
    search) in the first 512 bytes.  It should match this regex:

    ^\*Vertices[ \t]+[0-9]+
    """
    data = open(fname).read(512)
    if re.search(r'^\*Vertices[ \t]+[0-9]+', data, re.I|re.M):
        return True
def _test_gml(fname):
    """Does this look like a gml file?

    gml files have this regex in the first 512 bytes:

    '^[ \t]*graph\s[\s*$'
    """
    data = open(fname).read(512)
    if re.search(r'^[ \t]*graph\s\[[ \t]*$', data, re.I|re.M):
        return True
def _test_edgelist(fname):
    """Does this look like an edgelist file?

    edgelist files match this regex in the first 512 bytes:

    ^(?!#)\S+[ \t]+\S+([ \t]+[-0-9.e]+([ \t]+(.*?))?)?[ \t]*$

    translation:
    - beginning of line that doesn't start with '#':  ^(?!#)
    - two non-whitespace groups(no space in front):  \S+\s+\S+
    - optionally:
      - space, then a weight:   (\s+[-0-9.e]+ ....  )?
      - optionally, if weight exists:
        - space, any other text trailing:   (\s+(.*?))?
    - any trailing whitespace, end of line:  \s*$
    """
    data = open(fname, 'U').read(512)
    #lines = data.split('\n')
    #lines = [ l.strip() for l in lines if l.strip() and l[0]!='#']
    if re.search(r'^(?!#)\S+[ \t]+\S+([ \t]+[-0-9.e]+([ \t]+(.*?))?)?[ \t]*$', data, re.I|re.M):
        return True


import gzip
import bz2
def zopen(fname, mode='r', compress=None):
    """Automatic compression and decompression.

    If compress='gz' or 'bz2', write file with that compression.  If
    the file does not have the proper extension, add it to the
    filename.

    If filename ends in .gz or .bz2, automatically compress/decompress
    on opening.

    If filename doesn't exist, look for the same name with .gz or .bz2
    and open that if it exists (but only if exactly one of them
    exists).
    """
    # Find the true filename, if we are reading.  You can pass the
    # name 'x.txt' and it will find 'x.txt.gz' or 'x.txt.bz2'
    # automatically.
    if 'r' in mode and not os.path.exists(fname):
        gz_exists  = os.path.exists(fname+'.gz')
        bz2_exists = os.path.exists(fname+'.bz2')
        assert not (gz_exists and bz2_exists)
        if gz_exists:   fname = fname+'.gz'
        if bz2_exists:  fname = fname+'.bz2'
    # if compression specified, make sure filename ends in that.
    if compress:
        if compress == 'gz':
            if not fname.endswith('.gz'):
                fname = fname + '.gz'
        if compress == 'bz2':
            if not fname.endswith('.bz2'):
                fname = fname + '.bz2'
    # if compression not speecified, see if we can infer one from the
    # extension.
    else:
        if fname.endswith('.gz'):
            compress = 'gz'
        elif fname.endswith('.bz2'):
            compress = 'bz2'
    # Decide on the opening interface
    if compress == 'gz':
        compressor = gzip.GzipFile
    elif compress == 'bz2':
        compressor = bz2.BZ2File
    else:
        compressor = open

    if 'r' in mode:
        f = compressor(fname, mode=mode)
    elif 'w' in mode:
        f = compressor(fname, mode=mode)
    else:
        raise ValueError("We can't handle mode %s"%mode)
    return f
def _test_zopen(tmp='/tmp'):
    from os.path import join
    import random

    # Test writing compression.
    tmpname = random.uniform(1, 2**32-1)
    name1 = join(tmp, 'pcdiotestzopen-A-%d'%tmpname)
    zopen(name1+'.bz2', 'w').write('test data')
    assert zopen(name1, compress='bz2').read() == 'test data'
    os.unlink(name1+'.bz2')

    # Test automatically adding the extension
    tmpname = random.uniform(1, 2**32-1)
    name2 = join(tmp, 'pcdiotestzopen-B-%d'%tmpname)
    zopen(name2, 'w', compress='gz').write('test data2')
    assert zopen(name2+'.gz', compress='gz').read() == 'test data2'
    os.unlink(name2+'.gz')
def zexists(fname):
    """Test for existance of filename or compressed versions of it.

    If fname exists, return it.  If both fname.gz and fname.bz2
    exists, raise an error.  If just one of them exists, return that
    filename.

    This can be used as a boolean function, with True indicating that
    the given fname will be openable by zopen."""
    if os.path.exists(fname):
        return fname
    gz_exists  = os.path.exists(fname+'.gz')
    bz2_exists = os.path.exists(fname+'.bz2')
    if gz_exists and bz2_exists:
        raise ValueError("Both %s.gz and %s.bz2 exist"%(fname, fname))
    if gz_exists:   return fname+'.gz'
    if bz2_exists:  return fname+'.bz2'

