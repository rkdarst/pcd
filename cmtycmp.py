log2 = lambda x: log(x, 2)



def mutual_information(cmtys1, cmtys2):
    MI = 0.0
    N = len(cmtys1.nodes)
    assert N == len(cmtys2.nodes)
    for c1name, c1nodes in cmtys1.iteritems():
        for c2name, c2nodes in cmtys2.iteritems():
            n1 = len(c1nodes)
            n2 = len(c2nodes)
            n_shared = len(set(c1nodes) & set(c2nodes))
            if n_shared == 0:
                continue
            MI += (n_shared/float(N)) * log2(n_shared*N/float(n1*n2))
    return MI
def entropy(cmtys):
    H = 0.0
    N = float(len(cmtys.nodes))
    for cnodes in cmtys.itervalues():
        n = len(cnodes)
        if n == 0: continue
        H += n/N * log2(n/N)
    return -H
def In(cmtys1, cmtys2):
    I = mutual_information(cmtys1, cmtys2)
    Hs = (entropy(cmtys1) + entropy(cmtys2))
    if Hs == 0:
        return float('nan')
    In = 2.*I / Hs
    return In


def NMI_LF(cmtys1, cmtys2, check=True, use_existing=False):
    """Compute NMI using the overlap-including definition.

    This uses the external code 'mutual3/mutual' to calculate the
    overlap-using In.

    If the communities object is a pcd.cmty.CommunityFile, has a fname
    attribute, and it exists and doesn't end in .gz or .bz2, and the
    argument use_existing is true, then this file will be used for the
    input to the NMI code.  The NMI code is very dumb, and does not
    remove comments, and uses only floating-point or integer node IDs,
    and does not give any type of error if these conditions are not
    met.  Thus, only enable use_existing if you can vouch for the.

    check: bool, default True
        If true, check that all node names are either representable as
        integers or floating point numbers.  This can either be python
        integers or floats, or strings of those.
    use_existing: bool, default False
        If true, and community object has a '.fname' attribute, and
        that file exists, use that as the input filename instead of
        re-writing the file.  WARNING: if you use this, you must
        ensure that there are no comments within the file, or anything
        else which make cause the program to break.
        """
    from pcd.support.algorithms import _get_file
    from pcd.cmty import CommunityFile
    binary = _get_file('mutual3/mutual')

    def _is_float_str(x):
        """Is this a string representation of a float?"""
        try:
            float(x)
            return True
        except ValueError:
            return False
    # check that community files are valid for the program:
    if check:
        for nodes in cmtys1.itervalues():
            assert all(isinstance(x, (int, float)) or _is_float_str(x) for x in nodes )
        for nodes in cmtys2.itervalues():
            assert all(isinstance(x, (int, float)) or _is_float_str(x) for x in nodes )

    # We must use os.path.abspath *outside* of the tmpdir_context, or
    # else the absolute path will be wrong.
    args = [ binary ]
    if (use_existing
        and isinstance(cmtys1, CommunityFile)
        and hasattr(cmtys1, 'fname')
        and os.path.exists(cmtys1.fname)
        and not (cmtys1.fname.endswith('.bz2')
                 or cmtys1.fname.endswith('.gz'))):
        args.append(os.path.abspath(cmtys1.fname))
    else:
        args.append(None)

    if (use_existing
        and isinstance(cmtys1, CommunityFile)
        and hasattr(cmtys2, 'fname')
        and os.path.exists(cmtys2.fname)
        and not (cmtys2.fname.endswith('.bz2')
                 or cmtys2.fname.endswith('.gz'))):
        args.append(os.path.abspath(cmtys2.fname))
    else:
        args.append(None)

    with tmpdir_context(chdir=True, dir='.', prefix='tmp-nmi-'):
        # Write community files, if they do not already exist.  These
        # must be written inside of the tmpdir_context context because
        # only in here does it know the right
        if args[1] is None:
            cmtys1.write_clusters('cmtys1.txt', raw=True)
            args[1] = 'cmtys1.txt'
        #else:
        #    # Symlink it instead of run in-place (program can fail if
        #    # filenames are weird!)
        #    os.symlink(args[1], 'cmtys1.txt')
        #    args[1] = 'cmtys1.txt'
        if args[2] is None:
            cmtys2.write_clusters('cmtys2.txt', raw=True)
            args[2] = 'cmtys2.txt'
        #else:
        #    os.symlink(args[2], 'cmtys2.txt')
        #    args[1] = 'cmtys2.txt'
        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        ret = p.wait()
        stdout = p.stdout.read()
        if ret != 0:
            print stdout
            raise RuntimeError("The program '%s' returned non-zero: %s"%(
                args[0], ret))
        nmi = float(stdout.split(':', 1)[1])
    return nmi
ovIn_LF = NMI_LF
