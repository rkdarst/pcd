from nose.tools import assert_almost_equal
import unittest

from pcd.cmty import Communities

import cmtycmp

cmtys_one = Communities({0: list(range(100))})
# nodes 0-99 randomly distributed among 4 communities
cmtys_random_1A = Communities(
    {0: set([6, 15, 17, 20, 21, 22, 26, 28, 30, 33, 34, 36, 44, 47, 53,
             55, 58, 60, 67, 68, 79, 81, 83, 85, 95, 97]),
     1: set([64, 1, 2, 35, 98, 70, 96, 41, 42, 75, 12, 13, 14, 45, 88,
             50, 84, 78, 56, 90, 61]),
     2: set([0, 18, 19, 23, 24, 27, 29, 31, 37, 38, 39, 46, 48, 51, 52,
             62, 63, 65, 66, 69, 71, 72, 73, 74, 76, 77, 80, 86, 87, 89,
             94, 99]),
     3: set([32, 3, 4, 5, 7, 8, 9, 10, 11, 91, 59, 16, 40, 82, 43, 57, 25,
             49, 54, 92, 93])}
    )
cmtys_random_2A = Communities(
    {0: set([84, 33, 2, 99, 36, 5, 51, 39, 73, 74, 83, 17, 18, 19, 20,
             85, 23, 25, 58, 37]),
     1: set([0, 3, 4, 11, 15, 24, 26, 28, 42, 46, 48, 52, 55, 68, 80,
             81, 86, 88, 91, 92, 95, 96]),
     2: set([6, 7, 8, 14, 16, 22, 27, 30, 31, 34, 35, 41, 43, 44, 45, 47,
             49, 53, 54, 56, 57, 60, 61, 63, 65, 66, 71, 75, 78, 79, 82,
             89, 90, 97, 98]),
     3: set([1, 9, 10, 12, 13, 21, 29, 32, 38, 40, 50, 59, 62, 64, 67,
             69, 70, 72, 76, 77, 87, 93, 94])}
    )
cmtys_random_3A = Communities(
    {0: set([68, 73, 98, 20, 43, 17, 23, 27, 67, 28, 53, 89, 78, 95, 
             42, 62, 85, 33, 86, 46, 44, 4, 63, 90, 34, 47, 76, 31, 
             55, 30, 25, 56, 88, 72, 8, 35, 40, 36, 77, 38, 80, 14, 
             96, 1, 99, 3, 82, 51, 87, 41]),
     1: set([74, 50, 39, 83, 24, 37, 13, 21, 84, 18, 54, 16, 64, 9, 
             81, 69, 91, 93, 12, 57, 15, 2, 48, 52, 66, 60, 61, 58, 5,
             10, 11, 45, 71, 7, 22, 70, 65, 75, 92, 59, 79, 94, 6, 29,
             26, 49, 32, 97, 0, 19])}
    )


skipped_tests = [
    # LF code returns 0.0 instead of 1.0 for identical partitions.  I
    # think this is a bug.
    ('nmi_LFK_LF', id(cmtys_one), id(cmtys_one)),
    # LF also returns 0.0 when any one partition has only one
    # community.  Is this expected behavior or a bug?
    ('nmi_LFK_LF', id(cmtys_one), id(cmtys_random_1A)),
    ('nmi_LFK_LF', id(cmtys_random_1A), id(cmtys_one)),
    ('nmi_LFK_LF', id(cmtys_random_3A), id(cmtys_one)),
    ('nmi_LFK_LF', id(cmtys_one), id(cmtys_random_3A)),

    # igraph might might use natural log instead of log2.  In
    # cmtycmp.py, igraph functions are changed to return in bits
    # instead.  FIXME: what is the proper unit to return in?
    #('vi_igraph',         id(cmtys_one), id(cmtys_random_1A)),
    #('vi_igraph',         id(cmtys_random_1A), id(cmtys_random_2A)),

    # This returns nan.
    ('nmi_LFK_pcdpy', id(cmtys_one), id(cmtys_random_1A)),
    ('nmi_LFK_pcdpy', id(cmtys_random_1A), id(cmtys_one)),
    
    # returns >1 while other LFK measures return something somewhat intelligent
    ('nmi_LFK_pcdpy', id(cmtys_one), id(cmtys_random_3A)),
    ('nmi_LFK_pcdpy', id(cmtys_random_3A), id(cmtys_one)),

    # This is basically correct, but only to 2 places, not 6.  Or it
    # could be the _LF implementation is off.  FIXME.
    ('nmi_LFK_pcd', id(cmtys_random_1A), id(cmtys_random_2A)),
    ('nmi_LFK_pcd', id(cmtys_random_2A), id(cmtys_random_1A)),
    ('nmi_LFK_pcdpy', id(cmtys_random_1A), id(cmtys_random_2A)),
    ('nmi_LFK_pcdpy', id(cmtys_random_2A), id(cmtys_random_1A)),
    ('nmi_LFK_pcdpy', id(cmtys_random_2A), id(cmtys_random_3A)),
    
    # Returns nan, even though it should return 1.0
    ('adjusted_rand_igraph', id(cmtys_one), id(cmtys_one)),
    ]

def _do_test_same(func, cmtys):
    x = func(cmtys, cmtys)
    assert_almost_equal(x, 1.0)

def test_same():
    for name, implementations in cmtycmp.measures.iteritems():
        if name in ('vi', 'mutual_information', 'vi_norm'):
            continue
        for impl in implementations:
            if impl in ('nmi_LFK_LF', 'adjusted_rand_igraph',
                        'omega_index_python','gamma_coeff_python'):
                continue
            yield _do_test_same, getattr(cmtycmp, impl), cmtys_one
            yield _do_test_same, getattr(cmtycmp, impl), cmtys_random_1A

def _do_test_one(measure, cmtys=None):
    implementations = cmtycmp.measures[measure]
    funcs = [ getattr(cmtycmp, name) for name in implementations ]

    #cmtys = (cmtys_one, cmtys_one)
    #cmtys = (cmtys_random_1A, cmtys_random_1A)

    for i, name1 in enumerate(implementations):
        if (name1, id(cmtys[0]), id(cmtys[1])) in skipped_tests:
            continue
            #raise unittest.SkipTest
        func1 = getattr(cmtycmp, name1)
        v1 = func1(*cmtys)
        for j, name2 in enumerate(implementations[i+1:]):
            if (name2, id(cmtys[0]), id(cmtys[1])) in skipped_tests:
                continue
            func2 = getattr(cmtycmp, name2)
            v2 = func2(*cmtys)
            print measure, name1, name2, v1, v2
            print
            assert_almost_equal(v1, v2, places=6,
                      msg="%s %s != %s %s"%(name1, v1,
                                            name2, v2))
#def _assert_one(msg, v1, v2)

def test_one():
    for name in cmtycmp.measures:
        yield _do_test_one, name, (cmtys_one, cmtys_one)

def test_random1():
    for name in cmtycmp.measures:
        yield _do_test_one, name, (cmtys_random_1A, cmtys_random_1A)

def test_one_random1():
    for name in cmtycmp.measures:
        yield _do_test_one, name, (cmtys_one, cmtys_random_1A)

def test_one_random3():
    for name in cmtycmp.measures:
        yield _do_test_one, name, (cmtys_random_3A, cmtys_one)

def test_random2_random3():
    for name in cmtycmp.measures:
        yield _do_test_one, name, (cmtys_random_2A, cmtys_random_3A)

def test_distance1():
    A = Communities({0: set([1,2,3,4,5,6]),1: set([7,8,9])})
    B = Communities({0: set([1,2,4,5,7,8]),1: set([3,6]),2: set([9])})
    t = cmtycmp.distance_division_python(A,B)
    print
    assert_almost_equal(t, 2.0/3, places=6, msg="testdistance1 %f"%(t))
    
def test_distance2():
    A = Communities({0: set([1,2,3,4,5,6]),1: set([7,8,9])})
    B = Communities({0: set([1,2,4,5,7,8]),1: set([3,6]),2: set([9])})
    t = cmtycmp.distance_moved_python(A,B)
    print
    assert_almost_equal(t, 4.0/9, places=6, msg="testdistance2 %f"%(t))   
    #cmtys = (cmtys_one, cmtys_one)
    #cmtys = (cmtys_random_1A, cmtys_random_1A)
    
def test_distance3():
    t = cmtycmp.distance_moved_python(cmtys_one,cmtys_random_3A)
    print
    assert_almost_equal(t, 0.5, places=6, msg="testdistance3 %f"%(t))      
