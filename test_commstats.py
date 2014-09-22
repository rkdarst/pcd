
import unittest
from nose.tools import *

import numpy
from numpy.testing.utils import assert_allclose

from pcd.commstats import counts, norm
from pcd.commstats import lin_bin, lin_bin_width, log_bin, log_bin_width


class _test_bin(unittest.TestCase):
    range = numpy.arange(0, 100, .01)
    binargs = dict()
    def test_run(self):
        values = self.range
        binned = [self.func(x, **self.binargs) for x in values]
        count = counts(binned)
        #print count
        x, v = zip(*zip(*count))
        v = numpy.array(v)
        print x
        print v
        v = v / [self.width(_, **self.binargs) for _ in v]
        print v
        # exclude boundraies
        x = x[1:-1]
        v = v[1:-1]
        v = numpy.divide(v, max(v))
        print v
        assert_allclose(v, 1, rtol=.02)

class test_linbin(_test_bin):
    func = staticmethod(lin_bin)
    width = staticmethod(lin_bin_width)
    range = numpy.arange(0, 100, .01)

class test_logbin(_test_bin):
    #raise unittest.SkipTest
    func = staticmethod(log_bin)
    width = staticmethod(log_bin_width)
    binargs = dict(minlog=0)
    range = numpy.arange(.1, 100, .001)

class test_logbin2(_test_bin):
    #raise unittest.SkipTest
    func = staticmethod(log_bin)
    width = staticmethod(log_bin_width)
    binargs = dict(minlog=10)
    range = numpy.arange(10, 100, .001)

from nose.tools import *
def aeq(x, y):
    "approximately equal"
    return abs(x-y)/max(abs(x), abs(y)) < 1e-6
def test_logbin_float():
    # Everything less than 10 should be integers only.
    r = range(-10, 10+1)
    b = [log_bin(x, minlog=10) for x in r]
    [ assert_almost_equal(x, y) for x,y in zip(r, b) ]

    # Check that all bin widths are one
    r = range(-10, 9+1)  # 10 is different
    w = [log_bin_width(x, minlog=10) for x in r]
    #print r, w
    [ assert_almost_equal(1, y) for x,y in zip(r, w) ]

    # Long bin of 11 is no longer integer.
    assert_not_almost_equal(log_bin(11, minlog=10), 11)

    print 9.9, log_bin(9.9,minlog=10), log_bin_width(9.9,minlog=10)
    print 10 , log_bin(10, minlog=10), log_bin_width(10, minlog=10)
    print 11 , log_bin(11, minlog=10), log_bin_width(11, minlog=10)
    print 12 , log_bin(12, minlog=10), log_bin_width(12, minlog=10)

def test_logbin_int():
    r = range(0, 1000)
    counts = { }
    for x in r:
        v = log_bin(x, minlog=10, ints=True)
        print x, v
        counts[v] = counts.get(v, 0) + 1
    for v, w in counts.iteritems():
        if v > 900: continue
        print v, w
        assert_equal(w, log_bin_width(v, minlog=10, ints=True))


    r = range(0, 1000)
    counts = { }
    for x in r:
        v = log_bin(x, minlog=1, ints=True)
        print x, v, 'xx'
        counts[v] = counts.get(v, 0) + 1
    for v, w in counts.iteritems():
        if v > 900: continue
        print v, w, 'xx'
        assert_equal(w, log_bin_width(v, minlog=1, ints=True))
