
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
    raise unittest.SkipTest
    func = staticmethod(log_bin)
    width = staticmethod(log_bin_width)
    binargs = dict(minlog=0)
    range = numpy.arange(.1, 100, .001)

class test_logbin2(_test_bin):
    raise unittest.SkipTest
    func = staticmethod(log_bin)
    width = staticmethod(log_bin_width)
    binargs = dict(minlog=10)
    range = numpy.arange(10, 100, .001)

