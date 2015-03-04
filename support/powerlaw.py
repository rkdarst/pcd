# Richard Darst, May 2012

import itertools
import math
from math import log, exp
import numpy
import random


class PowerLaw(object):
    """

    = 1/N * x**power
    """
    tol = 1e-12
    def __str__(self):
        return '%s(power=%s, %s)'%(self.__class__.__name__, self.power,
                ", ".join("%s=%s"%(k,v) for k,v in sorted(self.orig_args.items())
                          if v is not None))
    @staticmethod
    def integrate(x, power):
        """Integral from zero to x"""
        if power == -1:
            return log(abs(x))
        else:
            return float(pow(x, power+1)) / (power+1.)
    @staticmethod
    def inverseintegrate(x, power):
        """Inverse of integral funcition."""
        if power == -1:
            return exp(x)
        else:
            return pow(x*(power+1.),  1./(power+1.))

    def __init__(self, power, xmin=None, xmax=None, xmean=None):
        self.power = power
        self.orig_args = dict(xmin=xmin, xmax=xmax, xmean=xmean)
        if xmin == 'inf': xmin = float('-inf')
        if xmax == 'inf': xmax = float('inf')
        if xmin is not None and xmax is not None:
            # if xmin and xmax given, xmean must not be given, and
            # then use the values directly.
            if xmean is not None:
                raise ValueError("Can't specify xmean with xmin and xmax.")
        elif xmin is not None and xmean is not None:
            # if xmin and xmean given, calculate xmax (and check
            # mean==xmean below)
            xmax = self._calc_xmax(xmin, xmean)
        elif xmax is not None and xmean is not None:
            # if xmax and xmean given, calculate xmin (and check
            # mean==xmean below)
            xmin = self._calc_xmin(xmax, xmean)

        self.norm = self._integrate(xmax)-self._integrate(xmin)
        self.xmin = xmin
        self.xmax = xmax
        if xmin:
            self._lowbound = self._integrate(xmin)
        # Double check that xmean was calculated properly.
        if xmean is not None:
            assert abs(self.mean() - xmean)*.5 < self.tol

    def _calc_xmax(self, xmin, xmean):
        """Given xmin and xmean, calculate xmax."""
        self.xmin = xmin
        self.xmax = xmean+.01
        for i in itertools.count():
            self.norm = self._integrate(self.xmax)-self._integrate(self.xmin)
            newmean = self.mean()
            delta = xmean - newmean
            #print self.xmin, newmean, self.xmax, delta
            if delta < self.tol:
                break
            self.xmax += delta
            if i > 1e3:
                raise ValueError("Power law parameter calculation did not converge: "
                                 "power=%s, xmin=%s, mean=%s"%(self.power, xmin, xmean))
        return self.xmax
    def _calc_xmin(self, xmax, xmean):
        """Given xmax and xmean, calculate xmin."""
        self.xmin = xmean-.01
        self.xmax = xmax
        for i in itertools.count():
            self.norm = self._integrate(self.xmax)-self._integrate(self.xmin)
            newmean = self.mean()
            delta = - newmean + xmean
            #print self.xmin, newmean, self.xmax, delta
            if abs(delta) < self.tol:
                break
            self.xmin += delta
            if i > 1e3:
                raise ValueError("Power law parameter calculation did not converge: "
                                 "power=%s, xmean=%s, xmax=%s"%(self.power, xmean, xmax))
        return self.xmin

    def _integrate(self, x):
        return self.integrate(x, self.power)
    def _inverseintegrate(self, x):
        return self.inverseintegrate(x, self.power)

    def pdf(self, x):
        return float(pow(x, self.power))  / self.norm
    def cdf(self, x):
        if x < self.xmin: raise ValueError("x < xmin")
        if x > self.xmax: raise ValueError("x > xmax")
        return (self._integrate(x) - self._lowbound) / self.norm
    def inversecdf(self, p):
        return self._inverseintegrate(p*self.norm + self._lowbound)
    def rv(self):
        return self.inversecdf(random.uniform(0, 1))
    def rvs(self, n=1):
        for i in xrange(n):
            yield self.inversecdf(random.uniform(0, 1))
    def mean(self):
        return (self.integrate(self.xmax, self.power+1)
                - self.integrate(self.xmin, self.power+1)) / self.norm
    def moment(self, m):
        return (self.integrate(self.xmax, self.power+m)
                - self.integrate(self.xmin, self.power+m)) / self.norm
    def var(self):
        return self.moment(2) - self.mean()**2
    def std(self):
        return math.sqrt(self.var())

class PowerLawInteger(PowerLaw):
    @staticmethod
    def integrate(self, n, power):
        return (1. - power**(n+1) ) / (1. - power)
    @staticmethod
    def integrate(self, n, power):
        pass



if __name__ == "__main__":
    import sys

    import pcd.util
    if len(sys.argv) > 1:
        data = sys.argv[1:]
        data = [x.split('=', 1) for x in data]
        kwargs = dict((x, pcd.util.leval(y)) for x,y in data)
        print kwargs

        p = PowerLaw(**kwargs)
        print p
        print 'norm:', p.norm
        print 'mean:', p.mean()
        print 'std: ', p.std()
        print 'min, max:', p.xmin, p.xmax
        print '3 RVs:', list(p.rvs(3))
        #xs = list(p.rvs(10))
        #print xs
        #print numpy.power(xs, -2+1.5)
        sys.exit()

    p = PowerLaw(-2, 1, 10)
    assert p.cdf(10) == 1
    assert p.cdf(1) == 0
    print p.inversecdf(0), p.inversecdf(1)
    from fitz import interactnow
