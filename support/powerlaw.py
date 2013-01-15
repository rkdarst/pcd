# Richard Darst, May 2012

import math
import numpy
import random


class PowerLaw(object):
    """

    = 1/N * x**power
    """
    @staticmethod
    def integrate(x, power):
        """Integral from zero to x"""
        if power == -1:
            return math.log(x)
        else:
            return float(pow(x, power+1)) / (power+1)
    @staticmethod
    def inverseintegrate(x, power):
        """Inverse of integral funcition."""
        if power == -1:
            return math.exp(x)
        else:
            return (power+1) * float(pow(x, 1./(power+1)))

    def __init__(self, power, xmin=None, xmax=None):
        self.power = power
        self.xmin = xmin
        self.xmax = xmax
        if xmin is not None and xmax is not None:
            self.norm = self.integrate(xmax, power)-self.integrate(xmin, power)
        if xmin:
            self._lowbound = self.integrate(xmin, power)
    def _integrate(self, x):
        return self.integrate(x, self.power)

    def pdf(self, x):
        return float(pow(x, self.power))  / self.norm
    def cdf(self, x):
        if x < self.xmin: raise ValueError("x < xmin")
        if x > self.xmax: raise ValueError("x > xmax")
        return (self._integrate(x) - self._lowbound)  / self.norm
    def inversecdf(self, p):
        return self.inverseintegrate(p*self.norm + self._lowbound, self.power)
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
        print 'norm:', p.norm
        print 'mean:', p.mean()
        print 'std: ', p.std()
        #xs = list(p.rvs(10))
        #print xs
        #print numpy.power(xs, -2+1.5)
        sys.exit()

    p = PowerLaw(-2, 1, 10)
    assert p.cdf(10) == 1
    assert p.cdf(1) == 0
    print p.inversecdf(0), p.inversecdf(1)
    from fitz import interactnow
