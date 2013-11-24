

from pcd.support import powerlaw


xmin=10
alpha = -2.5
pl = powerlaw.PowerLaw(-2.5, xmin=xmin, xmax=10000)

print pl.mean()

n = 10000
xs = list(pl.rvs(n))

# mangele the power  law.
n_low = sum(1 for x in xs if x<=2*xmin)
import random
xs.extend(random.uniform(0, xmin) for i in range(n_low)  )

import plfit

pl = plfit.plfit(xs,
                 quiet=False)
print pl._alpha,  pl._xmin
#print pl.plfit()
