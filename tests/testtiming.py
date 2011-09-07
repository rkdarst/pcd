# Richadr Darst, July 2011

import timeit

def setup():
    G = models.random_graph(size=5)
    G.check()


def make_setup(size=5, **kwargs):
    lines = [ ]
    lines.append('from graphs import random_graph')
    lines.append('G = random_graph(size=%s)'%size)
    return '\n'.join(lines)


def make_stmt(**kwargs):
    lines = [ ]
    lines.append('G.minimize(gamma=1.0)')
    return '\n'.join(lines)

results = { }
#for size in (5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70):
for size in (100, ):
    kwargs = {'size':size}
    T = timeit.Timer(stmt=make_stmt(**kwargs),
                     setup=make_setup(**kwargs), )
    times = T.repeat(repeat=3, number=1)
    results[size] = min(times)

for size, t in sorted(results.items()):
    print size, t
