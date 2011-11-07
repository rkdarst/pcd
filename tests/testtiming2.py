# Richard Darst, October 2011

import timeit

class Test(object):
    defaults = {'graph': 'random_graph',
                'size': 50,
                'trials':5,
                'gamma': 1,
                'gargs': dict(size=50),
                'verbosity':-2}
    def __init__(self,
                 setup=('from graphs import %(graph)s',
                        'G = %(graph)s(**%(gargs)r)',
                        'G.verbosity = %(verbosity)s'),
                 command=(('#print 11',
                           '#G.minimize(gamma=%(gamma)f)',
                           'G.trials(gamma=%(gamma)f, trials=%(trials)d)',
                           '#print 22',)),
                 name=None,
                 n=3,
                 **kwargs):
        #setup = [ ]
        #sa = setup.append
        #sa('from graphs import %(graph)s')
        #if gargs in kwargs:
        #    sa('G = %(graph)s(**%(gargs)')
        #else:
        #    sa('G = %(graph)s()')
        #sa('G.verbosity = -2')
        #command = [ ]
        #ca = command.append
        #ca('G.minimize(gamma=%(gamma)f)')


        defaults = self.defaults.copy()
        defaults.update(kwargs)
        self.kwargs = defaults
        self.setup = '\n'.join(setup)
        self.command = '\n'.join(command)
        self.name = name
        self.n = n

    def run(self):
        #print self.setup, self.command, self.kwargs
        T = timeit.Timer(setup=self.setup%self.kwargs,
                         stmt=self.command%self.kwargs)
        try:
            times = T.repeat(self.n, number=1)
        except:
            T.print_exc()
            import sys ; sys.exit()
        time = min(times)
        num = self.n * self.kwargs['trials']
        print ("%-15s %9.2f/real s  (%3.2e) (%3.1fs here)"%
               (self.name, num/time, time/num, sum(times)))



runs = (
    #Test(name="rg50*2"),
    Test(name="dolphins", graph='dolphins_G', gargs={}),
    Test(name="n5760", graph='bss2d_n5760_T040', gargs={}, n=1),
    )

for r in runs:
    r.run()
