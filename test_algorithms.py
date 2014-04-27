
import pcd.support.algorithms as algs
import pcd.graphs

from nose.plugins.attrib import attr

@attr('slow')
def test_algorithms():
    g = pcd.graphs.karate_club()

    exclude_methods = set(('CDMethod',
                           'BeliefPropogationZ',
                           'BeliefPropogationq',
                           'PCDmodSA',   # works but slow.
                           ))


    methods = sorted(
                name for name, cda in vars(algs).iteritems()
                if isinstance(cda, type)
                and issubclass(cda, algs.CDMethod)
                and name[0]!='_')

    methods = [x for x in methods if
               x not in exclude_methods
               and 'APM' not in x]
    method_args = {
        'BeliefPropogationZ':dict(q=2),
        }
    default_args = dict(verbosity=0)

    for method in methods:
        print method
        cda = getattr(algs, method)
        args = default_args.copy()
        args.update(method_args.get(method, {}))
        cda(g, **args)
