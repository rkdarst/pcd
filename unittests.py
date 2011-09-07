# Richard Darst, August 2011
import sys

execfile('tests/testmodels.py')
execfile('tests/testmeasures.py')
execfile('tests/testsubgraph.py')

if 'all' in sys.argv:
    execfile('tests/testpolopatribes.py')
    execfile('tests/test256node.py')
    execfile('tests/testdolphins.py')
    execfile('tests/testfractalsquares.py')
    execfile('tests/testkarateclub.py')
    execfile('tests/testlattice.py')
    execfile('tests/testtiming.py')


