# Richard Darst, August 2011

import os
import sys

def getglobals():
    return {'fast':True, }

if not os.access('tests-output/', os.F_OK):
    os.mkdir('tests-output/')

execfile('tests/testmodels.py', getglobals())
execfile('tests/testmeasures.py', getglobals())

execfile('tests/testoverlap.py', getglobals())
execfile('tests/testpolopatribes.py', getglobals())
execfile('tests/testsubgraph.py', getglobals())
execfile('tests/testtiming.py', getglobals())

if 'all' in sys.argv:
    execfile('tests/test256node.py', getglobals())
    execfile('tests/testamorphous.py', getglobals())
    execfile('tests/testdolphins.py', getglobals())
    execfile('tests/testfractalsquares.py', getglobals())
    execfile('tests/testgraphcolor.py', getglobals())
    execfile('tests/testkarateclub.py', getglobals())
    execfile('tests/testlattice.py', getglobals())


