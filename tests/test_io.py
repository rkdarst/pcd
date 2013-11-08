
import os
from os.path import join, dirname

import networkx
import pcd.ioutil

datadir = join(dirname(pcd.ioutil.__file__), 'data')
fname_gml = join(datadir, 'karate.gml')
fname_edgelist = join(datadir, 'test_edgelist')
fname_pajek = join(datadir, 'test_pajek')

assert pcd.ioutil._test_gml(fname_gml)
assert pcd.ioutil._get_reader('gml') == networkx.read_gml
assert pcd.ioutil._test_edgelist(fname_edgelist)
assert pcd.ioutil._get_reader('edgelist') == networkx.read_edgelist
assert pcd.ioutil._test_pajek(fname_pajek)
assert pcd.ioutil._get_reader('pajek') == networkx.read_pajek

assert isinstance(pcd.ioutil.read_any(fname_gml), networkx.Graph)
assert isinstance(pcd.ioutil.read_any('gml:'+fname_gml), networkx.Graph)

pcd.ioutil._test_zopen()
