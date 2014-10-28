
import os
from os.path import join, dirname

import networkx
import pcd.ioutil


def test():
    datadir = join(dirname(pcd.ioutil.__file__), 'data')
    fname_gml = join(datadir, 'karate.gml')
    fname_edgelist = join(datadir, 'test_edgelist')
    fname_pajek = join(datadir, 'test_pajek')

    assert pcd.ioutil._test_gml(open(fname_gml).read(512))
    assert pcd.ioutil._get_reader('gml') == networkx.read_gml
    assert pcd.ioutil._test_edgelist(open(fname_edgelist).read(512))
    assert pcd.ioutil._get_reader('edgelist') == networkx.read_edgelist
    assert pcd.ioutil._test_pajek(open(fname_pajek).read(512))
    assert pcd.ioutil._get_reader('pajek') == networkx.read_pajek

    data_gexf = '''<?xml version="1.0" encoding="UTF-8"?>
    <gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">
    <meta lastmodifieddate="2009-03-20">
    <creator>Gexf.net</creator>
    <description>A hello world! file</description>'''
    assert pcd.ioutil._test_gexf(data_gexf)

    data_graphml = '''<?xml version="1.0" encoding="UTF-8"?>
    <graphml xmlns="http://graphml.graphdrawing.org/xmlns"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
    http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
    <graph id="G" edgedefault="undirected">
    <node id="n0"/>
    '''
    assert pcd.ioutil._test_graphml(data_graphml)


    assert isinstance(pcd.ioutil.read_any(fname_gml), networkx.Graph)
    assert isinstance(pcd.ioutil.read_any('gml:'+fname_gml), networkx.Graph)

    pcd.ioutil._test_zopen()
