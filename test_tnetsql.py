from nose.tools import *

from tnetsql import TNet

T = TNet()
T.add_edges([
    # a, b, t
    (0, 1,  0),
    (0, 2,  1),
    (1, 2,  1),
    (2, 3,  2),
    (1, 3,  3),
    (0, 4,  3),
    ])
T.dump()

def test_has_edge():
    assert_true( T.has_edge(0, 1))
    assert_false(T.has_edge(0, 3))

def test_successors():
    assert_equal( list(T.successors_after(0, 2, 100)),
                   [(4, )])
    #assert_false(T.has_edge(0, 3))



