
import networkx

def assert_isomorphic(g1, g2, msg=None):
    """Assertion function for networkx isomorphism"""
    if not networkx.is_isomorphic(g1, g2):
        msg_ = ["%r and %r not isomorphic"%(g1, g2),
                #"A: %s"%networkx.to_dict_of_lists(g1),
                #"B: %s"%networkx.to_dict_of_lists(g2),
                ]
        n1 = set(g1.nodes_iter())
        n2 = set(g2.nodes_iter())
        if n1 != n2:
            msg_.append("Nodes in 1 only: %s"%(n1-n2))
            msg_.append("Nodes in 2 only: %s"%(n2-n1))
        e1 = set(frozenset((a,b)) for a,b in g1.edges_iter())
        e2 = set(frozenset((a,b)) for a,b in g2.edges_iter())
        if e1 != e2:
            msg_.append("Edges in 1 only: %s"%(' '.join('(%s,%s)'%(a,b) for a,b in e1-e2)))
            msg_.append("Edges in 2 only: %s"%(' '.join('(%s,%s)'%(a,b) for a,b in e2-e1)))
        if msg: msg_.insert(0, msg)
        raise AssertionError('\n'.join(msg_))

