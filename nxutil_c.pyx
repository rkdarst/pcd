import networkx as nx

def subgraph_density_iter(g, cmtys):
    pass

def n_edges_subgraph(adj, nodes):
    cdef int n_edges_subgraph = 0

    for n in nodes:
        for nbr in adj[n]:
            if nbr in nodes and nbr != n:
                n_edges_subgraph += 1
    return n_edges_subgraph


def n_edges_cmtys_iter(g, cmtys):
    cdef int n_edges_subgraph = 0
    adj = g.adj

    for cnodes in cmtys.itervalues():
        for n in cnodes:
            for nbr in adj[n]:
                if nbr in cnodes and nbr != n:
                    n_edges_subgraph += 1
        n_edges_subgraph /= 2
        yield n_edges_subgraph


def parse_line(it):
    pass

def read_edgelist_int(fname, create_using=nx.Graph):

    g = create_using()
    for line in open(fname):
        p = line.find('#')
        if p:
            line = line[:p]
        line = line.strip()
        if not line: continue
        s = line.split()
        u, v = int(s[0]), int(s[1])
        if len(s) == 2:
            undir_add_edge(g, u, v)
        if len(s) > 2:
            weight = float(s[2])
            undir_add_edge(g, u, v, weight=weight)
    return g

def undir_add_edge(g, u, v, weight=None):
    adj = g.adj
    node = g.node
    if u not in node:
        adj[u] = {}
        node[u] = {}
    if v not in node:
        adj[v] = {}
        node[v] = {}
    datadict = adj[u].get(v,{})
    # update attributes
    if weight:
        datadict['weight'] = weight
    adj[u][v] = datadict
    adj[v][u] = datadict
