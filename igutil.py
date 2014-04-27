

import igraph

def edgelist_stats(fname):
    g = igraph.Graph.Read_Edgelist(fname, directed=False)

    results = [ ]

    results.append("Nodes: %d"%len(g.vs))
    results.append("Edges: %d"%len(g.es))
    results.append("Density: %f"%g.density(loops=False))

    results.append("Mean-Degree: %f"%(2*len(g.es)/float(len(g.vs))))
    results.append("Degree-Assortitivity-Coefficient: %f"%g.assortativity_degree(directed=False))
    #results.append("Triangles: %f"%)
    results.append("Transitivity: %f"%g.transitivity_undirected(mode='nan'))
    results.append("Average-Clustering-Coefficient: %f"%
                   g.transitivity_avglocal_undirected(mode='zero'))
    results.append("Average-Clustering-Coefficient-Nonzero: %f"%
                   g.transitivity_avglocal_undirected(mode='nan'))
    #results.append(": %f"%)
    #results.append(": %f"%)
    #results.append(": %f"%)
    #results.append(": %f"%)

    return results
