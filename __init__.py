
# Networkx 2.0 removes nodes_iter and edges_iter to use the python3
# style of nodes and edges directly returning iterators.  This means
# compatibility issues.  This is a rather bad workaround for the
# problem.
import networkx
if not hasattr(networkx.Graph, 'nodes_iter'):
    networkx.Graph.nodes_iter = networkx.Graph.nodes
    networkx.Graph.edges_iter = networkx.Graph.edges
del networkx
