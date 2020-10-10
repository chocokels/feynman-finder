import networkx as nx
from FeynmanFinder import *

# Create graph
graph = nx.DiGraph()
graph.add_nodes_from(['0','i','j','2'])
graph.add_edge('0','i')
graph.add_edge('0','j')
graph.add_edge('i','2')
graph.add_edge('j','2')

eta = [-1,1,1,-1]

# Find diagrams
feyn = FeynmanFinderFunc('0','0',graph,eta)

# Print diagrams
PrintFeynman(feyn)
