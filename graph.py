import pickle
import networkx as nx
import numpy as np
import numpy.linalg
from scipy.sparse.csgraph import laplacian

vals = []
pkls = ["1", "2", "6", "10", "13", "16"]

for i in range(0, len(pkls)):
	fileObject = open("cuprotein_graph_6.pkl",'rb') 
	b = pickle.load(fileObject)
	these_vals = list(b.values())	
	vals += these_vals



graphs = list()

# returns list of eigenvalues of graph g
# organized greatest to least
def get_eigenvalues(g):
	adj = g["adj_contact_map"]
	for i in range(0, len(adj)):
		adj[i][i] = 0
	
	lap = laplacian(adj)
	e = numpy.linalg.eigvals(lap)
	e.sort()
	e = e[::-1]
	return e

def get_nx_graph(g):
	G = nx.Graph()
	
	for i in range(0, len(g["aa"])):
		G.add_node(i)
		
	for i in range(0, G.number_of_nodes()):
		for j in range(0, G.number_of_nodes()):
			if (not (i==j)) and g["adj_contact_map"][i][j] == 1:
				edge_length = 1 / g["adj_inverse_dcalpha"][i][j]
				G.add_edge(i, j, weight = edge_length)
				G.add_edge(j, i, weight = edge_length)
				
	print(G.number_of_nodes())
	print(G.number_of_edges())
				
	return G
				
				


def ISMAGS_method(g1, g2):
	gr1 = get_nx_graph(g1)
	gr2 = get_nx_graph(g2)
	ismags = nx.isomorphism.ISMAGS(gr1, gr2)
	lcs = list(ismags.largest_common_subgraph(symmetry=False))
	print(len(lcs))
	print(len(lcs[0]))
	

# performs eigenvalue method on graphs g1, g2
def eigenvalue_method(g1, g2):
	eig1 = get_eigenvalues(g1)
	sum1 = sum(eig1)
	k1 = len(eig1)
	for i in range(0, len(eig1)):
		if sum(eig1[:i]) >= 0.9 * sum1:
			k1 = i
			break
	
	eig2 = get_eigenvalues(g2)
	sum2 = sum(eig2)
	k2 = len(eig2)
	for i in range(0, len(eig2)):
		if sum(eig2[:i]) >= 0.9 * sum2:
			k2 = i
			break
			
	k = min(k1, k2)
	dist = 0
	for i in range(0, k):
		dist += pow(eig1[i] - eig2[i], 2)
		
	return dist


for i in range(0, len(vals)):
	for j in range(0, len(vals)):
		good = False
		if not (i == j):
			good = True
			if len(vals[i]["aa"]) == len(vals[j]["aa"]):
				for k in range(0, len(vals[i]["aa"])):
					if not np.array_equal(vals[i]["aa"][k], vals[j]["aa"][k]):
						good = False
						break
			else:
				good = False
						
		if good:
			print(str(i)+" "+str(j))




"""
for i in range(0, 1):
	for j in range(0, 1):
		res = ISMAGS_method(vals[75], vals[88])
"""

"""
for i in range(0, len(vals)):
	for j in range(0, len(vals)):
		eig_res = eigenvalue_method(vals[i], vals[j])
		print(eig_res)
"""







"""
for graph_index in range(0, 2):

	val = vals[graph_index]
	print(len(val["aa"]))
	print(len(val["aa"][0]))
	print(len(val["adj_contact_map"]))
	print(len(val["adj_contact_map"][0]))
	print(len(val["adj_inverse_dcalpha"]))
	print(len(val["adj_inverse_dcalpha"][0]))
	print()
	
	
	data = vals[graph_index]
	G = nx.Graph()
	
	for i in range(0, len(data["aa"])):
		G.add_node(i)
		
	for i in range(0, G.number_of_nodes()):
		for j in range(0, G.number_of_nodes()):
			if (not (i==j)) and data["adj_contact_map"][i][j] == 1:
				edge_length = 1 / data["adj_inverse_dcalpha"][i][j]
				G.add_edge(i, j, weight = edge_length)
	
	graphs.append(G)
"""