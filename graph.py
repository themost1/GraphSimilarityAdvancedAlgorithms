import pickle
import networkx as nx
import numpy as np
import numpy.linalg
from scipy.sparse.csgraph import laplacian
import glob
import math
import random
import copy
import time


vals = []
pkls = glob.glob("./*.pkl")

# Load our files of protein graphs
max_pkls = 1
num_pkls_load = min(max_pkls, len(pkls))
for i in range(0, num_pkls_load):
	fileObject = open(pkls[i],'rb') 
	b = pickle.load(fileObject)
	these_vals = list(b.values())	
	vals += these_vals
	
print("loaded: " + str(num_pkls_load))


curr_time = lambda: int(round(time.time() * 1000))

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

# used with NetworkX library -- unused, but here just in case
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
				
			
#NetworkX function that is too costly to use
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

# performs Euclidean method on graphs g1, g2
def euclidean_method(g1, g2):
	sum = 0
	minlen = min(len(g1["aa"]), len(g2["aa"]))
	for i in range(0, minlen):
		minisum = 0
		for j in range(0, 3):
			diff = g1["calpha_coord"][i][j] - g2["calpha_coord"][i][j]
			diff = diff*diff
			minisum += diff
		sum += minisum
	sum /= minlen
	sum = math.sqrt(sum)
	return sum

# performs RMSD on graphs g1, g2
def rmsd(g1, g2):
	dcalpha1 = g1['adj_inverse_dcalpha']
	dcalpha1 = np.reciprocal(dcalpha1)
	dcalpha1 = np.eye(np.shape(dcalpha1)[0]) - dcalpha1
	dcalpha2 = g2['adj_inverse_dcalpha']
	dcalpha2 = np.reciprocal(dcalpha2)
	dcalpha2 = np.eye(np.shape(dcalpha2)[0]) - dcalpha2
	total = np.shape(dcalpha2)[0]**2
	tmp = dcalpha2 - dcalpha1
	tmp = tmp**2
	squared = np.sum(tmp)
	mean_square = squared/total
	rmsd_i = np.sqrt(mean_square)
	return rmsd_i

# performs iterative similarity score method on graphs g1, g2
def sim_score_method(g1, g2):
	cutoff = 0.0001
	d = []
	row = []
	adj1 = g1["adj_contact_map"]
	adj2 = g2["adj_contact_map"]
	s = np.identity(len(adj1))
	for i in range(0, len(adj1)):
		row.append(0.5)
	for i in range(0, len(adj1)):
		d.append(row)
		
	while True:
		prod1 = np.matmul(adj1, s)
		prod1 = np.matmul(prod1, adj2)
		prod2 = np.matmul(adj1.transpose(), s)
		prod2 = np.matmul(prod2, adj2.transpose())
		sum = prod1 + prod2 + d
		norm = np.linalg.norm(sum)
		new_s = sum / norm
		
		diff = 0
		for i in range(0, len(s)):
			for j in range(0, len(s)):
				diff += abs(s[i][j] - new_s[i][j])
		s = new_s
		if (diff < cutoff):
			break

	return s.sum()


""" Quadratic time approx of GED using Hausdorff matching,
chosen just because it was recent,
and there is some testing of it with molecule graphs,
and it's quadratic not cubic
"""
"""
u,v are edges so they should contain some edge information
"""
def lower_vertex(u,v,g1,g2,cost):
	adj1 = g1["adj_contact_map"]
	adj2 = g2["adj_contact_map"]
	p = []
	q = []

	for i in range(0,len(adj1[u])):
		if adj1[u][i] == 1:
			p.append(i)
	for i in range(0,len(adj2[v])):
		if adj2[v][i] == 1:
			q.append(i)
	if len(p)>len(q):
		mincost = math.inf
		for i in range(0,len(p)):
			mincost = min(mincost,cost(p[i],none,g1,g2))
		return (len(p)-len(q))*mincost
	else:
		mincost = math.inf
		for i in range(0,len(q)):
			mincost = min(mincost,cost(none,q[i],g1,g2))
		return (len(q)-len(p))*mincost
def lower_graph(g1,g2,cost):
	v1 = len(g1("adj_contact_map"))
	v2 = len(g2("adj_contact_map"))
	
	
	if v1 > v2:
		mincost = math.inf
		for i in range(0,v1):
			mincost = min(mincost,cost(i,none,g1,g2))
		return (v1-v2)*mincost
	else:
		mincost = math.inf
		for i in range(0,v2):
			mincost = min(mincost,cost(none,i,g1,g2))
			return (v2-v1)*mincost
def cost(u,v,g1,g2):
	return 1
"""a and b are the sets of edges"""
def hed_set(a,b,g1,g2,cost):
	ca = []
	cb = []

	for i in range(0,len(a)):
		ca[x] = cost(a,none,g1,g2)
	for j in range(0,len(b)):
		cb[y] = cost(b,none,g1,g2)

	for x in range(0,len(a)):
		for y in range(0,len(b)):
			ca[x] = min(cost(a[x],b[y],g1,g2)/2,ca[x])
			cb[y] = min(cost(a[x],b[y],g1,g2)/2,cb[x])

	return sum(ca) + sum(cb)
def hed_graph(g1, g2, cost):
	sum = 0
	adj1 = g1["adj_contact_map"]
	adj2 = g2["adj_contact_map"]
	d1 = []
	edges1 = []
	edges2 = []
	for u in range(0,len(adj1)):
		for e in range(0,len(adj1)):
			if adj1[u][e] == 1:
				edges1.append(cost(j,none,g1,g2)/2)
		d1[i] = cost(i,none)+sum(edges)

	d2 = []
	for v in range(0,len(adj2)):
		edges = []
		for e in range(0,len(adj2)):
			if adj1[v][e] == 1:
				edges2.append(cost(j,none,g1,g2)/2)
		d2[v] = cost(v,none,g1,g2)+sum(edges)

	for u in range(0,len(adj1)):
		for v in range(0,len(adj2)):
			cost_edges = cost_hed(edges1,edges2,g1,g2,cost)
			cost_edges = max(lower_vertex(u,v,g1,g2),cost_edges)
			
			d1[u] = min((cost(u,v,g1,g2)+cost_edges/2)/2,d1[u])
			d2[v] = min((cost(u,v,g1,g2)+cost_edges/2)/2,d2[v])
	d = sum(d1)+sum(d2)
	d = max(lower_graph(g1,g2),d)

	return d

def cost_simple(i, j, g1, g2):
	return 1
	
def cost_hausdorff(i, j, g1, g2):
	return 1
	
def cost_weighted(i, j, g1, g2):
	adj1 = g1["adj_contact_map"]
	adj2 = g2["adj_contact_map"]
	sum = 0
	g_use = g1
	if g2["adj_contact_map"][i][j] == 1:
		g_use = g2
		
	# inverse of distance between points i and j
	return g_use["adj_inverse_dcalpha"][i][j]

""" We are looking at graphs with the same number of nodes,and
similar node ordering, so we can just do edge comparisons
to get a hopefully good edit distance.
"""
def edge_test(g1, g2, cost):
	sum = 0
	adj1 = g1["adj_contact_map"]
	adj2 = g2["adj_contact_map"]
	for i in range(0, len(adj1)):
		for j in range(0, len(adj1[i])):
			if adj1[i][j] != adj2[i][j]:
				sum += cost(i, j, g1, g2)
	
	return sum
	

def test_similarity(i1, i2):
	g1 = vals[i1]
	g2 = vals[i2]
	print("Indices: " + str(i1) + " " + str(i2))
	
	ct = curr_time()
	eig_res = eigenvalue_method(g1, g2)
	ct2 = curr_time()
	print("Eigenvalue method: " + str(eig_res))
	print("(^^ took " + str(ct2 - ct) + " ms)")
	ct = ct2
	
	
	euc_res = euclidean_method(g1, g2)
	ct2 = curr_time()
	print("Euclidean method: " + str(euc_res))
	print("(^^ took " + str(ct2 - ct) + " ms)")
	ct = ct2

	rmsd_res = rmsd(g1, g2)
	ct2 = curr_time()
	print("RMSD: " + str(rmsd_res))
	print("(^^ took " + str(ct2 - ct) + " ms)")
	ct = ct2

	sim_score_res = sim_score_method(g1, g2)
	ct2 = curr_time()
	print("Iterative sim. score: " + str(sim_score_res))
	print("(^^ took " + str(ct2 - ct) + " ms)")
	ct = ct2

	edge_res = edge_test(g1, g2, cost_simple)
	ct2 = curr_time()
	print("Edge score simple: " + str(edge_res))
	print("(^^ took " + str(ct2 - ct) + " ms)")
	ct = ct2

	edge_res = edge_test(g1, g2, cost_hausdorff)
	ct2 = curr_time()
	print("Edge score Hausdorff: " + str(edge_res))
	print("(^^ took " + str(ct2 - ct) + " ms)")
	ct = ct2

	edge_res = edge_test(g1, g2, cost_weighted)
	ct2 = curr_time()
	print("Edge score weighted: " + str(edge_res))
	print("(^^ took " + str(ct2 - ct) + " ms)")
	ct = ct2

	print()

	
	
	
# returns True if graphs have same residue but are different
def valid_comparison(i1, i2):
	g1 = vals[i1]
	g2 = vals[i2]

	if (len(g1["aa"]) != len(g2["aa"])):
		return False

	if np.array_equal(g1["calpha_coord"], g2["calpha_coord"]):
		return False
	
	return True
	

def compare_all():
	for i in range(0, len(vals_comp)):
		for j in range(0, len(vals_comp)):
			if valid_comparison(i, j):
				print("valid")
				test_similarity(i, j)


eig_scores = []
rmsd_scores = []
sim_scores = []
edge_simple_scores = []
edge_weighted_scores = []

eig_times = []
rmsd_times = []
sim_times = []
edge_simple_times = []
edge_weighted_times = []

def test_similarity_constructed(i):
	g1 = vals[i]
	g2 = vals_comp[i]
	print("Index: " + str(i))
	
	ct = curr_time()
	eig_res = eigenvalue_method(g1, g2)
	ct2 = curr_time()
	eig_scores.append(eig_res)
	eig_times.append(ct2 - ct)
	ct = ct2
	
	euc_res = euclidean_method(g1, g2)
	ct2 = curr_time()
	ct = ct2

	rmsd_res = rmsd(g1, g2)
	ct2 = curr_time()
	rmsd_scores.append(rmsd_res)
	rmsd_times.append(ct2 - ct)
	ct = ct2

	sim_score_res = sim_score_method(g1, g2)
	ct2 = curr_time()
	sim_scores.append(sim_score_res)
	sim_times.append(ct2 - ct)
	ct = ct2

	edge_res = edge_test(g1, g2, cost_simple)
	ct2 = curr_time()
	edge_simple_scores.append(edge_res)
	edge_simple_times.append(ct2 - ct)
	ct = ct2

	edge_res = edge_test(g1, g2, cost_hausdorff)
	ct2 = curr_time()
	ct = ct2

	edge_res = edge_test(g1, g2, cost_weighted)
	ct2 = curr_time()
	edge_weighted_scores.append(edge_res)
	edge_weighted_times.append(ct2 - ct)
	ct = ct2
	print(flush=True)

	
	
	
# returns True if graphs have same residue but are different
def valid_comparison_constructed(i):
	g1 = vals[i]
	g2 = vals_comp[i]

	if (len(g1["aa"]) != len(g2["aa"])):
		return False

	#if np.array_equal(g1["calpha_coord"], g2["calpha_coord"]):
		#return False
	
	return True
	

def compare_all_constructed():
	for i in range(0, len(vals_comp)):
		if valid_comparison_constructed(i):
				print("valid")
				test_similarity_constructed(i)
				
		print(eig_scores[i])
		print(rmsd_scores[i])
		print(sim_scores[i])
		print(edge_simple_scores[i])
		print(edge_weighted_scores[i])
		
		print(eig_times[i])
		print(rmsd_times[i])
		print(sim_times[i])
		print(edge_simple_times[i])
		print(edge_weighted_times[i])
		print()

vals_comp = []
def construct_comp_graphs():
	num_comp = 20
	for g in range(0, num_comp):
		to_add = copy.deepcopy(vals[g])
		for i in range(0, len(to_add["aa"])):
			for j in range(0, len(to_add["aa"])):
				curr_dist = 1/to_add["adj_inverse_dcalpha"][i][j]
				add_num = random.uniform(-1 * curr_dist/3, curr_dist/3)
				new_val = 1/(curr_dist + add_num)
				to_add["adj_inverse_dcalpha"][i][j] = new_val
				to_add["adj_inverse_dcalpha"][j][i] = new_val
				
				if curr_dist + add_num > 10:
					to_add["adj_contact_map"][i][j] = 0
					to_add["adj_contact_map"][j][i] = 0
				else:
					to_add["adj_contact_map"][i][j] = 1
					to_add["adj_contact_map"][j][i] = 1
		vals_comp.append(to_add)

construct_comp_graphs()
compare_all_constructed()