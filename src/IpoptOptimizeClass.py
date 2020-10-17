#!/bin/python

import sys
import numpy as np
import scipy.optimize
import pyipopt
import tqdm
import copy
import networkx as nx


class IpoptObject_split(object):
	def __init__(self, g, efflen_simple_paths, efflen_hyper_paths, eq_classes, forbidden_edges = []):
		# note that efflen_simple_paths contain a list of splice graph node list, prefix graph edges and the efflen
		# but efflen_hyper_paths contains the equivalent class mapping node lists and the efflen
		# forbidden_edges is a list of edge indices in the prefix graph constructed with pseudo eq classes corresponding to super short paths
		# graph info
		self.gid = g._name
		self.nodes = g.nodes
		self.edges = g.edges
		self.efflen = np.zeros(len(self.edges))
		self.edge_groups = []
		self.edge_group_releindicator = []
		self.weights = []
		self.weights_fixed = []
		self.weights_changing = []
		self.results = None
		forbidden_edges = set(forbidden_edges)
		# process effective length due to normalization constraint: sum c_p l_p = const
		# map to prefix graph edges: sum_p c_p l_p = sum_p (sum_{e in AS(p)} f_e) l_p = sum_e f_e (sum_{p: e in AS(p)} l_p)
		# prefix graph edge effective length is sum_{p: e in AS(p)} l_p
		# add efflen for simple edges
		processed_nl = []
		for info in efflen_simple_paths:
			(nodelist, edgelist, value) = info
			self.efflen[np.array([x for x in edgelist if not (x in forbidden_edges)])] += value
			processed_nl.append( tuple(nodelist) )
		processed_nl = set(processed_nl)
		# add efflen for hyper edges, get the prefix graph edge corresponding to the splice graph node list from eq_classes
		index_eq = {}
		for i in range(len(eq_classes)):
			eq = eq_classes[i]
			for j in range(len(eq)):
				if eq[j].gid == self.gid:
					index_eq[ tuple(eq[j].vpath) ] = (i,j)
		for info in efflen_hyper_paths:
			(nodelist, value) = info
			# find the corresponding edges by eq_classes
			if tuple(nodelist) in processed_nl:
				continue
			assert( tuple(nodelist) in index_eq )
			(i, j) = index_eq[ tuple(nodelist) ]
			assert( len(forbidden_edges & set(eq_classes[i][j].new_edges)) == 0 )
			self.efflen[ np.array(eq_classes[i][j].new_edges) ] += value
			processed_nl.add( tuple(nodelist) )
		# process weights
		# temporary map for edge groups and its weights
		egweight_fixed = {}
		egweight_changing = {}
		for eq in eq_classes:
			if len(eq) == 1:
				assert(eq[0].gid == self.gid)
				if eq[0].weights < 1e-6:
					continue
				if tuple(eq[0].new_edges) in egweight_fixed:
					egweight_fixed[tuple(eq[0].new_edges)] += eq[0].weights
				else:
					egweight_fixed[tuple(eq[0].new_edges)] = eq[0].weights
			else:
				for i in range(len(eq)):
					if eq[i].gid == self.gid:
						if eq[i].weights > 1e-6:
							if tuple(eq[i].new_edges) in egweight_changing:
								egweight_changing[tuple(eq[i].new_edges)] += eq[i].weights
							else:
								egweight_changing[tuple(eq[i].new_edges)] = eq[i].weights
					# if tuple(eq[i].new_edges) in egweight_changing:
					# 		egweight_changing[tuple(eq[i].new_edges)] += eq[i].weights
					# 	else:
					# 		egweight_changing[tuple(eq[i].new_edges)] = eq[i].weights
		# final edge groups and weights
		self.edge_groups = list(set(egweight_fixed.keys()) | set(egweight_changing.keys()))
		self.edge_groups.sort()
		for eg in self.edge_groups:
			tmp = np.zeros(len(self.edges))
			tmp[np.array(eg)] = 1
			self.edge_group_releindicator.append(tmp)
		for k in self.edge_groups:
			if k in egweight_fixed:
				self.weights_fixed.append( egweight_fixed[k] )
			else:
				self.weights_fixed.append( 0 )
			if k in egweight_changing:
				self.weights_changing.append(egweight_changing[k] )
			else:
				self.weights_changing.append( 0 )
		assert(len(self.weights_fixed) == len(self.weights_changing))
		# remove forbidden edges from the graph
		remaining_edges = np.array( list(set(list(range(len(self.edges)))) - forbidden_edges) )
		self.edges = [self.edges[i] for i in range(len(self.edges)) if not (i in forbidden_edges)]
		self.efflen = self.efflen[remaining_edges]
		self.edge_group_releindicator = [x[remaining_edges] for x in self.edge_group_releindicator]
		assert( len(self.efflen) == len(self.edges) )
		assert(len(self.edge_group_releindicator[0]) == len(self.edges))
		# update total weights_fixed
		self.weights = np.array([self.weights_fixed[i] + self.weights_changing[i] for i in range(len(self.weights_fixed))])
		# calculate coefficient matrix of constraints
		self.H = np.zeros( (len(self.nodes)-1, len(self.edges)) )
		# flow balance
		for i in range(len(self.edges)):
			e = self.edges[i]
			if e[0] != 0:
				self.H[e[0] - 1, i] = -1
			if e[1] != len(self.nodes) - 1:
				self.H[e[1] - 1, i] = 1
		# normalization to a constant
		self.H[-1, :len(self.edges)] = np.array(self.efflen)
		# ipopt info
		self.n_var = len(self.edges)
		self.n_cons = len(self.nodes) - 1 # flow, sum of reads, forbidden edges
		print("n_var = {}\tn_cons = {}".format(self.n_var, self.n_cons))
		self.nnzj = 2 * len(self.edges) - len([e for e in self.edges if e[0] == 0 or e[1]+1 == len(self.nodes)]) + self.n_var # number of non-zeros in the constraint jacobian
		self.nnzh = int(self.n_var * (self.n_var + 1) / 2) # number of non-zeros in the lagrangian hessian
		print("nnzj = {}\tnnzh = {}".format(self.nnzj, self.nnzh))

	def bound_variables(self):
		return np.zeros(self.n_var), np.inf * np.ones(self.n_var)

	def bound_constraints(self, numreads = 100):
		bnd = np.zeros(self.n_cons)
		bnd[-1] = numreads
		return bnd, bnd

	def eval_objective(self, flows):
		# sum of - (weight * np.log(relevance * flows))
		value = 0
		for i in range(len(self.edge_groups)):
			sum_relevant_flows = np.dot(self.edge_group_releindicator[i], flows)
			if sum_relevant_flows > 0:
				value -= self.weights[i] * np.log(sum_relevant_flows)
			else:
				value -= self.weights[i] * np.log(1e-50);
		return value

	def eval_grad_objective(self, flows):
		# sum of (weights / dot_product(relevance, flows) * relevance)
		grad = np.zeros(self.n_var)
		for i in range(len(self.edge_groups)):
			sum_relevant_flows = np.dot(self.edge_group_releindicator[i], flows)
			if sum_relevant_flows > 0:
				grad -= self.weights[i] / sum_relevant_flows * self.edge_group_releindicator[i]
			else:
				grad -= self.weights[i] / 1e-50 * self.edge_group_releindicator[i]
		return grad

	def eval_constraints(self, flows):
		return np.matmul(self.H, flows)

	def eval_jac_constraint(self, flows, flag):
		if flag:
			rows, cols = np.where(self.H[:-1,] != 0)
			rows = np.concatenate( (rows, np.ones(self.H.shape[1], dtype=np.int) * (self.H.shape[0]-1)) )
			cols = np.concatenate( (cols, np.array(list(range(self.H.shape[1])), dtype=np.int)) )
			return (rows, cols)
		else:
			H_nonzeros = self.H[np.where(self.H[:-1,] != 0)]
			H_nonzeros = np.concatenate( (H_nonzeros,self.H[-1,:]) )
			return H_nonzeros

	def eval_hess(self, flows, lagrange, obj_factor, flag):
		if flag:
			return np.tril_indices(self.n_var)
		else:
			hess = np.zeros( (self.n_var, self.n_var) )
			for i in range(len(self.edge_groups)):
				sum_relevant_flows = np.dot(self.edge_group_releindicator[i], flows)
				indicator = self.edge_group_releindicator[i].reshape((self.n_var, 1))
				if sum_relevant_flows > 0:
					hess += self.weights[i] / sum_relevant_flows / sum_relevant_flows * np.matmul(indicator, indicator.transpose())
				else:
					hess += self.weights[i] / 1e-50 / 1e-50 * np.matmul(indicator, indicator.transpose())
			hess *= obj_factor
			index = np.tril_indices(self.n_var)
			return hess[index]

	def initialflow(self, numreads = 100):
		if not(self.results is None):
			s = np.dot(self.efflen, self.results)
			return self.results * numreads / s
		else:
			node_inedges = [[] for i in range(len(self.nodes))]
			node_outedges = [[] for i in range(len(self.nodes))]
			for i in range(len(self.edges)):
				node_inedges[self.edges[i][1]].append(i)
				node_outedges[self.edges[i][0]].append(i)
			flow_edges = np.zeros(len(self.edges))
			node_visitied = np.zeros(len(self.nodes), dtype=np.int8)
			node_visitied[0] = 1
			flow_edges[np.array(node_outedges[0])] = 1.0 / len(node_outedges[0])
			while np.sum(flow_edges == 0) > 0:
				for idx_v in np.where(node_visitied == 0)[0]:
					if np.all(flow_edges[np.array(node_inedges[idx_v])]> 0):
						sum_in = np.sum(flow_edges[np.array(node_inedges[idx_v])])
						flow_edges[np.array(node_outedges[idx_v])] = sum_in / len(node_outedges[idx_v])
						node_visitied[idx_v] = 1
						break
			# scale initial flow
			s = np.dot(flow_edges, self.efflen)
			flow_edges *= numreads / s
			return flow_edges


def reinforce_flowbalance(opt, x):
	node_inedges = [[] for i in range(len(opt.nodes))]
	node_outedges = [[] for i in range(len(opt.nodes))]
	for i in range(len(opt.edges)):
		node_inedges[opt.edges[i][1]].append(i)
		node_outedges[opt.edges[i][0]].append(i)
	visited = np.array( [True] + [False] * (len(opt.nodes)-2) + [True] )
	while np.sum(visited) < len(opt.nodes):
		for idx_v in np.where(visited == False)[0]:
			in_nodes = [opt.edges[idx_e][0] for idx_e in node_inedges[idx_v]]
			if np.sum(visited[in_nodes]) == len(in_nodes):
				# make sure that the outflow sums up to the in flow
				sum_in = np.sum(x[node_inedges[idx_v]])
				sum_out = np.sum(x[node_outedges[idx_v]])
				if sum_out != 0:
					x[node_outedges[idx_v]] *= (sum_in / sum_out)
				elif sum_out == 0 and sum_in != 0:
					x[node_outedges[idx_v]] = 1.0 * sum_in / len(node_outedges[idx_v])
				visited[idx_v] = True
	return x


def optimizegraph(opt, max_iter = 300, max_cpu_time = 100):
	x_L, x_U = opt.bound_variables()
	g_L, g_U = opt.bound_constraints()
	nlp = pyipopt.create(opt.n_var, x_L, x_U, opt.n_cons, g_L, g_U, opt.nnzj, opt.nnzh, \
		opt.eval_objective, opt.eval_grad_objective, opt.eval_constraints, opt.eval_jac_constraint)
	nlp.num_option('max_cpu_time', max_cpu_time)
	nlp.int_option('print_level', 0)
	nlp.num_option('tol', 1e-12)
	nlp.num_option('acceptable_tol', 1e-12)
	x, zl, zu, constraint_multipliers, obj, status = nlp.solve( initialize_flow(opt) )
	nlp.close()
	# reinforce flow balance
	x = reinforce_flowbalance(opt, x)
	# scale x
	x = np.maximum(x, 0)
	s = np.sum(opt.weights) / np.dot(opt.efflen, x)
	opt.results = x * s
	return opt, status


def AdjustAssignment_ipopt(Opts, NameIndex, eq_classes):
	diff = 0
	for eq in eq_classes:
		flows = []
		if len(eq) > 1:
			for ali in eq:
				opt = Opts[NameIndex[ali.gid]]
				assert(not opt.results is None)
				flows.append( np.sum(opt.results[ali.new_edges]) )
			if np.sum(flows) != 0:
				flows /= np.sum(flows)
				readcount = np.sum([ali.weights for ali in eq])
				for i in range(len(eq)):
					diff += np.abs(eq[i].weights - flows[i] * readcount)
					eq[i] = eq[i]._replace(weights = flows[i] * readcount)
	print("Change in read assignment = {}".format(diff))
	return eq_classes


def UpdateOptWeight(Opts, NameIndex, eq_classes, Gene_Eq_Index):
	for opt in Opts:
		index_eq = Gene_Eq_Index[opt.gid]
		# only keep the eq_classes that have more than 1 alignments
		index_eq = [i for i in index_eq if len(eq_classes[i]) > 1]
		# update opt weights_changing
		opt.weights_changing = [0] * len(opt.edge_groups)
		for i in index_eq:
			eq = eq_classes[i]
			for ali in eq:
				if ali.gid == opt.gid:
					assert(tuple(ali.new_edges) in opt.edge_groups)
					idx_edge_group = opt.edge_groups.index( tuple(ali.new_edges) )
					opt.weights_changing[idx_edge_group] += ali.weights
		opt.weights = np.array([opt.weights_fixed[i] + opt.weights_changing[i] for i in range(len(opt.weights_fixed))])
	return Opts


def trace_back(opt, G, idx, flow_edges, map_pairnode_edge):
	# idx is the index of the node, NOT EDGE
	paths_in = [idx]
	while paths_in[0] != 0:
		parents = list(G.predecessors(paths_in[0]))
		this_flows = np.array([flow_edges[map_pairnode_edge[(x, paths_in[0])]] for x in parents])
		if np.sum(this_flows == 0) == 0:
			selection = parents[0]
		else:
			selection = parents[np.where(this_flows==0)[0][0]]
		paths_in = [selection] + paths_in
	return paths_in


def keep_searching(opt, G, idx, flow_edges, map_pairnode_edge):
	# idx is the index of the node, NOT EDGE
	paths_out = [idx]
	while paths_out[-1] != len(opt.nodes) - 1:
		children = list(G.successors(paths_out[-1]))
		this_flows = np.array([flow_edges[map_pairnode_edge[(paths_out[-1],x)]] for x in children])
		if np.sum(this_flows == 0) == 0:
			selection = children[0]
		else:
			selection = children[np.where(this_flows==0)[0][0]]
		paths_out += [selection]
	return paths_out


def initialize_flow(opt, numreads = 100):
	G = nx.DiGraph()
	G.add_nodes_from( list(range(len(opt.nodes))) )
	# edge capacity to be max of (e[0] out-degree, e[0] in-degree)
	G.add_edges_from( opt.edges )
	map_pairnode_edge = {opt.edges[i]:i for i in range(len(opt.edges))}
	flow_edges = np.zeros( len(opt.edges) )
	while np.sum(flow_edges == 0) > 0:
		idx = [i for i in range(len(opt.edges)) if flow_edges[i] == 0][0]
		if opt.edges[idx][0] == 0:
			path = keep_searching(opt, G, 0, flow_edges, map_pairnode_edge)
		elif opt.edges[idx][1] == len(opt.nodes) - 1:
			path = trace_back(opt, G, len(opt.nodes) - 1, flow_edges, map_pairnode_edge)
		else:
			paths_in = trace_back(opt, G, opt.edges[idx][0], flow_edges, map_pairnode_edge)
			paths_out = keep_searching(opt, G, opt.edges[idx][1], flow_edges, map_pairnode_edge)
			path = paths_in + paths_out
		edge_indexes = np.array([map_pairnode_edge[(path[i], path[i+1])] for i in range(len(path) - 1)])
		assert( len(edge_indexes) == len(path) - 1 )
		flow_edges[edge_indexes] += 1
	# normalize such that sum_e f_e l_e = sumreads
	s = flow_edges.dot(opt.efflen)
	flow_edges = flow_edges / s * numreads
	return flow_edges