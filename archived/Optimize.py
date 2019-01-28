#!/bin/python

import sys
import numpy as np
import pysam
import scipy.optimize
from GeneGraph import *
from TranscriptClass import *
from TranscriptAlignmentClass import *
import cvxopt
import tqdm
import copy


def ReadSalmonExpression(filename):
	Exp = {}
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		Exp[strs[0]] = float(strs[4])
	return Exp


def FlowLinearProgramming(g):
	n = len(g.vNodes) - 2 # number nodes except the source/sink node
	m = len(g.vEdges) # number edges
	# variable: coverage difference, estimated flow
	c = np.concatenate( (np.ones(n), np.zeros(m)) )
	# equality constraints for flow
	A = np.zeros((n, n+m)) # further not consider sink node, so n-1
	for i in range(n):
		for e in g.vNodes[i+1].InEdges:
			A[i, n+e] = 1
		for e in g.vNodes[i+1].OutEdges:
			A[i, n+e] = -1
	b = np.zeros((n,1))
	# inequality constraints
	G = np.zeros((2*n+m, n+m))
	# top left middle left : negative identity matrix
	for i in range(n):
		G[i,i] = -1
		G[i+n, i] = -1
	# top right: out edges; middle right: negative out edges
	for i in range(n):
		for e in g.vNodes[i+1].OutEdges:
			G[i, n+e] = g.vNodes[i+1].BiasMultiplier
			G[i+n, n+e] = -g.vNodes[i+1].BiasMultiplier
	# bottom left: 0; bottom right: negative identity (to ensure positive flow)
	for i in range(m):
		G[2*n+i, n+i] = -1
	h = np.zeros((2*n+m, 1))
	for i in range(n):
		h[i,0] = g.vNodes[i+1].Coverage
		h[i+n, 0] = -g.vNodes[i+1].Coverage
	# convert to cvxopt matrix
	c = matrix(c)
	A = matrix(A)
	b = matrix(b)
	G = matrix(G)
	h = matrix(h)
	sol = solvers.lp(c, G, h, A, b)
	return sol['x']


class OptimizationObject_hyper(object):
	def __init__(self, g, eqclass):
		self.graph = g
		self.PathGroups = []
		self.offset_hyper = 0
		self.Weights = []
		self.Weights_fixed = []
		self.Weights_changing = []
		self.Results = None
		# creating temporary path maps
		pathweight_fixed = {}
		pathweight_chaning = {}
		for eq in eqclass:
			assert(len(eq.GeneIDs) == len(eq.Aligned_Paths))
			if len(eq.GeneIDs) == 1:
				gname = eq.GeneIDs[0]
				if gname == self.graph.GeneID:
					if eq.Aligned_Paths[0] in pathweight_fixed:
						pathweight_fixed[eq.Aligned_Paths[0]] += eq.Weights[0]
					else:
						pathweight_fixed[eq.Aligned_Paths[0]] = eq.Weights[0]
			else:
				for i in range(len(eq.GeneIDs)):
					gname = eq.GeneIDs[i]
					if gname == self.graph.GeneID:
						if eq.Aligned_Paths[i] in pathweight_chaning:
							pathweight_chaning[eq.Aligned_Paths[i]] += eq.Weights[i]
						else:
							pathweight_chaning[eq.Aligned_Paths[i]] = eq.Weights[i]
		# merge paths and weights
		self.PathGroups = list( set(pathweight_fixed.keys()) | set(pathweight_chaning.keys()) )
		# first sort by length, then sort by node / edge index
		self.PathGroups.sort(key = lambda x : (len(x), x))
		self.offset_hyper = len(self.PathGroups)
		if len(self.PathGroups) > 0 and len(self.PathGroups[-1]) > 3:
			self.offset_hyper = [i for i in range(len(self.PathGroups)) if len(self.PathGroups[i]) > 3][0]
		for pg in self.PathGroups:
			if pg in pathweight_fixed:
				self.Weights_fixed.append( pathweight_fixed[pg] )
			else:
				self.Weights_fixed.append( 0 )
			if pg in pathweight_chaning:
				self.Weights_changing.append( pathweight_chaning[pg] )
			else:
				self.Weights_changing.append( 0 )
		assert(len(self.Weights_fixed) == len(self.Weights_changing))
		# update total weight
		self.Weights = [self.Weights_fixed[i] + self.Weights_changing[i] for i in range(len(self.Weights_fixed))]

	def update_weights_chaning(self, eqclass):
		# set existing Weights_changing to 0
		for i in range(len(self.Weights_changing)):
			self.Weights_changing[i] = 0
		Index = {self.PathGroups[i]:i for i in range(len(self.PathGroups))}
		for eq in eqclass:
			assert(len(eq.GeneIDs) > 1)
			for i in range(eq.GeneIDs):
				gname = eq.GeneIDs[i]
				if gname == self.graph.GeneID:
					assert(eq.Aligned_Paths[i] in Index)
					idx = Index[eq.Aligned_Paths[i]]
					self.Weights_changing[idx] += eq.Weights[i]
		# update total weight
		assert(len(self.Weights_fixed) == len(self.Weights_changing))
		for i in range(len(self.Weights_fixed)):
			self.Weights.append( self.Weights_fixed[i] + self.Weights_changing[i] )

	def objective(self, flows):
		# flow: first a few elements are normal edge flow, followed by hyper flow corresponding to each sub-path
		# maximization objective: weight_{eq} * np.log( hyper_flow / sum_{all node v} normal_flow_v * length_v )
		# rewritten as: weight_{eq} * np.log( hyper_flow ) - sum(weight_{eq}) * np.log( sum_{all node v} normal_flow_v * length_v )
		assert(len(flows) == len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper)
		value = 0
		for i in range(self.offset_hyper):
			assert(len(self.PathGroups[i]) == 1 or len(self.PathGroups[i]) == 3)
			# read aligned to single node
			if len(self.PathGroups[i]) == 1:
				v = self.PathGroups[i][0]
				node_flow = np.sum([flows[e] for e in self.graph.vNodes[v].OutEdges])
				# print(["node", v, self.Weights[i], node_flow, self.Weights[i] * np.log(node_flow)])
				if node_flow > 0:
					value -= self.Weights[i] * np.log(node_flow)
				else:
					value -= self.Weights[i] * (-50)
			# read occupies single junction / edge
			elif len(self.PathGroups[i]) == 3:
				edge_flow = flows[self.PathGroups[i][1]]
				# print(["edge", self.PathGroups[i][1], self.Weights[i], edge_flow, self.Weights[i] * np.log(edge_flow)])
				if edge_flow > 0:
					value -= self.Weights[i] * np.log(edge_flow)
				else:
					value -= self.Weights[i] * (-50)
		# print(value)
		for i in range(len(self.PathGroups) - self.offset_hyper):
			hyper_flow = flows[len(self.graph.vEdges) + i]
			if hyper_flow > 0:
				value -= self.Weights[self.offset_hyper + i] * np.log(hyper_flow)
			else:
				value -= self.Weights[self.offset_hyper + i] * (-50)
		sum_base = 0
		for e in self.graph.vEdges:
			sum_base += flows[e.ID] * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		value += np.sum(self.Weights) * np.log(sum_base)
		return value

	def derivative(self, flows):
		assert(len(flows) == len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper)
		der = np.zeros(len(flows))
		# normal flow
		for i in range(self.offset_hyper):
			assert(len(self.PathGroups[i]) == 1 or len(self.PathGroups[i]) == 3)
			# read aligned to single node
			if len(self.PathGroups[i]) == 1:
				v = self.PathGroups[i][0]
				node_flow = np.sum([flows[e] for e in self.graph.vNodes[v].OutEdges])
				if node_flow > 0:
					for e in self.graph.vNodes[v].OutEdges:
						der[e] -= self.Weights[i]/ node_flow
				else:
					for e in self.graph.vNodes[v].OutEdges:
						der[e] -= self.Weights[i] / 1e-50
			# read occupies single junction / edge
			elif len(self.PathGroups[i]) == 3:
				edge_flow = flows[self.PathGroups[i][1]]
				if edge_flow > 0:
					der[self.PathGroups[i][1]] -= self.Weights[i] / edge_flow
				else:
					der[self.PathGroups[i][1]] -= self.Weights[i] / 1e-50
		# hyper flow
		for i in range(len(self.PathGroups) - self.offset_hyper):
			hyper_flow = flows[len(self.graph.vEdges) + i]
			if hyper_flow > 0:
				der[len(self.graph.vEdges) + i] -= self.Weights[self.offset_hyper + i] / hyper_flow
			else:
				der[len(self.graph.vEdges) + i] -= self.Weights[self.offset_hyper + i] / 1e-50
		# number base pair normalizer
		sum_weights = np.sum(self.Weights)
		sum_base = 0
		for e in self.graph.vEdges:
			sum_base += flows[e.ID] * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		for e in self.graph.vEdges:
			der[e.ID] += sum_weights / sum_base * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		return der

	def hessian(self, flows):
		assert(len(flows) == len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper)
		hess = np.zeros((len(flows), len(flows)))
		# normal flow
		for i in range(self.offset_hyper):
			assert(len(self.PathGroups[i]) == 1 or len(self.PathGroups[i]) == 3)
			# read aligned to single node
			if len(self.PathGroups[i]) == 1:
				v = self.PathGroups[i][0]
				node_flow = np.sum([flows[e] for e in self.graph.vNodes[v].OutEdges])
				if node_flow > 0:
					for e1 in self.graph.vNodes[v].OutEdges:
						for e2 in self.graph.vNodes[v].OutEdges:
							hess[e1, e2] += self.Weights[i]/ node_flow / node_flow
				else:
					for e1 in self.graph.vNodes[v].OutEdges:
						for e2 in self.graph.vNodes[v].OutEdges:
							hess[e1, e2] += self.Weights[i] / 1e-50 / 1e-50
			# read occupies single junction / edge
			elif len(self.PathGroups[i]) == 3:
				edge_flow = flows[self.PathGroups[i][1]]
				if edge_flow > 0:
					hess[self.PathGroups[i][1], self.PathGroups[i][1]] += self.Weights[i] / edge_flow / edge_flow
				else:
					hess[self.PathGroups[i][1], self.PathGroups[i][1]] += self.Weights[i] / 1e-50 / 1e-50
		# hyper flow
		for i in range(len(self.PathGroups) - self.offset_hyper):
			hyper_flow = flows[len(self.graph.vEdges) + i]
			if hyper_flow > 0:
				hess[len(self.graph.vEdges) + i, len(self.graph.vEdges) + i] += self.Weights[self.offset_hyper + i] / hyper_flow / hyper_flow
			else:
				hess[len(self.graph.vEdges) + i, len(self.graph.vEdges) + i] += self.Weights[self.offset_hyper + i] / 1e-50 / 1e-50
		sum_weights = np.sum(self.Weights)
		sum_base = 0
		for e in self.graph.vEdges:
			sum_base += flows[e.ID] * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		for e1 in self.graph.vEdges:
			for e2 in self.graph.vEdges:
				hess[e1.ID, e2.ID] -= sum_weights / sum_base / sum_base * self.graph.vNodes[e1.Node_1].Length * self.graph.vNodes[e1.Node_1].BiasMultiplier *\
					self.graph.vNodes[e2.Node_1].Length * self.graph.vNodes[e2.Node_1].BiasMultiplier
		return hess

	def constraints(self, numreads = 1):
		# equality
		H = np.zeros((len(self.graph.vNodes)-1, len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper))
		# flow constraint
		for i in range(1, len(self.graph.vNodes)-1):
			v = self.graph.vNodes[i]
			H[i-1, np.array(v.InEdges)] = 1
			H[i-1, np.array(v.OutEdges)] = -1
		# number base pair
		H[-1, :len(self.graph.vEdges)] = np.array([self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier for e in self.graph.vEdges])
		# bound
		bound = np.zeros(len(self.graph.vNodes)-1)
		bound[-1] = numreads
		# inequality: normal flow of e/v >= sum of hyper flow that passes through e/v
		# A * flows >= lb
		A_node = np.zeros((len(self.graph.vNodes) - 1, len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper))
		A_edge = np.zeros((len(self.graph.vEdges), len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper))
		for v in self.graph.vNodes:
			if v.ID != len(self.graph.vNodes) - 1:
				A_node[v.ID, v.OutEdges] = 1
		for e in self.graph.vEdges:
			A_edge[e.ID, e.ID] = 1
		for i in range(len(self.PathGroups) - self.offset_hyper):
			path = self.PathGroups[self.offset_hyper + i]
			# for nodes in the hyper edge
			for j in range(0, len(path), 2):
				v = path[j]
				A_node[v, len(self.graph.vEdges) + i] = -1
			# for edges in the hyper edge:
			for j in range(1, len(path), 2):
				e = path[j]
				A_edge[e, len(self.graph.vEdges) + i] = -1
		A = np.vstack((A_node, A_edge))
		lb = np.zeros(len(self.graph.vNodes) - 1 + len(self.graph.vEdges))
		ub = np.ones(len(self.graph.vNodes) - 1 + len(self.graph.vEdges)) * np.inf
		return H, bound, A, lb, ub

	def initialflow(self, numreads = 1):
		flow_normal = np.zeros(len(self.graph.vEdges))
		flow_hyper = np.zeros(len(self.PathGroups) - self.offset_hyper)
		flow_normal[np.array(self.graph.vNodes[0].OutEdges)] = 1.0 / len(np.array(self.graph.vNodes[0].OutEdges))
		while len(np.where(flow_normal == 0)[0]) > 0:
			for e in np.where(flow_normal == 0)[0]:
				v = self.graph.vEdges[e].Node_1
				if np.all(np.array([flow_normal[e2] for e2 in self.graph.vNodes[v].InEdges]) > 0):
					inflow = np.sum([flow_normal[e2] for e2 in self.graph.vNodes[v].InEdges])
					flow_normal[ np.array(self.graph.vNodes[v].OutEdges) ] = inflow / len(self.graph.vNodes[v].OutEdges)
		# normalize the flow_normal to make sure sum of flow_normal * length = numreads
		s = np.sum([flow_normal[i] * self.graph.vNodes[self.graph.vEdges[i].Node_1].Length * self.graph.vNodes[self.graph.vEdges[i].Node_1].BiasMultiplier for i in range(len(flow_normal))])
		flow_normal *= numreads / s
		# initialize hyper flow
		min_normal_flow = np.min(flow_normal)
		flow_hyper[:] = min_normal_flow / len(self.PathGroups)
		return np.concatenate((flow_normal, flow_hyper))


class CVXOPTObject_split(object):
	def __init__(self, g, efflen_edgewise, eq_classes):
		self.GeneID = g._name
		self.EffLen = efflen_edgewise
		self.Nodes = g.nodes
		self.Edges = g.edges
		self.Aligned_Edge_Groups = []
		self.Weights = []
		self.Weights_fixed = []
		self.Weights_changing = []
		self.Results = None
		# temporary map for edge groups and its weight
		egweight_fixed = {}
		egweight_changing = {}
		for eq in eq_classes:
			if len(eq) == 1:
				assert(eq[0].gid == self.GeneID)
				if tuple(eq[0].new_edges) in egweight_fixed:
					egweight_fixed[tuple(eq[0].new_edges)] += eq[0].weights
				else:
					egweight_fixed[tuple(eq[0].new_edges)] = eq[0].weights
			else:
				for i in range(len(eq)):
					if eq[i].gid == self.GeneID:
						if tuple(eq[i].new_edges) in egweight_changing:
							egweight_changing[tuple(eq[i].new_edges)] += eq[i].weights
						else:
							egweight_changing[tuple(eq[i].new_edges)] = eq[i].weights
		self.Aligned_Edge_Groups = list(set(egweight_fixed.keys()) | set(egweight_changing.keys()))
		self.Aligned_Edge_Groups.sort()
		for k in self.Aligned_Edge_Groups:
			if k in egweight_fixed:
				self.Weights_fixed.append( egweight_fixed[k] )
			else:
				self.Weights_fixed.append( 0 )
			if k in egweight_changing:
				self.Weights_changing.append(egweight_changing[k] )
			else:
				self.Weights_changing.append( 0 )
		assert(len(self.Weights_fixed) == len(self.Weights_changing))
		# update total weight
		self.Weights = np.array([self.Weights_fixed[i] + self.Weights_changing[i] for i in range(len(self.Weights_fixed))])
		# misc attribute
		self.index_single = np.array([i for i in range(len(self.Aligned_Edge_Groups)) if len(self.Aligned_Edge_Groups[i]) == 1])
		self.index_multi = np.array([i for i in range(len(self.Aligned_Edge_Groups)) if len(self.Aligned_Edge_Groups[i]) > 1])
		self.index_single_flow = np.array([self.Aligned_Edge_Groups[i][0] for i in self.index_single])
		assert(len(self.index_single_flow) == len(set(list(self.index_single_flow))))

	def update_weights_changing(self, eq_classes):
		for i in range(len(self.Weights_changing)):
			self.Weights_changing[i] = 0
		Index = {self.Aligned_Edge_Groups[i]:i for i in range(len(self.Aligned_Edge_Groups))}
		for eq in eq_classes:
			assert(len(eq) > 1)
			for i in range(len(eq)):
				if eq[i].gid == self.GeneID:
					assert(eq[i].new_edges in Index)
					idx = Index[eq[i].new_edges]
					self.Weights_changing[idx] += eq[i].weights
		assert(len(self.Weights_fixed) == len(self.Weights_changing))
		# update total weight
		self.Weights = np.array([self.Weights_fixed[i] + self.Weights_changing[i] for i in range(len(self.Weights_fixed))])

	def objective(self, flows):
		# objective: maximize weight_{eq} * np.log( sum(flows_{eq})
		# subject to sum_{all nodes v }(flows_v * length_v) ) = const
		flows = np.array(flows).flatten()
		assert(len(flows) == len(self.Edges))
		value = 0
		for i in range(len(self.Aligned_Edge_Groups)):
			if np.sum(flows[np.array(list(self.Aligned_Edge_Groups[i]))]) > 0:
				value -= self.Weights[i] * np.log( np.sum(flows[np.array(list(self.Aligned_Edge_Groups[i]))]) )
			else:
				# value -= self.Weights[i] * (-50)
				return None
		return value

	def derivative(self, flows):
		flows = np.array(flows).flatten()
		assert(len(flows) == len(self.Edges))
		der = np.zeros(len(flows))
		# for reads that span multiple edges
		for i in self.index_multi:
			if np.sum(flows[np.array(list(self.Aligned_Edge_Groups[i]))]) > 0:
				der[np.array(list(self.Aligned_Edge_Groups[i]))] -= self.Weights[i] / np.sum(flows[np.array(list(self.Aligned_Edge_Groups[i]))])
			else:
				der[np.array(list(self.Aligned_Edge_Groups[i]))] -= self.Weights[i] / 1e-50
		# for reads that span single edge
		der[self.index_single_flow] -= self.Weights[self.index_single] / flows[self.index_single_flow]
		return der

	def hessian(self, flows):
		flows = np.array(flows).flatten()
		assert(len(flows) == len(self.Edges))
		hess = np.zeros((len(flows), len(flows)))
		# for reads that span multiple edges
		for i in self.index_multi:
			if np.sum(flows[np.array(list(self.Aligned_Edge_Groups[i]))]) > 0:
				hess[np.ix_(list(self.Aligned_Edge_Groups[i])), np.ix_(list(self.Aligned_Edge_Groups[i]))] += self.Weights[i] / np.sum(flows[np.array(list(self.Aligned_Edge_Groups[i]))]) / np.sum(flows[np.array(list(self.Aligned_Edge_Groups[i]))])
			else:
				hess[np.ix_(list(self.Aligned_Edge_Groups[i])), np.ix_(list(self.Aligned_Edge_Groups[i]))] += self.Weights[i] / 1e-50 / 1e-50
		# for reads that span single edge
		hess[self.index_single_flow, self.index_single_flow] += self.Weights[self.index_single] / flows[self.index_single_flow] / flows[self.index_single_flow]
		return hess

	def F(self, x = None, z=None):
		if x is None:
			return (0, cvxopt.matrix(self.initialflow()))
		elif not (x is None) and (z is None):
			obj = self.objective(x)
			der = self.derivative(x)
			if not (der is None):
				der = cvxopt.matrix(der).trans()
			return (obj, der)
		else:
			obj = self.objective(x)
			der = self.derivative(x)
			if not (der is None):
				der = cvxopt.matrix(der).trans()
			hess = self.hessian(x)
			if not (hess is None):
				hess = cvxopt.matrix(hess)
			return (obj, der, z[0] * hess)

	def constraints(self, numreads = 100):
		# equality
		H = np.zeros((len(self.Nodes)-1, len(self.Edges)))
		# flow constraint
		for i in range(len(self.Edges)):
			e = self.Edges[i]
			if e[0] != 0:
				H[e[0]-1, i] = -1
			if e[1] != len(self.Nodes)-1:
				H[e[1]-1, i] = 1
		# number base pair
		H[-1, :len(self.Edges)] = np.array(self.EffLen)
		# bound
		bound = np.zeros(len(self.Nodes)-1)
		bound[-1] = numreads
		# non-negative constraint
		A = np.eye(len(self.Edges))
		lb = np.zeros(len(self.Edges))
		return cvxopt.matrix(H), cvxopt.matrix(bound), cvxopt.matrix(A), cvxopt.matrix(lb)

	def initialflow(self, numreads = 100):
		if not (self.Results is None):
			s = np.dot(self.EffLen, self.Results)
			return self.Results * numreads / s
		else:
			node_inedges = [[] for i in range(len(self.Nodes))]
			node_outedges = [[] for i in range(len(self.Nodes))]
			for i in range(len(self.Edges)):
				node_inedges[self.Edges[i][1]].append(i)
				node_outedges[self.Edges[i][0]].append(i)
			flow_edges = np.zeros(len(self.Edges))
			node_visitied = np.zeros(len(self.Nodes), dtype=np.int8)
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
			s = np.dot(flow_edges, self.EffLen)
			flow_edges *= numreads / s
			return flow_edges


class CVXOPTObject_hyper(object):
	def __init__(self, g, eqclass):
		self.graph = g
		self.PathGroups = []
		self.offset_hyper = 0
		self.Weights = []
		self.Weights_fixed = []
		self.Weights_changing = []
		self.Results = None
		# creating temporary path maps
		pathweight_fixed = {}
		pathweight_chaning = {}
		for eq in eqclass:
			assert(len(eq.GeneIDs) == len(eq.Aligned_Paths))
			if len(eq.GeneIDs) == 1:
				gname = eq.GeneIDs[0]
				if gname == self.graph.GeneID:
					if eq.Aligned_Paths[0] in pathweight_fixed:
						pathweight_fixed[eq.Aligned_Paths[0]] += eq.Weights[0]
					else:
						pathweight_fixed[eq.Aligned_Paths[0]] = eq.Weights[0]
			else:
				for i in range(len(eq.GeneIDs)):
					gname = eq.GeneIDs[i]
					if gname == self.graph.GeneID:
						if eq.Aligned_Paths[i] in pathweight_chaning:
							pathweight_chaning[eq.Aligned_Paths[i]] += eq.Weights[i]
						else:
							pathweight_chaning[eq.Aligned_Paths[i]] = eq.Weights[i]
		# merge paths and weights
		self.PathGroups = list( set(pathweight_fixed.keys()) | set(pathweight_chaning.keys()) )
		# first sort by length, then sort by node / edge index
		self.PathGroups.sort(key = lambda x : (len(x), x))
		self.offset_hyper = len(self.PathGroups)
		if len(self.PathGroups) > 0 and len(self.PathGroups[-1]) > 3:
			self.offset_hyper = [i for i in range(len(self.PathGroups)) if len(self.PathGroups[i]) > 3][0]
		for pg in self.PathGroups:
			if pg in pathweight_fixed:
				self.Weights_fixed.append( pathweight_fixed[pg] )
			else:
				self.Weights_fixed.append( 0 )
			if pg in pathweight_chaning:
				self.Weights_changing.append( pathweight_chaning[pg] )
			else:
				self.Weights_changing.append( 0 )
		assert(len(self.Weights_fixed) == len(self.Weights_changing))
		# update total weight
		self.Weights = np.array([self.Weights_fixed[i] + self.Weights_changing[i] for i in range(len(self.Weights_fixed))])

	def update_weights_chaning(self, eqclass):
		# set existing Weights_changing to 0
		for i in range(len(self.Weights_changing)):
			self.Weights_changing[i] = 0
		Index = {self.PathGroups[i]:i for i in range(len(self.PathGroups))}
		for eq in eqclass:
			assert(len(eq.GeneIDs) > 1)
			for i in range(len(eq.GeneIDs)):
				gname = eq.GeneIDs[i]
				if gname == self.graph.GeneID:
					assert(eq.Aligned_Paths[i] in Index)
					idx = Index[eq.Aligned_Paths[i]]
					self.Weights_changing[idx] += eq.Weights[i]
		# update total weight
		assert(len(self.Weights_fixed) == len(self.Weights_changing))
		self.Weights = np.array([self.Weights_fixed[i] + self.Weights_changing[i] for i in range(len(self.Weights_fixed))])

	def objective(self, flows):
		# flow: first a few elements are normal edge flow, followed by hyper flow corresponding to each sub-path
		# maximization objective: weight_{eq} * np.log( hyper_flow / sum_{all node v} normal_flow_v * length_v )
		# rewritten as: weight_{eq} * np.log( hyper_flow ) - sum(weight_{eq}) * np.log( sum_{all node v} normal_flow_v * length_v )
		flows = np.array(flows).flatten()
		assert(len(flows) == len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper)
		value = 0
		for i in range(self.offset_hyper):
			assert(len(self.PathGroups[i]) == 1 or len(self.PathGroups[i]) == 3)
			# read aligned to single node
			if len(self.PathGroups[i]) == 1:
				v = self.PathGroups[i][0]
				node_flow = np.sum([flows[e] for e in self.graph.vNodes[v].OutEdges])
				# print(["node", v, self.Weights[i], node_flow, self.Weights[i] * np.log(node_flow)])
				if node_flow > 0:
					value -= self.Weights[i] * np.log(node_flow)
				else:
					return None
			# read occupies single junction / edge
			elif len(self.PathGroups[i]) == 3:
				edge_flow = flows[self.PathGroups[i][1]]
				# print(["edge", self.PathGroups[i][1], self.Weights[i], edge_flow, self.Weights[i] * np.log(edge_flow)])
				if edge_flow > 0:
					value -= self.Weights[i] * np.log(edge_flow)
				else:
					return None
		value -= np.sum(self.Weights[self.offset_hyper:] * np.log(np.maximum(flows[len(self.graph.vEdges):], 1e-50)))
		return value

	def derivative(self, flows):
		flows = np.array(flows).flatten()
		assert(len(flows) == len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper)
		der = np.zeros(len(flows))
		# normal flow
		for i in range(self.offset_hyper):
			assert(len(self.PathGroups[i]) == 1 or len(self.PathGroups[i]) == 3)
			# read aligned to single node
			if len(self.PathGroups[i]) == 1:
				v = self.PathGroups[i][0]
				node_flow = np.sum([flows[e] for e in self.graph.vNodes[v].OutEdges])
				der[np.array(self.graph.vNodes[v].OutEdges)] -= self.Weights[i]/ np.maximum(node_flow, 1e-50)
		# read occupies single junction / edge
		single_edge_paths = np.array([i for i in range(self.offset_hyper) if len(self.PathGroups[i]) == 3])
		single_edges = np.array([self.PathGroups[i][1] for i in single_edge_paths])
		assert(len(single_edge_paths) == len(single_edges))
		if len(single_edges) > 0:
			der[single_edges] -= self.Weights[single_edge_paths] / np.maximum(flows[single_edges], 1e-50)
		# hyper flow
		der[len(self.graph.vEdges):] -= self.Weights[self.offset_hyper:] / np.maximum(flows[len(self.graph.vEdges):], 1e-50)
		return der

	def hessian(self, flows):
		flows = np.array(flows).flatten()
		assert(len(flows) == len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper)
		hess = np.zeros((len(flows), len(flows)))
		# normal flow
		for i in range(self.offset_hyper):
			assert(len(self.PathGroups[i]) == 1 or len(self.PathGroups[i]) == 3)
			# read aligned to single node
			if len(self.PathGroups[i]) == 1:
				v = self.PathGroups[i][0]
				node_flow = np.sum([flows[e] for e in self.graph.vNodes[v].OutEdges])
				hess[np.ix_(self.graph.vNodes[v].OutEdges, self.graph.vNodes[v].OutEdges)] +=  self.Weights[i]/ np.square(np.maximum(node_flow, 1e-50))
		# read occupies single junction / edge
		hess_diag = np.diagonal(hess).copy()
		single_edge_paths = np.array([i for i in range(self.offset_hyper) if len(self.PathGroups[i]) == 3])
		single_edges = np.array([self.PathGroups[i][1] for i in single_edge_paths])
		assert(len(single_edge_paths) == len(single_edges))
		if len(single_edges) > 0:
			hess_diag[single_edges] += self.Weights[i] / np.square(np.maximum(flows[single_edges], 1e-50))
		# hyper flow
		hess_diag[len(self.graph.vEdges):] += self.Weights[self.offset_hyper:] / np.square(np.maximum(flows[len(self.graph.vEdges):], 1e-50))
		np.fill_diagonal(hess, hess_diag)
		return hess

	def F(self, x = None, z=None):
		if x is None:
			return (0, cvxopt.matrix(self.initialflow()))
		elif not (x is None) and (z is None):
			obj = self.objective(x)
			der = self.derivative(x)
			if not (der is None):
				der = cvxopt.matrix(der).trans()
			return (obj, der)
		else:
			obj = self.objective(x)
			der = self.derivative(x)
			if not (der is None):
				der = cvxopt.matrix(der).trans()
			hess = self.hessian(x)
			if not (hess is None):
				hess = cvxopt.matrix(hess)
			return (obj, der, z[0] * hess)

	def constraints(self, numreads = 1000):
		# equality
		H = np.zeros((len(self.graph.vNodes)-1, len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper))
		# flow constraint
		for i in range(1, len(self.graph.vNodes)-1):
			v = self.graph.vNodes[i]
			H[i-1, np.array(v.InEdges)] = 1
			H[i-1, np.array(v.OutEdges)] = -1
		# number base pair
		H[-1, :len(self.graph.vEdges)] = np.array([self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier for e in self.graph.vEdges])
		# bound
		bound = np.zeros(len(self.graph.vNodes)-1)
		bound[-1] = numreads
		# inequality: normal flow of e/v >= sum of hyper flow that passes through e/v
		# A * flows >= lb
		A_node = np.zeros((len(self.graph.vNodes) - 1, len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper))
		A_edge = np.zeros((len(self.graph.vEdges), len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper))
		A_nonnegative = np.eye( len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper )
		for v in self.graph.vNodes:
			if v.ID != len(self.graph.vNodes) - 1:
				A_node[v.ID, v.OutEdges] = 1
		for e in self.graph.vEdges:
			A_edge[e.ID, e.ID] = 1
		for i in range(len(self.PathGroups) - self.offset_hyper):
			path = self.PathGroups[self.offset_hyper + i]
			# for nodes in the hyper edge
			for j in range(0, len(path), 2):
				v = path[j]
				A_node[v, len(self.graph.vEdges) + i] = -1
			# for edges in the hyper edge:
			for j in range(1, len(path), 2):
				e = path[j]
				A_edge[e, len(self.graph.vEdges) + i] = -1
		A = np.vstack((A_nonnegative, A_node, A_edge))
		lb = np.zeros(len(self.graph.vEdges) + len(self.PathGroups) - self.offset_hyper + len(self.graph.vNodes) - 1 + len(self.graph.vEdges))
		return cvxopt.matrix(H), cvxopt.matrix(bound), cvxopt.matrix(A), cvxopt.matrix(lb)

	def initialflow(self, numreads = 1000):
		if not (self.Results is None):
			return self.Results
		else:
			flow_normal = np.zeros(len(self.graph.vEdges))
			flow_hyper = np.zeros(len(self.PathGroups) - self.offset_hyper)
			flow_normal[np.array(self.graph.vNodes[0].OutEdges)] = 1.0 / len(np.array(self.graph.vNodes[0].OutEdges))
			while len(np.where(flow_normal == 0)[0]) > 0:
				for e in np.where(flow_normal == 0)[0]:
					v = self.graph.vEdges[e].Node_1
					if np.all(np.array([flow_normal[e2] for e2 in self.graph.vNodes[v].InEdges]) > 0):
						inflow = np.sum([flow_normal[e2] for e2 in self.graph.vNodes[v].InEdges])
						flow_normal[ np.array(self.graph.vNodes[v].OutEdges) ] = inflow / len(self.graph.vNodes[v].OutEdges)
			# normalize the flow_normal to make sure sum of flow_normal * length = numreads
			s = np.sum([flow_normal[i] * self.graph.vNodes[self.graph.vEdges[i].Node_1].Length * self.graph.vNodes[self.graph.vEdges[i].Node_1].BiasMultiplier for i in range(len(flow_normal))])
			flow_normal *= numreads / s
			# initialize hyper flow
			edge_involved = []
			for i in range(self.offset_hyper, len(self.PathGroups)):
				edge_involved += [self.PathGroups[i][j] for j in range(1, len(self.PathGroups[i]), 2)]
			if len(edge_involved) > 0:
				edge_involved = np.array(list(set(edge_involved)))
				min_normal_flow = np.min(flow_normal[edge_involved])
				# min_normal_flow = np.min(flow_normal)
				flow_hyper[:] = min_normal_flow / max(len(self.PathGroups) - self.offset_hyper, 1)
			return np.concatenate((flow_normal, flow_hyper))


class OptimizationObject_read(object):
	def __init__(self, g, eqclass):
		self.graph = g
		self.FlowGroups = []
		self.Weights = []
		self.Weights_fixed = []
		self.Weights_changing = []
		# creating temporary maps
		flowweight_fixed = {}
		flowweight_changing = {}
		for eq in eqclass:
			assert(len(eq.GeneIDs) == len(eq.MajorFlowIndex))
			if len(eq.GeneIDs) == 1:
				gname = eq.GeneIDs[0]
				if gname == self.graph.GeneID:
					if eq.MajorFlowIndex[0] in flowweight_fixed:
						flowweight_fixed[eq.MajorFlowIndex[0]] += eq.Weights[0]
					else:
						flowweight_fixed[eq.MajorFlowIndex[0]] = eq.Weights[0]
			else:
				for i in range(len(eq.GeneIDs)):
					gname = eq.GeneIDs[i]
					if gname == self.graph.GeneID:
						if eq.MajorFlowIndex[i] in flowweight_changing:
							flowweight_changing[eq.MajorFlowIndex[i]] += eq.Weights[i]
						else:
							flowweight_changing[eq.MajorFlowIndex[i]] = eq.Weights[i]
		# merge flows and weights
		self.FlowGroups = list(set(flowweight_fixed.keys()) | set(flowweight_changing.keys()))
		self.FlowGroups.sort()
		for fg in self.FlowGroups:
			if fg in flowweight_fixed:
				self.Weights_fixed.append( flowweight_fixed[fg] )
			else:
				self.Weights_fixed.append( 0 )
		for fg in self.FlowGroups:
			if fg in flowweight_changing:
				self.Weights_changing.append( flowweight_changing[fg] )
			else:
				self.Weights_changing.append( 0 )
		# update total weight
		assert(len(self.Weights_fixed) == len(self.Weights_changing))
		for i in range(len(self.Weights_fixed)):
			self.Weights.append( self.Weights_fixed[i] + self.Weights_changing[i] )

	def update_weights_chaning(self, eqclass):
		# set existing Weights_changing to 0
		for i in range(len(self.Weights_changing)):
			self.Weights_changing[i] = 0
		Index = {self.FlowGroups[i]:i for i in range(len(self.FlowGroups))}
		for eq in eqclass:
			assert(len(eq.GeneIDs) > 1)
			for i in range(eq.GeneIDs):
				gname = eq.GeneIDs[i]
				if gname == self.graph.GeneID:
					assert(eq.MajorFlowIndex[i] in Index)
					idx = Index[eq.MajorFlowIndex[i]]
					self.Weights_changing[idx] += eq.Weights[i]
		# update total weight
		assert(len(self.Weights_fixed) == len(self.Weights_changing))
		for i in range(len(self.Weights_fixed)):
			self.Weights.append( self.Weights_fixed[i] + self.Weights_changing[i] )

	def objective(self, flows):
		# objective: maximize weight_{eq} * np.log( sum(flows_{eq}) / sum_{all nodes v }(flows_v * length_v) )
		# rewritten as: maximize weight_{eq} * np.log( sum(flows_{eq}) ) - sum(weight_{eq}) * np.log( sum_{all nodes v }(flows_v * length_v) )
		# rewritten as: maximize weight_{eq} * np.log( sum(flows_{eq}) ) - sum(weight_{eq}) * np.log( sum_{all edges e }(flows_e * length_{node 1 of e}) )
		assert(len(flows) == len(self.graph.vEdges))
		value = 0
		for i in range(len(self.FlowGroups)):
			if np.sum([flows[e] for e in self.FlowGroups[i]]) > 0:
				print([i, self.FlowGroups[i], self.Weights[i], np.sum([flows[e] for e in self.FlowGroups[i]]), self.Weights[i] * np.log(np.sum([flows[e] for e in self.FlowGroups[i]]))])
				value -= self.Weights[i] * np.log(np.sum([flows[e] for e in self.FlowGroups[i]]))
			else:
				# print(flows)
				# print(np.sum([flows[e] for e in self.FlowGroups[i]]))
				value -= self.Weights[i] * (-50)
		print(value)
		sum_base = 0
		for e in self.graph.vEdges:
			sum_base += flows[e.ID] * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		print(np.log(sum_base))
		value += np.sum(self.Weights) * np.log(sum_base)
		return value

	def derivative(self, flows):
		assert(len(flows) == len(self.graph.vEdges))
		der = np.zeros(len(flows))
		for i in range(len(self.FlowGroups)):
			for e in self.FlowGroups[i]:
				if np.sum([flows[e] for e in self.FlowGroups[i]]) > 0:
					der[e] -= self.Weights[i] / np.sum([flows[e] for e in self.FlowGroups[i]])
				else:
					# print(flows)
					# print(der)
					der[e] -= self.Weights[i] / 1e-50
		sum_weights = np.sum(self.Weights)
		sum_base = 0
		for e in self.graph.vEdges:
			sum_base += flows[e.ID] * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		for e in self.graph.vEdges:
			der[e.ID] += sum_weights / sum_base * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		return der

	def hessian(self, flows):
		assert(len(flows) == len(self.graph.vEdges))
		hess = np.zeros((len(flows), len(flows)))
		for i in range(len(self.FlowGroups)):
			if np.sum([flows[e] for e in self.FlowGroups[i]]) > 0:
				tmp_hess = self.Weights[i] / np.sum([flows[e] for e in self.FlowGroups[i]]) / np.sum([flows[e] for e in self.FlowGroups[i]])
			else:
				tmp_hess = self.Weights[i] / 1e-50 / 1e-50
			for e1 in self.FlowGroups[i]:
				for e2 in self.FlowGroups[i]:
					hess[e1, e2] += tmp_hess
		sum_weights = np.sum(self.Weights)
		sum_base = 0
		for e in self.graph.vEdges:
			sum_base += flows[e.ID] * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		for e1 in self.graph.vEdges:
			for e2 in self.graph.vEdges:
				hess[e1.ID, e2.ID] -= sum_weights / sum_base / sum_base * self.graph.vNodes[e1.Node_1].Length * self.graph.vNodes[e1.Node_1].BiasMultiplier *\
					self.graph.vNodes[e2.Node_1].Length * self.graph.vNodes[e2.Node_1].BiasMultiplier
		return hess

	def equconstraint(self, numreads = 1):
		H = np.zeros((len(self.graph.vNodes)-1, len(self.graph.vEdges)))
		# flow constraint
		for i in range(1, len(self.graph.vNodes)-1):
			v = self.graph.vNodes[i]
			H[i-1, np.array(v.InEdges)] = 1
			H[i-1, np.array(v.OutEdges)] = -1
		# flows sum up to 1
		H[-1, :] = np.array([self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier for e in self.graph.vEdges])
		# bound
		bound = np.zeros(len(self.graph.vNodes)-1)
		bound[-1] = numreads
		return H, bound

	def initialflow(self, numreads = 1):
		flow = np.zeros(len(self.graph.vEdges))
		flow[np.array(self.graph.vNodes[0].OutEdges)] = 1.0 / len(np.array(self.graph.vNodes[0].OutEdges))
		while len(np.where(flow == 0)[0]) > 0:
			for e in np.where(flow == 0)[0]:
				v = self.graph.vEdges[e].Node_1
				if np.all(np.array([flow[e2] for e2 in self.graph.vNodes[v].InEdges]) > 0):
					inflow = np.sum([flow[e2] for e2 in self.graph.vNodes[v].InEdges])
					flow[ np.array(self.graph.vNodes[v].OutEdges) ] = inflow / len(self.graph.vNodes[v].OutEdges)
		# normalize the flow to make sure sum of flow * length = numreads
		s = np.sum([flow[i] * self.graph.vNodes[self.graph.vEdges[i].Node_1].Length * self.graph.vNodes[self.graph.vEdges[i].Node_1].BiasMultiplier for i in range(len(flow))])
		flow *= numreads / s
		return flow

	def normalized_sumbase(self, flows, sumflow = 1):
		if np.abs(np.sum(flows[np.array(self.graph.vNodes[0].OutEdges)]) - sumflow) >= 1e-3:
			print( np.sum(flows[np.array(self.graph.vNodes[0].OutEdges)]) )
		assert( np.abs(np.sum(flows[np.array(self.graph.vNodes[0].OutEdges)]) - sumflow) < 1e-3 )
		sum_base = 0
		for e in self.graph.vEdges:
			sum_base += flows[e.ID] * self.graph.vNodes[e.Node_1].Length * self.graph.vNodes[e.Node_1].BiasMultiplier
		return sum_base


def optimizeflow_read(g, eqclass):
	# solve the subproblem: estimate the perfect flow inside each graph
	# assuming the sum of flows is 1. That is, the estimated flow is the proportion inside the given graph
	opt = OptimizationObject_read(g, eqclass)
	H, bound = opt.equconstraint()
	cons = [scipy.optimize.LinearConstraint(H, bound, bound)]
	bnds = scipy.optimize.Bounds(np.zeros(len(g.vEdges)), np.ones(len(g.vEdges)) * 1, keep_feasible=False)
	flows = np.array(opt.initialflow())
	res = scipy.optimize.minimize(opt.objective, flows, method='SLSQP', jac=opt.derivative, hess=opt.hessian, bounds=bnds, constraints=cons)
	for i in range(len(res.x)):
		g.vEdges[i].Flow = max(res.x[i], 0)
	# output normalized_sumbase and sum of read support
	# this two variables are used in calculating the proportion between graphs
	normalized_sumbase = opt.normalized_sumbase(res.x)
	sum_reads = np.sum(opt.Weights)
	return g, normalized_sumbase, sum_reads


def optimizeflow_hyper(g, eqclass, normalize_to_readcount = False):
	# solve the subproblem: estimate the perfect flow inside each graph
	opt = OptimizationObject_hyper(g, eqclass)
	H, bound, A, lb, ub = opt.constraints()
	cons = [scipy.optimize.LinearConstraint(H, bound, bound, keep_feasible=True), scipy.optimize.LinearConstraint(A, lb, ub, keep_feasible=True)]
	bnds = scipy.optimize.Bounds(np.zeros(H.shape[1]), np.ones(H.shape[1]) * 1, keep_feasible=True)
	flows = np.array(opt.initialflow())
	res = scipy.optimize.minimize(opt.objective, flows, method='trust-constr', jac=opt.derivative, hess=opt.hessian, bounds=bnds, constraints=cons)
	if not normalize_to_readcount:
		opt.Results = res.x
		for i in range(len(g.vEdges)):
			g.vEdges[i].Flow = max(res.x[i], 0)
	else:
		s = np.sum(opt.Weights)  / np.sum([res.x[e.ID] * g.vNodes[e.Node_1].Length * g.vNodes[e.Node_1].BiasMultiplier for e in g.vEdges])
		opt.Results = np.maximum(res.x, 0) * s
		for i in range(len(g.vEdges)):
			g.vEdges[i].Flow = s * max(res.x[i], 0)
	return g, opt


def CVXOPToptimize_split(opt, normalize_to_readcount = True):
	H, bound, A, lb = opt.constraints()
	cvxopt.solvers.options['show_progress'] = False
	cvxopt.solvers.options['maxiters'] = 50
	# try:
	sol = cvxopt.solvers.cp(opt.F, G=(-A), h=(-lb), A=H, b=bound)
	if not normalize_to_readcount:
		opt.Results = np.array(sol['x'].trans()).flatten()
	else:
		s = np.sum(opt.Weights)  / np.dot(np.maximum(np.array(sol['x']).flatten(),0), opt.EffLen)
		opt.Results = (np.maximum(np.array(sol['x'].trans()), 0) * s).flatten()
	# except:
	# 	opt.Results = None
	return opt


def CVXOPToptimize_hyper(g, eqclass, normalize_to_readcount = False):
	opt = CVXOPTObject_hyper(g, eqclass)
	H, bound, A, lb = opt.constraints()
	cvxopt.solvers.options['show_progress'] = False
	cvxopt.solvers.options['maxiters'] = 50
	try:
		sol = cvxopt.solvers.cp(opt.F, G=(-A), h=(-lb), A=H, b=bound)
		if not normalize_to_readcount:
			opt.Results = np.array(sol['x'].trans()).flatten()
			for i in range(len(g.vEdges)):
				g.vEdges[i].Flow = max(sol['x'][i], 0)
		else:
			s = np.sum(opt.Weights)  / np.sum([max(sol['x'][e.ID],0) * g.vNodes[e.Node_1].Length * g.vNodes[e.Node_1].BiasMultiplier for e in g.vEdges])
			opt.Results = (np.maximum(np.array(sol['x'].trans()), 0) * s).flatten()
			for i in range(len(g.vEdges)):
				g.vEdges[i].Flow = s * max(sol['x'][i], 0)
		return g, opt
	except:
		opt.Results = None
		return g, opt


def CVXOPToptimize_hyper_cont(g, opt, normalize_to_readcount = False):
	H, bound, A, lb = opt.constraints()
	cvxopt.solvers.options['show_progress'] = False
	cvxopt.solvers.options['maxiters'] = 50
	try:
		sol = cvxopt.solvers.cp(opt.F, G=(-A), h=(-lb), A=H, b=bound)
		if not normalize_to_readcount:
			opt.Results = np.array(sol['x'].trans()).flatten()
			for i in range(len(g.vEdges)):
				g.vEdges[i].Flow = max(sol['x'][i], 0)
		else:
			s = np.sum(opt.Weights)  / np.sum([max(sol['x'][e.ID],0) * g.vNodes[e.Node_1].Length * g.vNodes[e.Node_1].BiasMultiplier for e in g.vEdges])
			opt.Results = (np.maximum(np.array(sol['x'].trans()), 0) * s).flatten()
			for i in range(len(g.vEdges)):
				g.vEdges[i].Flow = s * max(sol['x'][i], 0)
		return g, opt
	except:
		opt.Results = None
		return g, opt


def MaxFlowSinglePath(vNodes, vEdges, capacities):
	# find the single path with max flow can be solved by following the max flow out edge
	path = []
	flow = np.inf
	v = 0
	while v < len(vNodes) - 1:
		idx = np.argmax([capacities[vNodes[v].OutEdges]])
		path.append(vNodes[v].OutEdges[idx])
		flow = np.minimum(capacities[path[-1]], flow)
		v = vEdges[path[-1]].Node_2
	assert(vEdges[path[0]].Node_1 == 0)
	assert(vEdges[path[-1]].Node_2 == len(vNodes)-1)
	return path, flow


def TopKPaths(g, k = 1):
	vNodes = copy.deepcopy(g.vNodes)
	vEdges = copy.deepcopy(g.vEdges)
	capacities = np.array([e.Flow for e in vEdges])
	paths = []
	flows = []
	for r in range(k):
		path, flow = MaxFlowSinglePath(vNodes, vEdges, capacities)
		paths.append(path)
		flows.append(flow)
		capacities[np.array(path)] -= flow
		assert(np.all(capacities >= -1e-6))
	assert(len(paths) == k)
	assert(len(flows) == k)
	return paths, flows


def ReadFlowFromSalmon(salmonquant, Graphs, GraphNameIndex, TransGeneMap):
	# initialize salmon flow vectors
	salmon_flow = []
	for g in Graphs:
		salmon_flow.append( np.zeros(len(g.vEdges)) )
	# read quantification
	fp = open(salmonquant, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		tmp_flow = float(strs[3])
		assert(strs[0] in TransGeneMap)
		gname = TransGeneMap[strs[0]]
		assert(gname in GraphNameIndex)
		idx = GraphNameIndex[gname]
		for e in Graphs[idx].vEdges:
			if strs[0] in e.IncidentTranscripts:
				salmon_flow[idx][e.ID] += tmp_flow
	fp.close()
	# normalize total flow to 1
	for i in range(len(Graphs)):
		sum_flow = np.sum(salmon_flow[i][Graphs[i].vNodes[0].OutEdges])
		if sum_flow > 0:
			salmon_flow[i] /= sum_flow
	return salmon_flow


def ReadHyperFlowFromSalmon(salmonquant, Graphs, GraphNameIndex, TransGeneMap, opts):
	# initialize salmon flow vectors
	salmon_flow = []
	for opt in opts:
		salmon_flow.append( np.zeros(len(opt.graph.vEdges) + len(opt.PathGroups) - opt.offset_hyper) )
	# read quantification
	TPM = {}
	fp = open(salmonquant, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		TPM[strs[0]] = float(strs[3])
	fp.close()
	for i in range(len(opts)):
		opt = opts[i]
		for e in opt.graph.vEdges:
			salmon_flow[i][e.ID] = np.sum([TPM[t] for t in e.IncidentTranscripts if t in TPM])
		for j in range(len(opt.PathGroups) - opt.offset_hyper):
			path = opt.PathGroups[opt.offset_hyper + j]
			assert(len(path) > 3)
			common_trans = set(opt.graph.vEdges[path[1]].IncidentTranscripts)
			for k in range(3, len(path), 2):
				common_trans = common_trans & set(opt.graph.vEdges[path[k]].IncidentTranscripts)
			assert(len(common_trans) > 0)
			salmon_flow[i][len(opt.graph.vEdges) + j] = np.sum([TPM[t] for t in common_trans if t in TPM])
	return salmon_flow


def AdjustAssignment(Opts, NameIndex, eq_classes):
	for eq in eq_classes:
		flows = []
		if len(eq) > 1:
			for ali in eq:
				opt = Opts[NameIndex[ali.gid]]
				if opt.Results is None:
					flows.append( 0 )
				else:
					flows.append( np.sum(opt.Results[eq.new_edges]) )
			if np.sum(flows) != 0:
				flows /= np.sum(flows)
				readcount = np.sum([ali.weights for ali in eq])
				for i in range(len(eq)):
					eq[i].weights = flows[i] * readcount
	return eq_classes