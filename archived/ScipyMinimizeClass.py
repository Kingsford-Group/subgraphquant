#!/bin/python

import sys
import numpy as np
import scipy.optimize
import pyipopt
import tqdm
import copy


class ScipyMinimizeObject_split(object):
	def __init__(self, g, efflen_edgewise, eq_classes):
		# graph info
		self.gid = g._name
		self.efflen = efflen_edgewise
		self.nodes = g.nodes
		self.edges = g.edges
		self.edge_groups = []
		self.edge_group_releindicator = []
		self.weights = []
		self.weights_fixed = []
		self.weights_changing = []
		self.results = None
		# ipopt info
		self.n_var = len(self.edges)
		self.n_cons = len(self.nodes) - 1 # flow, sum of reads
		print("n_var = {}\tn_cons = {}".format(self.n_var, self.n_cons))
		self.nnzj = 2 * len(self.edges) - len([e for e in self.edges if e[0] == 0 or e[1]+1 == len(self.nodes)]) + self.n_var # number of non-zeros in the constraint jacobian
		self.nnzh = int(self.n_var * (self.n_var + 1) / 2) # number of non-zeros in the lagrangian hessian
		# print("nnzj = {}\tnnzh = {}".format(self.nnzj, self.nnzh))
		# process weights
		# temporary map for edge groups and its weights
		egweight_fixed = {}
		egweight_changing = {}
		for eq in eq_classes:
			if len(eq) == 1:
				assert(eq[0].gid == self.gid)
				if tuple(eq[0].new_edges) in egweight_fixed:
					egweight_fixed[tuple(eq[0].new_edges)] += eq[0].weights
				else:
					egweight_fixed[tuple(eq[0].new_edges)] = eq[0].weights
			else:
				for i in range(len(eq)):
					if eq[i].gid == self.gid:
						if tuple(eq[i].new_edges) in egweight_changing:
							egweight_changing[tuple(eq[i].new_edges)] += eq[i].weights
						else:
							egweight_changing[tuple(eq[i].new_edges)] = eq[i].weights
		# final edge groups and weights
		self.edge_groups = list(set(egweight_fixed.keys()) | set(egweight_changing.keys()))
		for eg in self.edge_groups:
			tmp = np.zeros(self.n_var)
			tmp[np.array(eg)] = 1
			self.edge_group_releindicator.append(tmp)
		self.edge_groups.sort()
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
		# update total weights_fixed
		self.weights = np.array([self.weights_fixed[i] + self.weights_changing[i] for i in range(len(self.weights_fixed))])
		# calculate coefficient matrix of constraints
		self.H = np.zeros( (len(self.nodes)-1, self.n_var) )
		for i in range(len(self.edges)):
			e = self.edges[i]
			if e[0] != 0:
				self.H[e[0] - 1, i] = -1
			if e[1] != len(self.nodes) - 1:
				self.H[e[1] - 1, i] = 1
		self.H[-1, :len(self.edges)] = np.array(self.efflen)
		assert(len(np.where(self.H[:-1,] != 0)[0]) == self.nnzj - self.n_var)

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
				value -= self.weights[i] * 1e-50;
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
			rows = np.concatenate( (rows, np.ones(self.H.shape[1]) * (self.H.shape[0]-1)) )
			cols = np.concatenate( (cols, np.array(list(range(self.H.shape[1])))) )
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


def optimizegraph(opt):
	x_L, x_U = opt.bound_variables()
	g_L, g_U = opt.bound_constraints()
	# construct bounds and constraints
	cons = [scipy.optimize.LinearConstraint(opt.H, g_L, g_U, keep_feasible=True)]
	bnds = scipy.optimize.Bounds(x_L, x_U, keep_feasible=True)
	# initial value
	initialf = np.array(opt.initialflow())
	# optimize
	try:
		res = scipy.optimize.minimize(opt.eval_objective, initialf, method='trust-constr', jac=opt.eval_grad_objective, bounds=bnds, constraints=cons)
		# scale to the number of reads
		opt.results = np.maximum(res.x, 0)
		s = np.sum(opt.weights) / np.dot(opt.efflen, opt.results)
		opt.results *= s
		return opt
	except:
		print("error in optimizing gene " + opt.gid)
		return None
