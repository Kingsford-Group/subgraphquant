#!/bin/python

import sys
import numpy as np
from pathlib import Path
from TranscriptClass import *
from TranscriptAlignmentClass import *
from IpoptOptimizeClass import *
from utils import *
from trie_conversion import *
from EffLen_VLMM import *
from flow_graph import FlowGraph
import tqdm
import copy
import pickle
import concurrent.futures
import time
from scipy.sparse import csr_matrix
import scipy.sparse.linalg


def loss_main(opt, t, flows):
	# negative log likelihood
	# f(x) = - sum_eq n_eq * log(sum_e flow_e)
	# n_eq is the number of read support of each equivalent class; flow_e is the prefix graph edge flows that constribute to the hyperedge (read alignment)
	value = 0
	for i in range(len(opt.edge_groups)):
		sum_relevant_flows = np.dot(opt.edge_group_releindicator[i], flows)
		value -= opt.weights[i] * np.log(sum_relevant_flows)
	return t * value


def loss_barrier(opt, flows):
	# - 1 / t * sum_i log(-h_i(x))
	# where t is a large positive number, h_i(x) <= 0 is the original inequality constraint, and h_i(x) = -x
	return -1.0 * np.sum(np.log(flows))


def gradient_main(opt, t, flows):
	# df(x) = - sum_eq n_eq / (sum_e flow_e) * \mathbbm{1}
	grad = np.zeros( len(flows) )
	for i in range(len(opt.edge_groups)):
		sum_relevant_flows = np.dot(opt.edge_group_releindicator[i], flows)
		grad -= opt.weights[i] / sum_relevant_flows * opt.edge_group_releindicator[i]
	return t * grad


def gradient_barrier(opt, flows):
	# - 1 / t * 1 / x, where x is a vector
	return -1.0 / flows


def hessian_main(opt, t, flows):
	# d^2f(x) = 
	hess = np.zeros( (len(flows), len(flows)) )
	for i in range(len(opt.edge_groups)):
		sum_relevant_flows = np.dot(opt.edge_group_releindicator[i], flows)
		indicator = opt.edge_group_releindicator[i].reshape((len(flows), 1))
		hess += opt.weights[i] / sum_relevant_flows / sum_relevant_flows * np.matmul(indicator, indicator.transpose())
	return t * hess


def hessian_barrier(opt, flows):
	# 1 / t * diagonal(1 / x_1^2, 1 / x_2^2, ..., 1 / x_n^2)
	return 1.0 * np.diag(1/flows**2)


def newton_eqconstraint_direction(grad, hess, A):
	# Newton update direction v for x^+ = x + m * v
	# v is solved by argmin_{Az=0} grad^T * z + 0.5 * z^T * hess * z
	# From KKT condition, the optimal v satisfies the linear system:
	# [hess A^T] * [v] = [-grad]
	# [A     0 ]   [w]   [  0  ]
	lhs_row1 = np.hstack( (hess, A.transpose()) )
	lhs_row2 = np.hstack( (A, np.zeros( (A.shape[0], hess.shape[1] + A.shape[0] - A.shape[1]) )) )
	lhs = np.vstack( (lhs_row1, lhs_row2) )
	lhs = csr_matrix(lhs)
	rhs = np.concatenate( (-grad, np.zeros(hess.shape[0] + A.shape[0] - grad.shape[0])) )
	# solve the linear system
	res = scipy.sparse.linalg.spsolve(lhs, rhs)
	return res[:grad.shape[0]]


def Newton_bt(opt, x0, t, m0 = 1, shrinkage = 0.8, alpha = 0.5, stop_criteria = 1e-8, max_iter = 100):
	# newton update with equality constraint x^+ = x + m * v
	# v is solved by newton_eqconstraint_direction
	# m is set to m0, but get shrinked if f(x + mv) > f(x) + alpha * m * grad(x)^T * v
	# initialize flow
	last_x = copy.copy(x0)
	# keep track of loss
	loss_values = [loss_main(opt, t, last_x) + loss_barrier(opt, last_x)]
	while (len(loss_values) < 2 or np.abs(loss_values[-1] - loss_values[-2]) >= stop_criteria) and len(loss_values) < max_iter:
		grad = gradient_main(opt, t, last_x) + gradient_barrier(opt, last_x)
		hess = hessian_main(opt, t, last_x) + hessian_barrier(opt, last_x)
		# newton direction
		v = newton_eqconstraint_direction(grad, hess, opt.H)
		# backtracking line search
		m = m0
		x = last_x + m * v
		while np.any(x < 0) or loss_main(opt, t, x) + loss_barrier(opt, x) > loss_values[-1] - alpha * m * grad.dot(v):
			m *= shrinkage
			x = last_x + m * v
		last_x = x
		loss_values.append( loss_main(opt, t, last_x) + loss_barrier(opt, last_x) )
	return last_x, loss_values


def barrier_method(opt, x0, t0 = 1, m0 = 1, mu = 5, shrinkage = 0.8, alpha = 0.5, stop_criteria = 1e-8, max_iter = 100):
	# solve the objective of min f(x) - 1 / t sum_i log(-h_i(x)) with initial t = t0
	# for each fixed t, solve the equality-constrained problem using Newton method with back tracking line search
	# after each Newton solving, enlarge t by the multiplicative factor mu
	n_inequalities = len(x0)
	loss_blocks = []
	# initialize parameter
	t = t0
	last_x = x0
	while n_inequalities / t >= stop_criteria:
		# with the current t, apply Newton's method
		x, loss_values = Newton_bt(opt, last_x, t, m0, shrinkage, alpha)
		# update x and t and losses
		last_x = x
		t *= mu
		loss_blocks.append(loss_values)
		if len(loss_blocks) > max_iter:
			break
	return last_x, loss_blocks


if __name__=="__main__":
	salmon_dir = "/home/congm1/savanna/savannacong33/GraphSalmonSimulation/run01/salmon_01"
	prefix_graph = "/home/congm1/savanna/savannacong33/GraphSalmonSimulation/run01/salmon_01/graphsalmon/gs"

	graphs, eq_classes = load_data(prefix_graph)
	efflen_hyper = ReadEfflen_hyper(prefix_graph)
	efflen_simple = ReadEfflen_simple(prefix_graph)

	Gene_Eq_Index = {gname:[] for gname in graphs.keys()}
	for i in range(len(eq_classes)):
		for g in [ali.gid for ali in eq_classes[i]]:
			Gene_Eq_Index[g].append(i)
	Gene_Eq_Index = {g:list(set(v)) for g,v in Gene_Eq_Index.items() if len(v) > 0}

	Names = list(Gene_Eq_Index.keys())
	Names.sort()
	NameIndex = {Names[i]:i for i in range(len(Names))}

	# index = 0
	gname = "ENSG00000000003.14"
	opt = IpoptObject_split(graphs[gname], efflen_simple[gname], efflen_hyper[gname], [eq_classes[j] for j in Gene_Eq_Index[gname]])
	x0 = opt.initialflow( np.sum(opt.weights) )

	# index = 10
	gname = "ENSG00000001460.17"
	opt = IpoptObject_split(graphs[gname], efflen_simple[gname], efflen_hyper[gname], [eq_classes[j] for j in Gene_Eq_Index[gname]])
	x0 = opt.initialflow( np.sum(opt.weights) )
	x, loss_blocks = barrier_method(opt, x0)
