#!/bin/python

import sys
import numpy as np
from pathlib import Path
import cvxopt
from GeneGraph import *
from trie_conversion import *


def linearprogramming_lb(old_g, new_g, res, tname):
	# old_g is the original gene graph, assuming each edge contain the information of its incident transcripts
	# new_g is the prefix graph of the corresponding gene
	# tname is the name of the query transcripts
	trans = sum([e.IncidentTranscripts for e in old_g.vEdges], [])
	trans = list(set(trans))
	# min: c^Tx -- c is zero for not tname and one for tname; x is the varialbe for transcript abundance
	# st: Gx <= h -- G is negative identity matrix; h is zero vector
	#     Ax = b -- A (num_edge * num_trans) is the bibnary indicator of whether the transcript contain the edge; b is the edge abundance
	c = np.zeros(len(trans))
	c[trans.index(tname)] = 1
	G = np.eye(len(trans)) * -1
	h = np.zeros(len(trans))
	A = np.zeros( (len(new_g.edges), len(trans)) )
	for i in range(len(trans)):
		old_edges = [e.ID for e in old_g.vEdges if trans[i] in e.IncidentTranscripts]
		old_nodes = [old_g.vEdges[e_id].Node_1 for e_id in old_edges] + [old_g.vEdges[old_edges[-1]].Node_2]
		assert(old_nodes[0] == old_g.vNodes[0].ID and old_nodes[-1] == old_g.vNodes[-1].ID)
		new_nodes, new_edges = new_g.walk(old_nodes)
		A[new_edges,i] = 1
	b = res
	# convert to cvxopt matrice
	c = cvxopt.matrix(c)
	G = cvxopt.matrix(G)
	h = cvxopt.matrix(h)
	A = cvxopt.matrix(A)
	b = cvxopt.matrix(b)
	sol = cvxopt.solvers.lp(c = c, G = G, h = h, A = A, b = b, solver='glpk', options={'glpk':{'msg_lev':'GLP_MSG_OFF'}})
	return max(sol['primal objective'], 0)


def linearprogramming_ub(old_g, new_g, res, tname):
	# old_g is the original gene graph, assuming each edge contain the information of its incident transcripts
	# new_g is the prefix graph of the corresponding gene
	# tname is the name of the query transcripts
	trans = sum([e.IncidentTranscripts for e in old_g.vEdges], [])
	trans = list(set(trans))
	# min: c^Tx -- c is zero for not tname and NEGATIVE one for tname; x is the varialbe for transcript abundance
	# st: Gx <= h -- G is negative identity matrix; h is zero vector
	#     Ax = b -- A (num_edge * num_trans) is the bibnary indicator of whether the transcript contain the edge; b is the edge abundance
	c = np.zeros(len(trans))
	c[trans.index(tname)] = -1
	G = np.eye(len(trans)) * -1
	h = np.zeros(len(trans))
	A = np.zeros( (len(new_g.edges), len(trans)) )
	for i in range(len(trans)):
		old_edges = [e.ID for e in old_g.vEdges if trans[i] in e.IncidentTranscripts]
		old_nodes = [old_g.vEdges[e_id].Node_1 for e_id in old_edges] + [old_g.vEdges[old_edges[-1]].Node_2]
		assert(old_nodes[0] == old_g.vNodes[0].ID and old_nodes[-1] == old_g.vNodes[-1].ID)
		new_nodes, new_edges = new_g.walk(old_nodes)
		A[new_edges,i] = 1
	b = res
	# convert to cvxopt matrice
	c = cvxopt.matrix(c)
	G = cvxopt.matrix(G)
	h = cvxopt.matrix(h)
	A = cvxopt.matrix(A)
	b = cvxopt.matrix(b)
	sol = cvxopt.solvers.lp(c = c, G = G, h = h, A = A, b = b, solver='glpk', options={'glpk':{'msg_lev':'GLP_MSG_OFF'}})
	return max(-sol['primal objective'], 0)


def GetTransGeneMap(old_graphs):
	trans_gene_map = {}
	for old_g in old_graphs:
		trans = sum([e.IncidentTranscripts for e in old_g.vEdges], [])
		trans = set(trans)
		for tname in trans:
			trans_gene_map[tname] = old_g.GeneID
	return trans_gene_map


def WriteTranscriptBounds(bounds, outputfile):
	fp = open(outputfile, 'w')
	fp.write("# Name\tlb\tub\n")
	for t,v in bounds.items():
		fp.write("{}\t{}\t{}\n".format(t, v[0], v[1]))
	fp.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python lp_fixed_transcripts.py <graph_prefix> <result_object> <outputfile>")
	else:
		graph_prefix = sys.argv[1]
		result_object = sys.argv[2]
		outputfile = sys.argv[3]

		if not Path(outputfile).exists():
			old_graphs = ReadGraphFile(graph_prefix + "_graph_fragstart.txt")
			old_name_index = {old_graphs[i].GeneID:i for i in range(len(old_graphs))}
			new_graphs, eq_classes = load_data(graph_prefix)
			results = pickle.load(open(result_object, 'rb'))
			assert(len(old_graphs) == len(results) and len(old_graphs) == len(new_graphs))

			trans_gene_map = GetTransGeneMap(old_graphs)
			bounds = {}
			for tname,gname in trans_gene_map.items():
				lb = linearprogramming_lb(old_graphs[old_name_index[gname]], new_graphs[gname], results[gname], tname)
				ub = linearprogramming_ub(old_graphs[old_name_index[gname]], new_graphs[gname], results[gname], tname)
				bounds[tname] = (lb,ub)
			WriteTranscriptBounds(bounds, outputfile)
