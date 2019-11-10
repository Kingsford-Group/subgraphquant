#!/bin/python

import sys
import numpy as np
import pickle
from pathlib import Path
from GeneGraph import *
from utils import *
from TranscriptClass import *
from trie_conversion import *


def ProcessEdgeFlow_newgraph(old_graphs, new_graphs, Exp):
	# old_graphs are lists of graphs, of which the edges contain the attribute of the incident transcripts
	# new graphs are maps from gene ID to prefix graphs
	Edge_Flow = {}
	for gname,g in new_graphs.items():
		Edge_Flow[gname] = np.zeros(len(g.edges))
	# process transcript, add abundance
	for old_g in old_graphs:
		new_g = new_graphs[old_g.GeneID]
		# get the list of transcripts
		trans = sum([e.IncidentTranscripts for e in old_g.vEdges], [])
		trans = set(trans)
		# project each transcripts to the edges in prefix graph (new_graphs)
		for tname in trans:
			if not (tname in Exp):
				continue
			old_edges = [e.ID for e in old_g.vEdges if tname in e.IncidentTranscripts]
			old_nodes = [old_g.vEdges[e_id].Node_1 for e_id in old_edges] + [old_g.vEdges[old_edges[-1]].Node_2]
			assert(old_nodes[0] == old_g.vNodes[0].ID and old_nodes[-1] == old_g.vNodes[-1].ID)
			new_nodes, new_edges = new_g.walk(old_nodes)
			# add the expression to the corresponding edge flow in prefix graph
			Edge_Flow[old_g.GeneID][new_edges] += Exp[tname]
	return Edge_Flow


def ReadSalmonQuant(filename):
	Exp = {}
	fp = open(filename, 'r')
	linecount = 0
	col_len = -1
	col_reads = -1
	for line in fp:
		linecount += 1
		if linecount == 1:
			strs = line.strip().split("\t")
			col_len = strs.index("EffectiveLength")
			col_reads = strs.index("NumReads")
			continue
		strs = line.strip().split("\t")
		Exp[strs[0]] = float(strs[col_reads]) / float(strs[col_len])
	fp.close()
	return Exp


def ReadKallistoQuant(filename):
	Exp = {}
	fp = open(filename, 'r')
	linecount = 0
	col_len = -1
	col_reads = -1
	for line in fp:
		linecount += 1
		if linecount == 1:
			strs = line.strip().split("\t")
			col_len = strs.index("eff_length")
			col_reads = strs.index("est_counts")
			continue
		strs = line.strip().split("\t")
		Exp[strs[0]] = float(strs[col_reads]) / float(strs[col_len])
	fp.close()
	return Exp


def ReadGroundTruth(filename):
	Exp = {}
	fp = open(filename, 'r')
	col_exp = -1
	for line in fp:
		if line[0] == '#':
			strs = line[1:].strip().split("\t")
			col_exp = strs.index("CopyNumber")
			continue
		strs = line.strip().split("\t")
		Exp[strs[0]] = float(strs[col_exp])
	fp.close()
	return Exp


def WriteEdgeFlow_new(outputfile, new_Edge_Flow, new_graphs, NameOrder):
	assert(len(new_graphs) == len(Edge_Flow))
	fp = open(outputfile, 'w')
	for gname in NameOrder:
		assert(gname in new_graphs)
		g = new_graphs[gname]
		flows = new_Edge_Flow[gname]
		assert(len(flows) == len(g.edges))
		for j in range(len(flows)):
			e = g.edges[j]
			fp.write("{}\t{}\t{}\t{}\t{}\n".format(g._name, j, flows[j], e[0], e[1]))
	fp.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python GetFlow_quantifier.py <graph_prefix> <salmon_quant> <output_pickle>")
	else:
		graph_prefix = sys.argv[1]
		salmon_quant = sys.argv[2]
		output_pickle = sys.argv[3]

		if not Path(output_pickle).exists():
			old_graphs = ReadGraphFile(graph_prefix + "_graph_fragstart.txt")
			new_graphs, eq_classes = load_data(graph_prefix)
			Exp = ReadSalmonQuant(salmon_quant)
			edge_flow = ProcessEdgeFlow_newgraph(old_graphs, new_graphs, Exp)
			pickle.dump( edge_flow, open(output_pickle, 'wb') )
