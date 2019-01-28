#!/bin/python

import sys
import numpy as np
from GeneGraph import *
from utils import *
from TranscriptClass import *
from trie_conversion import *


def ProcessNodeFlow(TransPaths, Graphs, GraphNameIndex, Exp):
	# initialize Node_Flow for each graph
	Node_Flow = []
	for g in Graphs:
		Node_Flow.append( np.zeros(len(g.vNodes)) )
	# process each transcript
	for info in TransPaths:
		(gname,tname,edges) = info
		g = Graphs[GraphNameIndex[gname]]
		if tname in Exp:
			this_flow = Exp[tname]
			nodes = np.array([g.vEdges[e].Node_1 for e in edges] + [g.vEdges[edges[-1]].Node_2])
			Node_Flow[GraphNameIndex[gname]][nodes] += this_flow
		else:
			continue
	return Node_Flow


def ProcessEdgeFlow(TransPaths, Graphs, GraphNameIndex, Exp):
	# initialize Edge_Flow for each graph
	Edge_Flow = []
	for g in Graphs:
		Edge_Flow.append( np.zeros(len(g.vEdges)) )
	# process each transcript
	for info in TransPaths:
		(gname,tname,edges) = info
		g = Graphs[GraphNameIndex[gname]]
		if tname in Exp:
			this_flow = Exp[tname]
			Edge_Flow[GraphNameIndex[gname]][edges] += this_flow
		else:
			continue
	return Edge_Flow


def ProcessEdgeFlow_newgraph(TransNodes, new_graphs, Exp):
	Edge_Flow = {}
	for gname,g in new_graphs.items():
		Edge_Flow[gname] = np.zeros(len(g.edges))
	# process transcript, add abundance
	for info in TransNodes:
		(gname,tname,old_nodes) = info
		g = new_graphs[gname]
		new_nodes, new_edges = g.walk(old_nodes)
		if tname in Exp:
			this_flow = Exp[tname]
			Edge_Flow[gname][new_edges] += this_flow
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


def WriteNodeFlow(outputfile, Node_Flow, Graphs):
	assert(len(Graphs) == len(Node_Flow))
	fp = open(outputfile, 'w')
	fp.write("# GeneName\tNodeID\tFlow\tChr\tStartPos\tEndPos\tStrand\n")
	for i in range(len(Graphs)):
		g = Graphs[i]
		flows = Node_Flow[i]
		assert(len(flows) == len(g.vNodes))
		for j in range(len(flows)):
			v = g.vNodes[j]
			fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(g.GeneID, j, flows[j], v.Chr, v.StartPos, v.EndPos, v.Strand))
	fp.close()


def WriteEdgeFlow(outputfile, Edge_Flow, Graphs):
	assert(len(Graphs) == len(Edge_Flow))
	fp = open(outputfile, 'w')
	fp.write("# GeneName\tEdgeID\tFlow\tNode_1\tNode_2\n")
	for i in range(len(Graphs)):
		g = Graphs[i]
		flows = Edge_Flow[i]
		assert(len(flows) == len(g.vEdges))
		for j in range(len(flows)):
			e = g.vEdges[j]
			fp.write("{}\t{}\t{}\t{}\t{}\n".format(g.GeneID, j, flows[j], e.Node_1, e.Node_2))
	fp.close()


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


# def WriteEntropy(outputfile, Exp, GeneTransMap):
# 	fp = open(outputfile, 'w')
# 	fp.write("# GeneID\tEntropy\tNormalizedEntropy\n")
# 	for g,trans in GeneTransMap.items():
# 		this_exp = np.array([Exp[t] for t in trans if t in Exp])
# 		if len(this_exp) <= 1 or np.sum(this_exp) < 1e-6:
# 			fp.write("{}\t0\t0\n".format(g))
# 		else:
# 			this_exp[np.where(this_exp == 0)] = 1e-10
# 			this_exp /= np.sum(this_exp)
# 			entropy = -np.dot(this_exp, np.log(this_exp))
# 			max_potential = np.log(len(this_exp))
# 			fp.write("{}\t{}\t{}\n".format(g, entropy, entropy / max_potential))
# 	fp.close()


def WriteEntropy(outputfile, Edge_Flow, Graphs):
	fp = open(outputfile, 'w')
	fp.write("# GeneID\tEntropy\n")
	for i in range(len(Graphs)):
		g = Graphs[i]
		this_flow = Edge_Flow[i]
		entropy = 0
		for v in g.vNodes:
			involved_flow = np.array([this_flow[idx_e] for idx_e in v.OutEdges])
			if len(involved_flow) == 0 or np.sum(involved_flow) == 0:
				continue
			involved_flow[np.where(involved_flow == 0)] += 1e-20
			involved_flow /= np.sum(involved_flow)
			entropy -= np.dot(involved_flow, np.log(involved_flow))
		fp.write("{}\t{}\n".format(g.GeneID, entropy))
	fp.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python GetFlow_NodeEdge.py <mode 0 (ground truth) or 1 (salmon) or 2 (kallisto)> <gtf_file> <graph_file> <prefixtrie> <seq_file> <exp_file> <output_prefix> (<assembly_gtf>)")
	else:
		gtf_file = sys.argv[2]
		graph_file = sys.argv[3]
		prefixtrie = sys.argv[4]
		seq_file = sys.argv[5]
		exp_file = sys.argv[6]
		output_prefix = sys.argv[7]

		if sys.argv[1] == "0":
			Exp = ReadGroundTruth(exp_file)
		elif sys.argv[1] == "1":
			Exp = ReadSalmonQuant(exp_file)
		elif sys.argv[1] == "2":
			Exp = ReadKallistoQuant(exp_file)
		else:
			print("Input mode does not exist")
			sys.exit()

		Transcripts = ReadGTF(gtf_file)
		[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)

		Graphs = ReadGraphFile(graph_file)
		GraphNameIndex = {Graphs[i].GeneID:i for i in range(len(Graphs))}

		new_graphs, eq_classes = load_data(prefixtrie)

		TransEdges = ReadTransEdges(seq_file, Graphs, GraphNameIndex, TransGeneMap)
		TransNodes = ConvertEdges2Nodes(TransEdges, Graphs, GraphNameIndex)
		if len(sys.argv) > 8:
			assembly_gtf = sys.argv[-1]
			asm_transcripts = ReadGTF(assembly_gtf)
			newTransPaths = PathsFromAssembly(asm_transcripts, Graphs)
			TransEdges += newTransPaths

		Node_Flow = ProcessNodeFlow(TransEdges, Graphs, GraphNameIndex, Exp)
		Edge_Flow = ProcessEdgeFlow(TransEdges, Graphs, GraphNameIndex, Exp)

		new_Edge_Flow = ProcessEdgeFlow_newgraph(TransNodes, new_graphs, Exp)

		WriteNodeFlow(output_prefix + "_node_flow.txt", Node_Flow, Graphs)
		WriteEdgeFlow(output_prefix + "_edge_flow.txt", Edge_Flow, Graphs)
		WriteEdgeFlow_new(output_prefix + "_new_edge_flow.txt", new_Edge_Flow, new_graphs, [g.GeneID for g in Graphs])
		# WriteEntropy(output_prefix + "_entropy.txt", Exp, GeneTransMap)
		WriteEntropy(output_prefix + "_entropy.txt", Edge_Flow, Graphs)