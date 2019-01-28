#!/bin/python

import sys
import numpy as np
from GeneGraph import *
from utils import *
from TranscriptClass import *
from trie_conversion import *
from IpoptOptimizeClass import *


def ReadNodeFlow(filename):
	Node_Flow = []
	GeneIDs = []
	tmp_flows = []
	tmp_geneid = ""
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[0] != tmp_geneid:
			if tmp_geneid != "":
				# make s and t having flow 0
				tmp_flows[0] = 0
				tmp_flows[-1] = 0
				Node_Flow.append( np.array(tmp_flows) )
				GeneIDs.append( tmp_geneid )
			tmp_geneid = strs[0]
			tmp_flows = []
		tmp_flows.append( float(strs[2]) )
	fp.close()
	# add the last graph
	if tmp_geneid != "":
		tmp_flows[0] = 0
		tmp_flows[-1] = 0
		Node_Flow.append( np.array(tmp_flows) )
		GeneIDs.append( tmp_geneid )
	return Node_Flow, GeneIDs


def ReadEdgeFlow(filename):
	Edge_Flow = []
	GeneIDs = []
	tmp_flows = []
	tmp_geneid = ""
	# read file
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[0] != tmp_geneid:
			if tmp_geneid != "":
				Edge_Flow.append( np.array(tmp_flows) )
				GeneIDs.append(tmp_geneid)
			tmp_geneid = strs[0]
			tmp_flows = []
		tmp_flows.append( float(strs[2]) )
	fp.close()
	# add the last graph
	if tmp_geneid != "":
		Edge_Flow.append( np.array(tmp_flows) )
		GeneIDs.append(tmp_geneid)
	return Edge_Flow, GeneIDs


def normalize_num_node_flow(Node_Flow, GeneIDs, Graphs, GraphNameIndex, normalizer = 1e6):
	cur_sum = 0
	for i in range(len(Node_Flow)):
		g = Graphs[GraphNameIndex[GeneIDs[i]]]
		s_out_node = [e.Node_2 for e in g.vEdges if e.Node_1 == 0]
		cur_sum += np.sum([Node_Flow[i][idx_v] for idx_v in s_out_node])
	print("current sum = {}".format(cur_sum))
	for nf in Node_Flow:
		nf *= normalizer / cur_sum
	return Node_Flow


def normalize_num_edge_flow(Edge_Flow, GeneIDs, Graphs, GraphNameIndex, normalizer = 1e6):
	cur_sum = 0
	for i in range(len(Edge_Flow)):
		g = Graphs[GraphNameIndex[GeneIDs[i]]]
		cur_sum += np.sum([Edge_Flow[i][e.ID] for e in g.vEdges if e.Node_1 == 0])
	print("current sum = {}".format(cur_sum))
	for ef in Edge_Flow:
		ef *= normalizer / cur_sum
	return Edge_Flow


def is_5p_covered(opt):
	nodes_p5 = [e.Node_2 for e in opt.graph.vEdges if e.Node_1 == 0]
	is_covered = False
	for idx_v in nodes_p5:
		if idx_v in [p[0] for p in opt.PathGroups]:
			is_covered = True
	return is_covered


def is_5p_covered(opt, edge_flow_truth):
	s_out_flow = [(e.ID, e.Node_2, edge_flow_truth[e.ID]) for e in opt.graph.vEdges if e.Node_1 == 0]
	s_out_flow.sort(key = lambda x:x[2], reverse = True)
	max_initial_node = s_out_flow[0][1]
	if max_initial_node in [p[0] for p in opt.PathGroups]:
		return True
	return False


def ReadSimuRatio(filename, TransGeneMap, GeneIDs):
	known_copynum = {g:0 for g in GeneIDs}
	simu_copynum = {g:0 for g in GeneIDs}
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[0] in TransGeneMap:
			known_copynum[TransGeneMap[strs[0]]] += float(strs[2])
		else:
			gname = strs[0].split("_")[0]
			simu_copynum[gname] += float(strs[2])
	fp.close()
	ratio = [simu_copynum[gname] / max(known_copynum[gname] + simu_copynum[gname], 1e-8) for gname in GeneIDs]
	return ratio


def WeightedDistance_unknown2known(truth_exp_file, distance_file, TransGeneMap, GeneIDs):
	known_copynum = {g:[] for g in GeneIDs}
	simu_copynum = {g:[] for g in GeneIDs}
	fp = open(truth_exp_file, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[0] in TransGeneMap:
			known_copynum[TransGeneMap[strs[0]]].append( (strs[0], float(strs[2])) )
		else:
			gname = strs[0].split("_")[0]
			simu_copynum[gname].append( (strs[0], float(strs[2])) )
	fp.close()
	distance = {}
	fp = open(distance_file, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		distance[strs[1]] = float(strs[2])
	fp.close()
	weighted_dist = []
	for gname in GeneIDs:
		this_known = known_copynum[gname]
		this_simu = simu_copynum[gname]
		total_copynum = sum([x[1] for x in this_known]) + sum([x[1] for x in this_simu])
		if total_copynum == 0:
			weighted_dist.append(0)
		else:
			for x in this_simu:
				if not (x[0] in distance):
					print(x[0])
			weighted_dist.append( np.sum([1.0 * x[1]/total_copynum * distance[x[0]] for x in this_simu if x[0] in distance]) )
	fp.close()
	return weighted_dist


def CalculateReadARD(reads_est, reads_truth):
	assert(len(reads_est) == len(reads_truth))
	numer = np.abs(reads_est - reads_truth)
	denom = reads_est + reads_truth
	denom[np.where(denom == 0)[0]] = 1
	ARD =  numer / denom 
	return ARD


def GeneLevelCalculatingFlowARD(flow_est, flow_truth, graphs):
	assert(len(flow_est) == len(flow_truth))
	assert(len(flow_est) == len(graphs))
	GeneLevelARD = []
	for i in range(len(graphs)):
		g = graphs[i]
		assert(len(flow_est[i]) == len(g.vEdges))
		assert(len(flow_truth[i]) == len(flow_est[i]))
		# normalize for each gene that the flow sums to 10
		s_est = np.sum([flow_est[i][e.ID] for e in g.vEdges if e.Node_1 == 0])
		s_truth = np.sum([flow_truth[i][e.ID] for e in g.vEdges if e.Node_1 == 0])
		adj_flow_est = flow_est[i]
		if s_est != 0:
			adj_flow_est *= 10 / s_est
		adj_flow_truth = flow_truth[i]
		if s_truth != 0:
			adj_flow_truth *= 10 / s_truth
		numer = np.abs(adj_flow_est - adj_flow_truth)
		denom = adj_flow_est + adj_flow_truth
		denom[np.where(denom == 0)[0]] = 1
		GeneLevelARD.append( numer / denom )
	return GeneLevelARD


def ReadSalmonGeneReads(filename, TransGeneMap, Names, NameIndex, mode="read"):
	GeneReads = np.zeros(len(Names))
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		numreads = float(strs[4])
		gname = TransGeneMap[strs[0]]
		if mode == "read":
			GeneReads[NameIndex[gname]] += numreads
		elif mode == "molecule":
			GeneReads[NameIndex[gname]] += numreads / float(strs[2])
		else:
			print("invalid mode")
			return
	fp.close()
	return GeneReads


def ReadKallistoGeneReads(filename, TransGeneMap, Names, NameIndex, mode="read"):
	GeneReads = np.zeros(len(Names))
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		numreads = float(strs[3])
		gname = TransGeneMap[strs[0]]
		if mode == "read":
			GeneReads[NameIndex[gname]] += numreads
		elif mode == "molecule":
			GeneReads[NameIndex[gname]] += numreads / float(strs[2])
		else:
			print("invalid mode")
			return
	fp.close()
	return GeneReads


def ReadTruthGeneReads(filename, TransGeneMap, Names, NameIndex, mode="read"):
	GeneReads = np.zeros(len(Names))
	col_read = -1
	col_len = -1
	fp = open(filename, 'r')
	for line in fp:
		strs = line.strip().split("\t")
		if line[0] == '#':
			col_read = strs.index("NumReads")
			col_len = strs.index("Length")
			continue
		if not (strs[0] in TransGeneMap):
			continue
		gname = TransGeneMap[strs[0]]
		if not (gname in NameIndex):
			continue
		if mode == "read":
			GeneReads[NameIndex[gname]] += float(strs[col_read])
		elif mode == "molecule":
			GeneReads[NameIndex[gname]] += float(strs[col_read]) / float(strs[col_len])
		else:
			print("invalid mode")
			return
	fp.close()
	return GeneReads


def ComparingReadsARD(ARD, mode="mean"):
	if mode == "mean":
		return np.mean(ARD)
	elif mode == "median":
		return np.median(ARD)


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python CompareARD.py <mode 0 (salmon) 1 (kallisto)> <gtffile> <graphfile> \
			<truth_exp_file> <truth_edge_file> <quant_exp_file> <quant_edge_file> <distance_file> <prefixtrie> <optobjects> <outputfile>")
	else:
		mode = int(sys.argv[1])
		gtffile = sys.argv[2]
		graphfile = sys.argv[3]
		truth_exp_file = sys.argv[4]
		truth_edge_file = sys.argv[5]
		quant_exp_file = sys.argv[6]
		quant_edge_file = sys.argv[7]
		distance_file = sys.argv[8]
		prefixtrie = sys.argv[9]
		optobjects = sys.argv[10]
		outputfile = sys.argv[11]

		Transcripts = ReadGTF(gtffile)
		[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)
		TransLength = GetTransLength(Transcripts)

		old_graphs = ReadGraphFile(graphfile)
		old_NameIndex = {old_graphs[i].GeneID:i for i in range(len(old_graphs))}

		Reads_truth = ReadTruthGeneReads(truth_exp_file, TransGeneMap, [old_g.GeneID for old_g in old_graphs], old_NameIndex, mode="read")
		Edge_Flow_truth, GeneIDs_truth = ReadEdgeFlow(truth_edge_file)
		SumSFlow = np.array( [ np.sum([Edge_Flow_truth[i][e.ID] for e in old_graphs[i].vEdges if e.Node_1 == 0]) for i in range(len(old_graphs)) ] )
		ratio = ReadSimuRatio(truth_exp_file, TransGeneMap, GeneIDs_truth)
		weighted_dist = WeightedDistance_unknown2known(truth_exp_file, distance_file, TransGeneMap, GeneIDs_truth)
		if mode == 0:
			Reads_quant = ReadSalmonGeneReads(quant_exp_file, TransGeneMap, [old_g.GeneID for old_g in old_graphs], old_NameIndex, mode="read")
		elif mode == 1:
			Reads_quant = ReadKallistoGeneReads(quant_exp_file, TransGeneMap, [old_g.GeneID for old_g in old_graphs], old_NameIndex, mode="read")
		else:
			print("Invalid mode {}".format(mode))
			sys.exit()
		Edge_Flow_quant, GeneIDs_quant = ReadEdgeFlow(quant_edge_file)

		graphs, eq_classes = load_data(prefixtrie)
		Opts = pickle.load(open(optobjects, 'rb'))
		NameIndex = {Opts[i].gid:i for i in range(len(Opts))}
		Edge_Flow_gs = []
		Reads_gs = []
		for i in range(len(GeneIDs_truth)):
			gname = GeneIDs_truth[i]
			g = graphs[gname]
			old_g = old_graphs[old_NameIndex[gname]]
			opt = Opts[NameIndex[gname]]
			this_edge_flow = []
			for e in old_g.vEdges:
				this_edge_flow.append( np.sum(opt.results[g.match([e.Node_1, e.Node_2])]) )
			Edge_Flow_gs.append( np.array(this_edge_flow) )
			Reads_gs.append( np.sum(opt.weights) )

		# gene level flow ard
		FlowARD_quant = GeneLevelCalculatingFlowARD(Edge_Flow_quant, Edge_Flow_truth, old_graphs)
		FlowARD_gs = GeneLevelCalculatingFlowARD(Edge_Flow_gs, Edge_Flow_truth, old_graphs)
		# read ard per gene
		ReadARD_quant = CalculateReadARD(Reads_quant, Reads_truth)
		ReadARD_gs = CalculateReadARD(Reads_gs, Reads_truth)

		# write output file
		fp = open(outputfile, 'w')
		fp.write("# GeneID\tMethod\tMeanFlowARD\tReadARD\tExternalRatio\tweighted_distance\tGS_numreads\tGS_effective_vars\tNumTrans\tTrueSFlow\n")
		quantifier = "Salmon"
		if mode == 1:
			quantifier = "Kallisto"
		for i in range(len(FlowARD_quant)):
			opt = Opts[NameIndex[GeneIDs_truth[i]]]
			gs_numreads = np.sum(opt.weights)
			gs_effective_vars = opt.n_var - opt.n_cons
			numtrans = len(GeneTransMap[GeneIDs_truth[i]])
			true_flow = SumSFlow[i]
			fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(GeneIDs_truth[i], quantifier, np.mean(FlowARD_quant[i]), \
				ReadARD_quant[i], ratio[i], weighted_dist[i], gs_numreads, gs_effective_vars, numtrans, true_flow))
			fp.write("{}\tGraphSalmon\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(GeneIDs_truth[i], np.mean(FlowARD_gs[i]), \
				ReadARD_gs[i], ratio[i], weighted_dist[i], gs_numreads, gs_effective_vars, numtrans, true_flow))
		fp.close()
		# print to screen some comparison
		index = np.array([i for i in range(len(ratio)) if ratio[i] > 0.3])
		print("30% quantifier mean flow = {}".format( np.mean(np.concatenate([FlowARD_quant[i] for i in index])) ))
		print("30% quantifier median flow = {}".format( np.median(np.concatenate([FlowARD_quant[i] for i in index])) ))
		print("30% graph salmon mean flow = {}".format( np.mean(np.concatenate([FlowARD_gs[i] for i in index])) ))
		print("30% graph salmon median flow = {}".format( np.median(np.concatenate([FlowARD_gs[i] for i in index])) ))
		index = np.array([i for i in range(len(ratio)) if ratio[i] > 0.4])
		print("40% quantifier mean flow = {}".format( np.mean(np.concatenate([FlowARD_quant[i] for i in index])) ))
		print("40% quantifier median flow = {}".format( np.median(np.concatenate([FlowARD_quant[i] for i in index])) ))
		print("40% graph salmon mean flow = {}".format( np.mean(np.concatenate([FlowARD_gs[i] for i in index])) ))
		print("40% graph salmon median flow = {}".format( np.median(np.concatenate([FlowARD_gs[i] for i in index])) ))
