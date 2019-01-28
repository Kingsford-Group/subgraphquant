#!/bin/python

import sys
import numpy as np
import pysam
import scipy.optimize
from TranscriptClass import *
from TranscriptAlignmentClass import *
from IpoptOptimizeClass import *
from Optimize import *
from utils import *
from trie_conversion import *
from EffLen_VLMM import *
from flow_graph import FlowGraph
import pyipopt
import cvxopt
import tqdm
import copy
import pickle
import concurrent.futures
import time


def SamplePaths(g, opt):
	fg = FlowGraph(g.edges, opt.results, 0, len(g.nodes)-1)
	tmp_edge_map = {}
	existing_edges = set()
	for e in range(len(g.edges)):
		if e in existing_edges:
			continue
		tmp_edge_list, tmp_log_weights = fg.sample_centered_path(1, e, weighted=True, seed=0)
		tmp_edge_list = tuple(tmp_edge_list[0])
		if not tmp_edge_list in tmp_edge_map:
			tmp_edge_map[tmp_edge_list] = tmp_log_weights[0]
			existing_edges = existing_edges | set(list(tmp_edge_list))
	edge_list = []
	log_weights = []
	for k,v in tmp_edge_map.items():
		edge_list.append(k)
		log_weights.append(v)
	return edge_list, log_weights


def WriteSampleSTPaths(graphs, Opts, old_graphs, old_NameIndex, genomefasta, outputfile):
	genome = ReadGenome(genomefasta)
	fp = open(outputfile, 'w')
	total_paths = 0
	with tqdm.tqdm(total=len(Opts)) as pbar:
		for i in range(len(Opts)):
			opt = Opts[i]
			if np.sum(opt.weights) != 0:
				g = graphs[opt.gid]
				old_g = old_graphs[old_NameIndex[opt.gid]]
				edge_list, log_weights = SamplePaths(g, opt)
				total_paths += len(edge_list)
				assert(len(edge_list) == len(log_weights))
				for j in range(len(edge_list)):
					edges = edge_list[j]
					w = log_weights[j]
					nodes = [g.nodes[g.edges[e][0]][-1] for e in edges] + [g.nodes[g.edges[edges[-1]][1]][-1]]
					node_str =  ""
					seq = ""
					for idx_v in nodes:
						if old_g.vNodes[idx_v].Strand:
							seq += genome[old_g.vNodes[idx_v].Chr][old_g.vNodes[idx_v].StartPos:old_g.vNodes[idx_v].EndPos]
						else:
							seq += ReverseComplement( genome[old_g.vNodes[idx_v].Chr][old_g.vNodes[idx_v].StartPos:old_g.vNodes[idx_v].EndPos] )
						node_str += "(" + str(idx_v) +","+ str(old_g.vNodes[idx_v].Chr) +","+ str(old_g.vNodes[idx_v].StartPos) +","+ str(old_g.vNodes[idx_v].EndPos) +","+ str(old_g.vNodes[idx_v].Strand) +")"
					fp.write(">" + opt.gid + "_sample_" + format(j, '04d') +" "+ str(w) +" "+ node_str + "\n")
					count = 0
					while count < len(seq):
						fp.write(seq[count:min(count+70, len(seq))] + "\n")
						count += 70
			pbar.update(1)
	print("Total paths = {}".format(total_paths))
	fp.close()


folder = "/home/congm1/savanna/savannacong33/GraphSalmonSimu/salmon/nb_Full_GEU_new"
prefix = folder + "/gs"
genomefasta = "/home/congm1/savanna/savannacong33/NCBI/GRCh38.primary_assembly.genome.fa"

corrections = ReadCorrectedLength(prefix + "_corrections_fragstart.dat")
Exp = ReadSalmonQuant(folder + "/quant.sf") 

graphs, eq_classes = load_data(prefix)
old_node_efflens = CalculateOldNodeEffLen(prefix + "_graph_fragstart.txt", Exp, corrections)
efflens = CalculateEdgewiseEffLen(graphs, old_node_efflens)

Gene_Eq_Index = {gname:[] for gname in graphs.keys()}
for i in range(len(eq_classes)):
	for g in [ali.gid for ali in eq_classes[i]]:
		Gene_Eq_Index[g].append(i)
for g in Gene_Eq_Index.keys():
	Gene_Eq_Index[g] = list(set(Gene_Eq_Index[g]))


Names = list(graphs.keys())
Names.sort()
NameIndex = {Names[i]:i for i in range(len(Names))}
Opts_old = [CVXOPTObject_split(graphs[gname], efflens[gname], [eq_classes[j] for j in Gene_Eq_Index[gname]]) for gname in Names]
Opts = [IpoptObject_split(graphs[gname], efflens[gname], [eq_classes[j] for j in Gene_Eq_Index[gname]]) for gname in Names]
Status = {}

pool = concurrent.futures.ProcessPoolExecutor(8)
index = [i for i in range(len(Opts)) if np.sum(Opts[i].weights) > 0]
future_to_opt = [pool.submit(optimizegraph, Opts[i], max_iter = 500, max_cpu_time = 200) for i in index]
with tqdm.tqdm(total=len(index)) as pbar:
	count = 0
	for future in concurrent.futures.as_completed(future_to_opt):
		opt, status = future.result()
		idx = NameIndex[opt.gid]
		assert( len(Opts[idx].edges) == len(opt.results) )
		Opts[idx].results = opt.results
		Status[opt.gid] = status
		count += 1
		pbar.set_postfix(gname=opt.gid, n_var = opt.n_var)
		pbar.update(1)
for i in range(len(Opts)):
	if np.sum(Opts[i].weights) == 0:
		Opts[i].results = np.zeros(Opts[i].n_var)
pickle.dump( Opts, open(prefix + "_opt_ipopt_round0.pkl", 'wb') )
pickle.dump( {opt.gid:opt.results for opt in Opts}, open(prefix + "_result_ipopt_round0.pkl", 'wb') )

# after round 0, sample paths, and update the effective length according to the sample paths
old_graphs = ReadGraphFile(prefix + "_graph_fragstart.txt")
old_NameIndex = {old_graphs[i].GeneID:i for i in range(len(old_graphs))}
WriteSampleSTPaths(graphs, Opts, old_graphs, old_NameIndex, genomefasta, prefix + "_sampled_paths.fa")
# after calculating the new bias, update the efflens
transnodes = ReadTransNodes(prefix + "_sampled_paths.fa", TransGeneMap)
corrections = ReadCorrectedLength(prefix + "_sampled_paths_corrections.dat")
log_weights = ReadTransLogWeights(prefix + "_sampled_paths.fa")
efflens = UpdateEdgewiseEfflen(graphs, old_graphs, corrections, transnodes, log_weights, efflens)
for opt in Opts:
	opt.efflens = efflens[opt.gid]

r = 1
while True:
	diff = 0
	# re-assign eqclass based on flow
	eq_classes = AdjustAssignment_ipopt(Opts, NameIndex, eq_classes)
	Opts = UpdateOptWeight(Opts, NameIndex, eq_classes, Gene_Eq_Index)
	# re-estimate weights
	future_to_opt = [pool.submit(optimizegraph, Opts[i]) for i in range(len(Opts)) if np.sum(Opts[i].weights) > 0]
	with tqdm.tqdm(total=len(future_to_opt)) as pbar:
		for future in concurrent.futures.as_completed(future_to_opt):
			opt, status = future.result()
			idx = NameIndex[opt.gid]
			assert( len(Opts[idx].edges) == len(opt.results) )
			# calculate difference of flow
			diff += np.sum(np.abs( opt.results - Opts[idx].results ))
			# update new
			Opts[idx].results = opt.results
			Status[opt.gid] = status
			pbar.set_postfix(gname=opt.gid, n_var = opt.n_var)
			pbar.update(1)
	# handle the ones with weight 0
	for i in range(len(Opts)):
		if np.sum(Opts[i].weights) == 0:
			Opts[i].results = np.zeros(Opts[i].n_var)
	pickle.dump( Opts, open(prefix + "_opt_ipopt_round" + str(r) + ".pkl", 'wb') )
	pickle.dump( {opt.gid:opt.results for opt in Opts}, open(prefix + "_result_ipopt_round" + str(r) + ".pkl", 'wb') )
	r += 1


