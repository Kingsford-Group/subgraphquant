#!/bin/python

import sys
import numpy as np
import pysam
import scipy.optimize
from TranscriptClass import *
from TranscriptAlignmentClass import *
from Optimize import *
from utils import *
import cvxopt
import tqdm
import copy
import pickle
import concurrent.futures
import time

folder = "/home/congm1/savanna/savannacong33/GraphSalmonSimu/salmon/Full_GEU_new"
prefix = folder + "/gs"

Transcripts = ReadGTF("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf")
[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)
TransLength = GetTransLength(Transcripts)

Graphs = ReadGraphFile(prefix + "_graph_fragstart.txt")
GraphNameIndex = {Graphs[i].GeneID:i for i in range(len(Graphs))}

eqclass = ProcessEquivalentClass(folder + "/mapping.bam", Transcripts, Graphs, GraphNameIndex)
WriteEquivalentClass(folder + "/eq_classes.txt", eqclass)
# eqclass = ReadEquivalentClass(folder + "/eq_classes.txt")

Gene_Eq_Index = {}
for i in range(len(eqclass)):
	for g in eqclass[i].GeneIDs:
		if g in Gene_Eq_Index:
			Gene_Eq_Index[g].append(i)
		else:
			Gene_Eq_Index[g] = [i]

# start = time.time()
# Opts = [None] * len(Graphs)
# pool = concurrent.futures.ProcessPoolExecutor(20)
# future_to_opt = [pool.submit(CVXOPToptimize_hyper, g, [eqclass[j] for j in Gene_Eq_Index[g.GeneID]], True) for g in Graphs if g.GeneID in Gene_Eq_Index]
# for future in concurrent.futures.as_completed(future_to_opt):
# 	g, opt = future.result()
# 	idx_g = GraphNameIndex[g.GeneID]
# 	Graphs[idx_g] = g
# 	Opts[idx_g] = opt
# end = time.time()
# print(end - start)

Opts = []
with tqdm.tqdm(total=len(Graphs)) as pbar:
	for i in range(len(Graphs)):
		g = Graphs[i]
		if g.GeneID in Gene_Eq_Index:
			g, opt = CVXOPToptimize_hyper(g, [eqclass[j] for j in Gene_Eq_Index[g.GeneID]], normalize_to_readcount = True)
			Opts.append(opt)
			Graphs[i] = g
		else:
			opt = None
			Opts.append(opt)
		pbar.update(1)
assert(len(Opts) == len(Graphs))
for i in range(len(Graphs)):
	assert(Graphs[i].GeneID == Opts[i].graph.GeneID)
WriteGraphFile(prefix + "_graph_withrawcoverage_round0", Graphs)
pickle.dump( (Graphs, Opts) , open(prefix + "_opt_object_round0.pkl", 'wb'))

# only keep eqclass that are multi-mapped
eqclass = [eq for eq in eqclass if len(eq.GeneIDs) > 1]
for eq in eqclass:
	eq.assignweights(Graphs, GraphNameIndex, Opts)
# update the index of eqclass of each gene
Gene_Eq_Index = {}
for i in range(len(eqclass)):
	for g in eqclass[i].GeneIDs:
		if g in Gene_Eq_Index:
			Gene_Eq_Index[g].append(i)
		else:
			Gene_Eq_Index[g] = [i]

r = 1
diff = 100
while diff > 1:
	diff = 0
	prev_flows = [np.array([e.Flow for e in g.vEdges]) for g in Graphs]
	prev_numreads = np.sum([np.sum(opt.Weights) for opt in Opts])
	# update read assignment
	for i in range(len(Graphs)):
		g = Graphs[i]
		opt = Opts[i]
		if not (opt is None) and (g.GeneID in Gene_Eq_Index):
			opt.update_weights_chaning( [eqclass[j] for j in Gene_Eq_Index[g.GeneID]] )
	current_numreads = np.sum([np.sum(opt.Weights) for opt in Opts])
	# estimate flow with new read assignment
	pool = concurrent.futures.ThreadPoolExecutor(8)
	future_to_opt = [pool.submit(CVXOPToptimize_hyper_cont, Graphs[i], Opts[i], True) for i in range(len(Graphs)) if (Graphs[i].GeneID in Gene_Eq_Index)]
	with tqdm.tqdm(total=len(future_to_opt)) as pbar:
		for future in concurrent.futures.as_completed(future_to_opt):
			g, opt = future.result()
			idx_g = GraphNameIndex[g.GeneID]
			Graphs[idx_g] = g
			Opts[idx_g] = opt
			pbar.update(1)
	# calculate difference of flow
	print("number of optimization failure = {}".format(len([i for i in range(len(Opts)) if Opts[i].Results is None])))
	current_flow = [np.array([e.Flow for e in g.vEdges]) for g in Graphs]
	diff = np.sum([np.sum(np.abs(pf[i] - cf[i])) for i in range(len(prev_flow))])
	WriteGraphFile(prefix + "_graph_withrawcoverage_round" + str(r), Graphs)
	pickle.dump( (Graphs, Opts) , open(prefix + "_opt_object_round"+str(r)+".pkl", 'wb'))
	r += 1


# # write paths with top flows
# Genome = ReadGenome("/home/congm1/savanna/savannacong33/NCBI/GRCh38.primary_assembly.genome.fa")
# WriteTopExpressedPaths("trans_top_flows.fa", [Graphs[i] for i in gene_with_del], Genome)

# TransFlow = ReadTruth("/home/congm1/savanna/savannacong33/SADstandardsimu/salmon/Full_GEU_500/expression_truth.txt", TransLength)
# SalmonTransFlow = ReadSalmon("/home/congm1/savanna/savannacong33/GraphSalmonSimu/salmon/Full_GEU_500/quant.sf")
# normalized_numreads = {Opts[i].graph.GeneID : np.sum(Opts[i].Weights) for i in range(len(Opts))}
# NodeFlow = GetTrueNodeFlow(Graphs, GraphNameIndex, Transcripts, TransFlow, normalized_numreads)
# SalmonNodeFlow = GetTrueNodeFlow(Graphs, GraphNameIndex, Transcripts, SalmonTransFlow, normalized_numreads)
# result = []
# for i in range(len(Graphs)):
# 	g = Graphs[i]
# 	x = np.vstack( (NodeFlow[g.GeneID], g.get_node_flow(), SalmonNodeFlow[g.GeneID]) ).transpose()
# 	result.append( [np.mean(np.abs(x[:,1] - x[:,0]) / np.maximum(x[:,1] + x[:,0], 1e-10)), np.mean(np.abs(x[:,2] - x[:,0]) / np.maximum(x[:,2] + x[:,0], 1e-10))] )
# pickle.dump(result, open("result_round_0.pkl", 'rb'))

# # XXX: calculate bias correction
# corrections = ReadCorrection("trans_top_correction.dat")
# Graphs = UpdateNodeBiasMultiplier("trans_top_flows.fa", corrections, Graphs, GraphNameIndex, damping = 0.1) 
# Opts2 = []
# with tqdm.tqdm(total=len(gene_with_del)) as pbar:
# 	for i in gene_with_del:
# 		try:
# 			g = Graphs[i]
# 			g, opt = optimizeflow_hyper(g, eqclass, normalize_to_readcount = True)
# 			Opts2.append(opt)
# 			Graphs[i] = g
# 		except:
# 			print("error in gene "+Graphs[i].GeneID)
# 		pbar.update(1)

# result2 = []
# for i in gene_with_del:
# 	g = Graphs[i]
# 	x = np.vstack( (NodeFlow[g.GeneID], g.get_node_flow(), SalmonNodeFlow[g.GeneID]) ).transpose()
# 	result2.append( [np.mean(np.abs(x[:,1] - x[:,0]) / np.maximum(x[:,1] + x[:,0], 1e-10)), np.mean(np.abs(x[:,2] - x[:,0]) / np.maximum(x[:,2] + x[:,0], 1e-10))] )
