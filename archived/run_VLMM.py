#!/bin/python

import sys
import numpy as np
import pysam
import scipy.optimize
from TranscriptClass import *
from TranscriptAlignmentClass import *
from Optimize import *
from utils import *
from trie_conversion import *
from EffLen_VLMM import *
import cvxopt
import tqdm
import copy
import pickle
import concurrent.futures
import time

folder = "/home/congm1/savanna/savannacong33/GraphSalmonSimu/salmon/Full_GEU_new"
prefix = folder + "/gs"

corrections = ReadCorrectedLength(prefix + "_corrections_fragstart.dat")
Exp = ReadSalmonQuant(folder + "/quant.sf") 

graphs, eq_classes = load_data("/home/hongyuz/src/gsalmon_fork/sim2/test")
old_node_efflens = CalculateOldNodeEffLen(prefix + "_graph_fragstart.txt", Exp, corrections)
efflens = CalculateEdgewiseEffLen(graphs, old_nodes_efflen)

Gene_Eq_Index = {gname:[] for gname in graphs.keys()}
for i in range(len(eq_classes)):
	for g in [ali.gid for ali in eq_classes[i]]:
		Gene_Eq_Index[g].append(i)
for g in Gene_Eq_Index.keys():
	Gene_Eq_Index[g] = list(set(Gene_Eq_Index[g]))


Opts = [CVXOPTObject_split(g, efflens[gname], [eq_classes[j] for j in Gene_Eq_Index[g._name]]) for gname,g in graphs.items()]
NameIndex = {Opts[i].GeneID:i for i in range(len(Opts))}

with tqdm.tqdm(total=len(graphs)) as pbar:
	for i in range(len(Opts)):
		pbar.set_postfix(gname=Opts[i].GeneID)
		pbar.update(1)
		if len(Opts[i].Edges) > 10:
			opt = CVXOPToptimize_split(Opts[i], True)
		if not (opt is None):
			Opts[i] = opt


pool = concurrent.futures.ThreadPoolExecutor(4)
future_to_opt = [pool.submit(CVXOPToptimize_split, Opts[i], True) for i in range(len(Opts)) if len(Gene_Eq_Index[Opts[i].GeneID]) > 0]
with tqdm.tqdm(total=len(future_to_opt)) as pbar:
	for future in concurrent.futures.as_completed(future_to_opt):
		try:
			opt = future.result(timeout = 10)
			idx = NameIndex[opt.GeneID]
			Opts[idx] = opt
		except concurrent.futures.TimeoutError:
			print("Timeout.")
		pbar.update(1)
