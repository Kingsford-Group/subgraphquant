#!/bin/python

import sys
import numpy as np
import pysam
import scipy.optimize
from TranscriptClass import *
from TranscriptAlignmentClass import *
from Optimize import *
import cvxopt
import tqdm
import copy


Transcripts = ReadGTF("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf")
[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)

Graphs = ReadGraphFile("/home/congm1/savanna/savannacong33/Code/graphsalmon/all_graph_fragstart.txt")
GraphNameIndex = {Graphs[i].GeneID:i for i in range(len(Graphs))}

eqclass = ReadEquivalentClass("/home/congm1/savanna/savannacong33/Code/graphsalmon/scripts/eq_classes.txt")

# read removed transcripts and find corresponding genes
gene_with_del = []
fp = open("/home/congm1/savanna/savannacong33/GraphSalmonSimu/simu_Full_GEU_500_removelist.txt", 'r')
for line in fp:
	gname = TransGeneMap[line.strip().split()[0]]
	if gname in GraphNameIndex:
		gene_with_del.append( GraphNameIndex[gname] )
	else:
		print("All transcripts of {} are deleted.".format(gname))
fp.close()
print("There are {} genes with simulated deletion.".format(len(gene_with_del)))

salmon_flows = ReadFlowFromSalmon("/home/congm1/savanna/savannacong33/GraphSalmonSimu/salmon/Full_GEU_500/quant.sf", Graphs, GraphNameIndex, TransGeneMap)

count_worse = 0
for i in gene_with_del:
		try:
			g = Graphs[i]
			g, normalized_sumbase, sum_reads = optimizeflow_read(g, eqclass)
		except:
			print("error in gene "+Graphs[i].GeneID)
		f_s = salmon_flows[i]
		f_opt = np.array([e.Flow for e in g.vEdges])
		opt = OptimizationObject_read(g, eqclass)
		print("{}: optimized objective = {}\tsalmon objective = {}".format(g.GeneID, opt.objective(f_opt), opt.objective(f_s)))
		if opt.objective(f_opt) > opt.objective(f_s):
			count_worse += 1
