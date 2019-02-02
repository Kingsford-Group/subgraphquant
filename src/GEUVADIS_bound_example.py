#!/bin/python

import sys
import numpy as np
import pickle
from pathlib import Path
from GeneGraph import *
from trie_conversion import *
from CompareARD import *


def ReadCenterMetadata(filename):
	CenterSampleMap = {}
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[1] in CenterSampleMap:
			CenterSampleMap[strs[1]].append(strs[-1])
		else:
			CenterSampleMap[strs[1]] = [strs[-1]]
	fp.close()
	return CenterSampleMap


def ReadAbundanceGap(filename):
	AbundanceGap = {}
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		lb = float(strs[1])
		ub = float(strs[2])
		if lb < 0:
			lb = 0
		if ub < 0:
			ub = 0
		if lb > ub:
			ub = lb
		AbundanceGap[strs[0]] = (lb, ub)
	fp.close()
	return AbundanceGap


def ReadTPM(filename):
	TPM = {}
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		TPM[strs[0]] = float(strs[3])
	fp.close()
	return TPM


def ReadGraphSalmonAbundance(graphfile, prefixtrie, resultfile):
	old_graphs = ReadGraphFile(graphfile)
	old_NameIndex = {old_graphs[i].GeneID:i for i in range(len(old_graphs))}
	graphs, eq_classes = load_data(prefixtrie)
	results = pickle.load(open(resultfile, 'rb'))
	# convert prefix trie edge to splicing edge and get flow
	Edge_Flow_gs = []
	for i in range(len(old_graphs)):
		old_g = old_graphs[i]
		gname = old_g.GeneID
		g = graphs[gname]
		res = results[gname]
		this_edge_flow = []
		for e in old_g.vEdges:
			this_edge_flow.append( np.sum(res[g.match([e.Node_1, e.Node_2])]) )
		Edge_Flow_gs.append( np.array(this_edge_flow) )
	return Edge_Flow_gs, [old_g.GeneID for old_g in old_graphs]


def CalculateIntervalSimilarity(intervals1, intervals2):
	breaks = list(set([x[0] for x in intervals1+intervals2] + [x[1] for x in intervals1+intervals2]))
	breaks.sort()
	breaks = np.array(breaks)
	dist1 = np.zeros(len(breaks))
	dist2 = np.zeros(len(breaks))
	for e in intervals1:
		idx_s = np.where(breaks >= e[0])[0][0]
		idx_t = np.where(breaks >= e[1])[0][0]
		dist1[idx_s:idx_t] += 1
	for e in intervals2:
		idx_s = np.where(breaks >= e[0])[0][0]
		idx_t = np.where(breaks >= e[1])[0][0]
		dist2[idx_s:idx_t] += 1
	union = 0
	for i in range(len(breaks)-1):
		union += max(dist1[i], dist2[i]) * (breaks[i+1] - breaks[i])
	intersect = 0
	for i in range(len(breaks)-1):
		intersect += min(dist1[i], dist2[i]) * (breaks[i+1] - breaks[i])
	if union == 0:
		return 0
	else:
		return intersect / union


def CalculateDist(intervals1, intervals2):
	breaks = list(set([x[0] for x in intervals1+intervals2] + [x[1] for x in intervals1+intervals2]))
	breaks.sort()
	breaks = np.array(breaks)
	dist1 = np.zeros(len(breaks))
	dist2 = np.zeros(len(breaks))
	for e in intervals1:
		idx_s = np.where(breaks >= e[0])[0][0]
		idx_t = np.where(breaks >= e[1])[0][0]
		dist1[idx_s:idx_t] += 1
	for e in intervals2:
		idx_s = np.where(breaks >= e[0])[0][0]
		idx_t = np.where(breaks >= e[1])[0][0]
		dist2[idx_s:idx_t] += 1
	return breaks, dist1, dist2


def ProcessCenterwiseSimilarity(folder, CenterSampleMap):
	Gaps = {c:[] for c in CenterSampleMap.keys()}
	TransNames = set()
	for c,samples in CenterSampleMap.items():
		for s in samples:
			gapfile = folder + "/salmon_Full_" + s + "/prefixgraph/gs_bound_round0.txt"
			AbundanceGap = ReadAbundanceGap(gapfile)
			Gaps[c].append(AbundanceGap)
			if len(TransNames) == 0:
				TransNames = set(AbundanceGap.keys())
			else:
				assert(TransNames == set(AbundanceGap.keys()))
			print("finish reading {}:{}".format(c, s))
	TransSimilarity = {t:0 for t in TransNames}
	centers = list(CenterSampleMap.keys())
	assert(len(centers) == 2)
	for t in TransNames:
		intervals1 = [ AbundanceGap[t] for AbundanceGap in Gaps[centers[0]] ]
		intervals2 = [ AbundanceGap[t] for AbundanceGap in Gaps[centers[1]] ]
		TransSimilarity[t] = CalculateIntervalSimilarity(intervals1, intervals2)
	return TransSimilarity


def ProcessAllSampleEdgeARD(folder, CenterSampleMap):
	SampleGeneARD = {} # key is (sample, gene) tuple
	GeneIDs = []
	for c,samples in CenterSampleMap.items():
		for s in samples:
			salmonedgefile = folder + "/salmon_Full_" + s + "/prefixgraph/salmon_edge_flow.txt"
			graphfile = folder + "/salmon_Full_" + s + "/prefixgraph/gs_graph_fragstart.txt"
			prefixtrie = folder + "/salmon_Full_" + s + "/prefixgraph/gs"
			resultfile = folder + "/salmon_Full_" + s + "/prefixgraph/gs_result_ipopt_round0.pkl"
			old_graphs = ReadGraphFile(graphfile)
			Edge_Flow_salmon, tmpgene_salmon = ReadEdgeFlow(salmonedgefile)
			Edge_Flow_gs, tmpgene_gs = ReadGraphSalmonAbundance(graphfile, prefixtrie, resultfile)
			assert(len(tmpgene_gs) == len(tmpgene_salmon))
			assert(np.all([tmpgene_gs[i] == tmpgene_salmon[i] for i in range(len(tmpgene_gs))]))
			# calculate ARD
			ARD_pergene = GeneLevelCalculatingFlowARD(Edge_Flow_salmon, Edge_Flow_gs, old_graphs)
			assert(len(ARD_pergene) == len(old_graphs))
			for i in range(len(old_graphs)):
				SampleGeneARD[(s,old_graphs[i].GeneID)] = ARD_pergene[i]
			if len(GeneIDs) == 0:
				GeneIDs = [old_g.GeneID for old_g in old_graphs]
			else:
				assert(np.all([GeneIDs[i] == old_graphs[i].GeneID]))
			print("finish reading {}:{}".format(c, s))
			print("sanity check: size of SampleGeneARD = {}; size of first key = {}; size of second key = {}".format(len(SampleGeneARD), len(set([x[0] for x in SampleGeneARD.keys()])), len(set([x[1] for x in SampleGeneARD.keys()]))))
	# for each gene, calculate mean ARD across all samples
	ARD = {}
	samples = sum(list(CenterSampleMap.values()), [])
	for g in GeneIDs:
		ard_allsample = [SampleGeneARD[(s,g)] for s in samples]
		ard_allsample = np.concatenate(ard_allsample)
		ARD[g] = np.mean(ard_allsample)
	return ARD


def WriteTransSimilarity(outputfile, TransSimilarity):
	fp = open(outputfile, 'w')
	for t,v in TransSimilarity.items():
		fp.write("{}\t{}\n".format(t, v))
	fp.close()


def WriteAllSampleARD(outputfile, ARD):
	fp = open(outputfile, 'w')
	fp.write("# GeneID\tAvgARD\n")
	for g,v in ARD.items():
		fp.write("{}\t{}\n".format(g, v))
	fp.close()


def WriteExample(folder, CenterSampleMap, t, output_example):
	fp = open(output_example, 'w')
	fp.write("# Name\tSample\tCenter\tlb\tub\tTPM\n")
	for c,samples in CenterSampleMap.items():
		for s in samples:
			gapfile = folder + "/salmon_Full_" + s + "/prefixgraph/gs_bound_round0.txt"
			AbundanceGap = ReadAbundanceGap(gapfile)
			quantfile = folder + "/salmon_Full_" + s + "/quant.sf"
			TPM = ReadTPM(quantfile)
			fp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(t, s, c, AbundanceGap[t][0], AbundanceGap[t][1], TPM[t]))
	fp.close()


def WriteExampleDist(folder, CenterSampleMap, t, output_example):
	Gaps = {c:[] for c in CenterSampleMap.keys()}
	fp = open(output_example, 'w')
	fp.write("# Name\tCenter\tposition\tcount\n")
	for c,samples in CenterSampleMap.items():
		for s in samples:
			gapfile = folder + "/salmon_Full_" + s + "/prefixgraph/gs_bound_round0.txt"
			AbundanceGap = ReadAbundanceGap(gapfile)
			Gaps[c].append(AbundanceGap[t])
	centers = list(CenterSampleMap.keys())
	assert(len(centers) == 2)
	intervals1 = Gaps[centers[0]]
	intervals2 = Gaps[centers[1]]
	breaks, dist1, dist2 = CalculateDist(intervals1, intervals2)
	for i in range(len(dist1)-1):
		fp.write("{}\t{}\t{}\t{}\n".format(t, centers[0], breaks[i], dist1[i]))
	for i in range(len(dist2)-1):
		fp.write("{}\t{}\t{}\t{}\n".format(t, centers[1], breaks[i], dist2[i]))
	fp.close()


if __name__=="__main__":
	metafile = "../data/Metadata.txt"
	folder = "GEUVADIS/"
	outputprefix = "GEUVADIS/bounds"
	output_example = outputprefix + "_example"

	if len(sys.argv) == 1:
		print("python GEUVADIS_bound_example.py (<data meta file>) (<GEUVADIS folder>) (<output prefix>)")
		print("Using the following default parameters:")
		print("\tdata meta file = {}".format(metafile))
		print("\tGEUVADIS folder = {}".format(folder))
		print("\toutput prefix = {}".format(outputprefix))
	else:
		metafile = sys.argv[1]
		folder = sys.argv[2]
		outputprefix = sys.argv[3]
		output_example = outputprefix + "_example"

	CenterSampleMap = ReadCenterMetadata(metafile)
	TransSimilarity = ProcessCenterwiseSimilarity(folder, CenterSampleMap)
	WriteTransSimilarity(outputprefix + "_weightedJaccard.txt", TransSimilarity)

	ARD = ProcessAllSampleEdgeARD(folder, CenterSampleMap)
	WriteAllSampleARD(outputprefix + "_meanflowARD.txt", ARD)

	t = "ENST00000483767.5"
	WriteExample(folder, CenterSampleMap, t, output_example + t + ".txt")
	WriteExampleDist(folder, CenterSampleMap, t, output_example + "_dist_" + t + ".txt")
	t = "ENST00000375754.8"
	WriteExample(folder, CenterSampleMap, t, output_example + t + ".txt")
	WriteExampleDist(folder, CenterSampleMap, t, output_example + "_dist_" + t + ".txt")
	t = "ENST00000422247.6"
	WriteExample(folder, CenterSampleMap, t, output_example + t + ".txt")
	WriteExampleDist(folder, CenterSampleMap, t, output_example + "_dist_" + t + ".txt")
	t = "ENST00000265187.4"
	WriteExample(folder, CenterSampleMap, t, output_example + t + ".txt")
	WriteExampleDist(folder, CenterSampleMap, t, output_example + "_dist_" + t + ".txt")
