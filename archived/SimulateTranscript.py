#!/bin/python

import sys
import numpy as np
from GeneGraph import *
from utils import *


def RandomTraverseGraph(g, numpaths = 1):
	assert(len(g.vNodes[0].OutEdges) * len(g.vNodes[-1].InEdges) >= numpaths)
	simu_paths = []
	# random select branch
	count_trials = 0
	while len(simu_paths) < numpaths:
		this_path = [g.vNodes[0].OutEdges[np.random.randint(0, len(g.vNodes[0].OutEdges))]]
		while g.vEdges[this_path[-1]].Node_2 != g.vNodes[-1].ID:
			idx_v = g.vEdges[this_path[-1]].Node_2
			if len(g.vNodes[idx_v].OutEdges) <= 0:
				print([len(g.vNodes), idx_v, g.vNodes[idx_v].OutEdges])
			this_path.append( g.vNodes[idx_v].OutEdges[np.random.randint(0, len(g.vNodes[idx_v].OutEdges))] )
		if not (this_path) in simu_paths:
			simu_paths.append( (this_path) )
		count_trials += 1
		if count_trials > numpaths * 5:
			break
	if len(simu_paths) == numpaths:
		# convert each path from tuple to list
		simu_paths = [list(x) for x in simu_paths]
	else:
		# enumerating all paths
		all_paths = [[e.ID] for e in g.vEdges if e.Node_1 == 0]
		ending_properly = (np.sum([g.vEdges[x[-1]].Node_2 == g.vNodes[-1].ID for x in all_paths]) == len(all_paths))
		while not ending_properly:
			extended_all_paths = []
			for i in range(len(all_paths)):
				if g.vEdges[all_paths[i][-1]].Node_2 != g.vNodes[-1].ID:
					this_path = all_paths[i]
					# extending this path
					v = g.vEdges[this_path[-1]].Node_2
					for e in g.vNodes[v].OutEdges:
						extended_all_paths.append( this_path + [e] )
				else:
					extended_all_paths.append( all_paths[i] )
			all_paths = extended_all_paths
			ending_properly = (np.sum([g.vEdges[x[-1]].Node_2 == g.vNodes[-1].ID for x in all_paths]) == len(all_paths))
		index = np.arange(len(all_paths))
		np.random.shuffle(index)
		simu_paths = [all_paths[o] for o in index[:numpaths]]
	return simu_paths


def SimulateTranscripts(Graphs, multiplier = 2):
	print("Simulating new paths from gene graphs")
	TransPaths = [] # (geneID, transID, path)
	for g in Graphs:
		trans_existing = sum([e.IncidentTranscripts for e in g.vEdges], [])
		trans_existing = list(set(trans_existing))
		path_existing = [([e.ID for e in g.vEdges if tname in e.IncidentTranscripts]) for tname in trans_existing]
		if len(g.vNodes[0].OutEdges) * len(g.vNodes[-1].InEdges) == 1:
			num_path_simulate = 1
		else:
			num_path_simulate = np.random.randint(1, min(multiplier * len(trans_existing), len(g.vNodes[0].OutEdges) * len(g.vNodes[-1].InEdges)) )
		path_new_tmp = RandomTraverseGraph(g, num_path_simulate)
		# convert to tuple
		path_existing = [tuple(p) for p in path_existing]
		path_new_tmp = [tuple(p) for p in path_new_tmp]
		path_new = []
		for p in path_new_tmp:
			if not (p in path_existing):
				path_new.append( p )
		# convert back to list
		path_existing = [list(p) for p in path_existing]
		path_new = [list(p) for p in path_new]
		for i in range(len(path_existing)):
			tname = trans_existing[i]
			p = path_existing[i]
			TransPaths.append( (g.GeneID, tname, p) )
		for i in range(len(path_new)):
			tname = g.GeneID + "_simu_" + str(i)
			p = path_new[i]
			TransPaths.append( (g.GeneID, tname, p) )
	print("Done. There are {} paths in total.".format(len(TransPaths)))
	return TransPaths


def EstimatingExpDistribution(salmonquant):
	SalmonCopyNumber = {}
	fp = open(salmonquant, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		if float(strs[4]) > 0:
			SalmonCopyNumber[strs[0]] = float(strs[4]) / float(strs[2])
	fp.close()
	LogCopyNumber = np.log(np.array(list(SalmonCopyNumber.values())))
	logExp_mean = np.mean(LogCopyNumber)
	logExp_std = np.std(LogCopyNumber)
	return SalmonCopyNumber, logExp_mean, logExp_std


def SimulateExpression(TransPaths, Graphs, GraphNameIndex, logExp_mean, logExp_std, SalmonCopyNumber = None):
	LogCopyNumber = np.random.normal(loc=logExp_mean, scale=logExp_std, size=len(TransPaths))
	CopyNumber = np.array([np.exp(x) for x in LogCopyNumber])
	ReadNumber = []
	for i in range(len(TransPaths)):
		(gname, tname, path) = TransPaths[i]
		g = Graphs[GraphNameIndex[gname]]
		if not (SalmonCopyNumber is None) and (tname in SalmonCopyNumber):
			CopyNumber[i] = SalmonCopyNumber[tname]
		nodes = [g.vEdges[e].Node_1 for e in path] + [g.vEdges[path[-1]].Node_2]
		assert(nodes[0] == 0 and nodes[-1] == g.vNodes[-1].ID)
		length = np.sum([g.vNodes[idx_v].Length for idx_v in nodes])
		copy_number = CopyNumber[i]
		ReadNumber.append( copy_number * length )
	return CopyNumber, ReadNumber


def WriteSequence(outputfile, TransPaths, Graphs, GraphNameIndex, Genome):
	fp = open(outputfile, 'w')
	for i in range(len(TransPaths)):
		(gname, tname, path) = TransPaths[i]
		g = Graphs[GraphNameIndex[gname]]
		nodes = [g.vEdges[e].Node_1 for e in path] + [g.vEdges[path[-1]].Node_2]
		node_str = ""
		seq = ""
		for idx_v in nodes:
			node_str += "(" + str(idx_v) +","+ str(g.vNodes[idx_v].Chr) +","+ str(g.vNodes[idx_v].StartPos) +","+ str(g.vNodes[idx_v].EndPos) +","+ str(g.vNodes[idx_v].Strand) +")"
			if g.vNodes[idx_v].Strand:
				seq += Genome[g.vNodes[idx_v].Chr][g.vNodes[idx_v].StartPos:g.vNodes[idx_v].EndPos]
			else:
				seq += ReverseComplement( Genome[g.vNodes[idx_v].Chr][g.vNodes[idx_v].StartPos:g.vNodes[idx_v].EndPos] )
		assert(len(seq) == np.sum([g.vNodes[idx_v].Length for idx_v in nodes]) - 2)
		# write transcript name and sequence
		fp.write(">" + tname + "\n")
		count = 0
		while count < len(seq):
			fp.write(seq[count:min(count+70, len(seq))] + "\n")
			count += 70
	fp.close()


def WriteSequenceHeader(outputfile, TransPaths, Graphs, GraphNameIndex):
	fp = open(outputfile, 'w')
	for i in range(len(TransPaths)):
		(gname, tname, path) = TransPaths[i]
		g = Graphs[GraphNameIndex[gname]]
		nodes = [g.vEdges[e].Node_1 for e in path] + [g.vEdges[path[-1]].Node_2]
		node_str = ""
		for idx_v in nodes:
			node_str += "(" + str(idx_v) +","+ str(g.vNodes[idx_v].Chr) +","+ str(g.vNodes[idx_v].StartPos) +","+ str(g.vNodes[idx_v].EndPos) +","+ str(g.vNodes[idx_v].Strand) +")"
		# write transcript name and sequence
		fp.write(">" + tname + " " + node_str + "\n")
	fp.close()


def WriteExpression(outputfile, TransPaths, CopyNumber, ReadNumber):
	fp = open(outputfile, 'w')
	fp.write("# Name\tLength\tCopyNumber\tNumReads\n")
	for i in range(len(TransPaths)):
		(gname, tname, path) = TransPaths[i]
		fp.write("{}\t{}\t{}\t{}\n".format(tname, int(ReadNumber[i]/CopyNumber[i]), CopyNumber[i], ReadNumber[i]))
	fp.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python <graph_file> <genome_file> <salmonquant_example> <output_prefix>")
	else:
		graph_file = sys.argv[1]
		genome_file = sys.argv[2]
		salmonquant_example = sys.argv[3]
		output_prefix = sys.argv[4]

		Graphs = ReadGraphFile(graph_file)
		GraphNameIndex = {Graphs[i].GeneID:i for i in range(len(Graphs))}
		TransPaths = SimulateTranscripts(Graphs)

		SalmonCopyNumber, logExp_mean, logExp_std = EstimatingExpDistribution(salmonquant_example)
		CopyNumber, ReadNumber = SimulateExpression(TransPaths, Graphs, GraphNameIndex, logExp_mean, logExp_std, SalmonCopyNumber)

		Genome = ReadGenome(genome_file)
		WriteSequence(output_prefix+"_sequence.fa", TransPaths, Graphs, GraphNameIndex, Genome)
		WriteSequenceHeader(output_prefix+"_sequence_header.fa", TransPaths, Graphs, GraphNameIndex)
		WriteExpression(output_prefix+"_expression.txt", TransPaths, CopyNumber, ReadNumber)