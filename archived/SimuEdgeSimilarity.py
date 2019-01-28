#!/bin/python

import sys
import numpy as np
from GeneGraph import *
from utils import *
from TranscriptClass import *


def CalculateJaccardDistance(info1, info2):
	(gid1, tname1, edges1) = info1
	(gid2, tname2, edges2) = info2
	assert(gid1 == gid2)
	union = set(edges1) | set(edges2)
	intercept = set(edges1) & set(edges2)
	return 1 - 1.0 * len(intercept) / len(union)


def WriteNewPathSimilarity(outputfile, TransPaths, TransGeneMap):
	fp = open(outputfile, 'w')
	fp.write("# GeneID\tName\tJaccardDistance\tCorrespondingKnownTranscript\n")
	TransPaths.sort(key = lambda x:x[0])
	paths = [TransPaths[0]]
	for info in TransPaths[1:]:
		if info[0] != paths[-1][0]:
			# for the last gene, calculate the similarity of novel paths to the existing paths
			paths_existing = [p for p in paths if p[1] in TransGeneMap]
			paths_new = [p for p in paths if not (p[1] in TransGeneMap)]
			for p1 in paths_new:
				min_jaccard = 1
				min_existing = None
				for p2 in paths_existing:
					jac = CalculateJaccardDistance(p1, p2)
					if jac < min_jaccard:
						min_jaccard = jac
						min_existing = p2[1]
				fp.write("{}\t{}\t{}\t{}\n".format(p1[0], p1[1], min_jaccard, min_existing))
			paths = [info]
		else:
			paths.append(info)
	# for the last group
	paths_existing = [p for p in paths if p[1] in TransGeneMap]
	paths_new = [p for p in paths if not (p[1] in TransGeneMap)]
	for p1 in paths_new:
		min_jaccard = 1
		min_existing = None
		for p2 in paths_existing:
			jac = CalculateJaccardDistance(p1, p2)
			if jac < min_jaccard:
				min_jaccard = jac
				min_existing = p2[1]
		fp.write("{}\t{}\t{}\t{}\n".format(p1[0], p1[1], min_jaccard, min_existing))
	fp.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python SimuEdgeSimilarity.py <gtf_existing_file> <graphfile> <sequence_header> <outputfile>")
	else:
		gtf_existing_file = sys.argv[1]
		graphfile = sys.argv[2]
		sequence_header = sys.argv[3]
		outputfile = sys.argv[4]

		Transcripts = ReadGTF(gtf_existing_file)
		[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)

		Graphs = ReadGraphFile(graphfile)
		GraphNameIndex = {Graphs[i].GeneID:i for i in range(len(Graphs))}

		TransPaths = ReadTransPaths(sequence_header, Graphs, GraphNameIndex, TransGeneMap)

		WriteNewPathSimilarity(outputfile, TransPaths, TransGeneMap)