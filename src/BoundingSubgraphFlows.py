#!/bin/python

import sys
import numpy as np
from pathlib import Path
from GeneGraph import *
from IpoptOptimizeClass import *
from utils import *
from trie_conversion import *
from flow_graph import FlowGraph
import pyipopt
import tqdm
import copy
import pickle
import concurrent.futures
import time


def FindUncertainAdjacentEdges(old_g):
	edgepairs = []
	for v in old_g.vNodes:
		if len(v.InEdges) > 1 and len(v.OutEdges) > 1:
			for e1 in v.InEdges:
				for e2 in v.OutEdges:
					edgepairs.append( [e1,e2] )
	return edgepairs


def GetEdgePairBound(g, old_g, res, edgepairs):
	bounds = []
	fg = FlowGraph(g.edges, res, 0, len(g.nodes)-1)
	for p in edgepairs:
		[idx_e1, idx_e2] = p
		edges = [ [old_g.vEdges[idx_e1].Node_1, old_g.vEdges[idx_e1].Node_2], [old_g.vEdges[idx_e2].Node_1, old_g.vEdges[idx_e2].Node_2] ]
		qlist = list(g.match(e) for e in edges)
		qlist = sum(qlist, [])
		lb = fg.diff_flow(qlist)
		ub = fg.split_flow(qlist)
		bounds.append( (lb,ub) )
	return bounds


def WriteEdgePairBound(graphs, old_graphs, results, outputfile):
	fp = open(outputfile, 'w')
	fp.write("# GeneID\toldedge1\toldedge2\tlb\tub\tsum_gene_flow\toldnodes_edge1\toldnodes_edge2\n")
	with tqdm.tqdm(total=len(old_graphs)) as pbar:
		for old_g in old_graphs:
			g = graphs[old_g.GeneID]
			res = results[old_g.GeneID]
			edgepairs = FindUncertainAdjacentEdges(old_g)
			bounds = GetEdgePairBound(g, old_g, res, edgepairs)
			# additional output sum_gene_flow
			s_sum = sum([res[j] for j in range(len(res)) if g.edges[j][0] == 0])
			for i in range(len(edgepairs)):
				[e1,e2] = edgepairs[i]
				(lb,ub) = bounds[i]
				fp.write("{}\t{}\t{}\t{}\t{}\t{}\t({},{})\t({},{})\n".format(old_g.GeneID, e1, e2, lb, ub, s_sum,\
					old_g.vEdges[e1].Node_1, old_g.vEdges[e1].Node_2, old_g.vEdges[e2].Node_1, old_g.vEdges[e2].Node_2))
			pbar.update(1)
	fp.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python BoundingSubgraphFlows.py <graph_prefix> <resobject> <outputfile>")
	else:
		graph_prefix = sys.argv[1]
		resobject = sys.argv[2]
		outputfile = sys.argv[3]

		if not Path(outputfile).exists():
			old_graphs = ReadGraphFile(graph_prefix + "_graph_fragstart.txt")
			graphs, eq_classes = load_data(graph_prefix)
			results = pickle.load(open(resobject, 'rb'))
			assert(len(old_graphs) == len(results) and len(old_graphs) == len(graphs))

			WriteEdgePairBound(graphs, old_graphs, results, outputfile)