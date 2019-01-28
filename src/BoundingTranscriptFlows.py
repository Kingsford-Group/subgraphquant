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


def GetTranscriptBounds(old_graphs, graphs, Opts):
	bounds = {}
	with tqdm.tqdm(total=len(old_graphs)) as pbar:
		for i in range(len(old_graphs)):
			old_g = old_graphs[i]
			g = graphs[old_g.GeneID]
			opt = Opts[i]
			# get a list of transcripts
			tnames = sum([e.IncidentTranscripts for e in old_g.vEdges], [])
			tnames = list(set(tnames))
			# construct flow graph
			fg = FlowGraph(g.edges, opt.results, 0, len(g.nodes)-1)
			# sum of s-out flow
			s_sum = sum([opt.results[j] for j in range(len(opt.results)) if g.edges[j][0] == 0])
			for tn in tnames:
				old_nodes = [e.Node_1 for e in old_g.vEdges if tn in e.IncidentTranscripts] + [old_g.vNodes[-1].ID]
				new_nodes, new_edges = g.walk(old_nodes)
				qlist = [[x] for x in new_edges]
				lb = fg.split_block_flow(qlist)
				ub = fg.block_flow(qlist)
				bounds[tn] = (lb, ub, s_sum)
			pbar.update(1)
	return bounds


def WriteTranscriptBounds(bounds, outputfile):
	fp = open(outputfile, 'w')
	fp.write("# Name\tlb\tub\tsum_gene_flow\n")
	for t,v in bounds.items():
		fp.write("{}\t{}\t{}\t{}\n".format(t, v[0], v[1], v[2]))
	fp.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python BoundingTranscriptFlows.py <graph_prefix> <optobject> <outputfile>")
	else:
		graph_prefix = sys.argv[1]
		optobject = sys.argv[2]
		outputfile = sys.argv[3]

		if not Path(outputfile).exists():
			old_graphs = ReadGraphFile(graph_prefix + "_graph_fragstart.txt")
			graphs, eq_classes = load_data(graph_prefix)
			Opts = pickle.load(open(optobject, 'rb'))
			assert(len(old_graphs) == len(Opts) and len(old_graphs) == len(graphs))
			assert( np.all([old_graphs[i].GeneID == Opts[i].gid for i in range(len(old_graphs))]) )

			bounds = GetTranscriptBounds(old_graphs, graphs, Opts)
			WriteTranscriptBounds(bounds, outputfile)