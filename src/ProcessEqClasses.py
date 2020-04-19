#!/bin/python

import sys
import numpy as np
import pysam
import tqdm
import copy
from TranscriptClass import *
from TranscriptAlignmentClass import *
from GeneGraph import *
from trie_conversion import *


def write_paths_prefixedge(prefix_graphs, outprefix_trie):
	fp = open(outprefix_trie + "_paths_of_trie.txt", 'w')
	fp.write("# gene\tpath\n")
	for gname, g in prefix_graphs.items():
		for e in g.edges:
			path = g.nodes[e[0]] + g.nodes[e[1]]
			path = list(set(path))
			path.sort()
			path = [str(x) for x in path]
			fp.write("{}\t{}\n".format(gname, ",".join(path)) )
	fp.close()


def write_mapping_simpleedge_prefixedge(Graphs, prefix_graphs, outprefix_trie):
	# one simple edge of original splice graph, in theory, should correspond to the sum of abundances of several prefix graph edges.
	# simple edge means either single node or one pair of connected node
	fp = open(outprefix_trie + "_simpleedge_to_path.txt", 'w')
	fp.write("# gene\tsimple_node_list\tprefix_edge_indexes\n")
	for g in Graphs:
		gname = g.GeneID
		new_g = prefix_graphs[gname]
		# single node
		for v in g.vNodes:
			if v.ID == 0 or v.ID == len(g.vNodes) - 1:
				continue
			new_edges = new_g.match([v.ID])
			fp.write("{}\t{}\t{}\n".format(gname, v.ID, ",".join([str(x) for x in new_edges])))
		for e in g.vEdges:
			if e.Node_1 == 0 or e.Node_2 == len(g.vNodes) - 1:
				continue
			old_nodes = [e.Node_1, e.Node_2]
			new_edges = new_g.match(old_nodes)
			fp.write("{}\t{}\t{}\n".format(gname, ",".join([str(x) for x in old_nodes]), ",".join([str(x) for x in new_edges])))
	fp.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python ProcessEquivalentClass.py <gtffile> <graphfile> <salmon_folder> <out_eqclass> <outprefix_trie>")
	else:
		gtffile = sys.argv[1]
		graphfile = sys.argv[2]
		salmon_folder = sys.argv[3]
		out_eqclass = sys.argv[4]
		outprefix_trie = sys.argv[5]

		Transcripts = ReadGTF(gtffile)
		Graphs = ReadGraphFile(graphfile)
		GraphNameIndex = {Graphs[i].GeneID:i for i in range(len(Graphs))}

		eqclasses = ProcessEquivalentClass(salmon_folder, Transcripts, Graphs, GraphNameIndex, out_eqclass)

		work(graphfile, out_eqclass, outprefix_trie, debug = False)

		# write the paths corresponding to each prefix graph edge
		prefix_graphs, eqclasses = load_data(outprefix_trie)
		write_paths_prefixedge(prefix_graphs, outprefix_trie)