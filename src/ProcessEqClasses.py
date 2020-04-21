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


def write_mapping_simpleedge_prefixedge(Graphs, prefix_graphs, outprefix_trie, maxfraglen = 300):
	# one simple edge of original splice graph, in theory, should correspond to the sum of abundances of several prefix graph edges.
	# simple edge means either single node or one pair of connected node
	fp = open(outprefix_trie + "_simpleedge_to_path.txt", 'w')
	fp.write("# gene\tsimple_node_list\tprefix_edge_indexes\n")
	for g in Graphs:
		gname = g.GeneID
		new_g = prefix_graphs[gname]
		# collect subpaths in reference transcripts where the sum of length of inner nodes < maxfraglen
		all_trans = list(set(sum([e.IncidentTranscripts for e in g.vEdges], [])))
		all_paths = []
		for t in all_trans:
			node_list = [e.Node_1 for e in g.vEdges if t in e.IncidentTranscripts] + [g.vNodes[-1].ID]
			for i in range(1, len(node_list) - 1):
				for j in range(i+2, len(node_list)):
					inner_length = np.sum([g.vNodes[k].EndPos - g.vNodes[k].StartPos for k in node_list[(i+1):j]])
					if inner_length < maxfraglen:
						all_paths.append( tuple(node_list[i:j]) )
						assert(all_paths[-1][0] != 0 and all_paths[-1][-1] != len(g.vNodes) - 1)
		all_paths = list(set(all_paths))
		# single node
		for v in g.vNodes:
			if v.ID == 0 or v.ID == len(g.vNodes) - 1:
				continue
			new_edges = new_g.match([v.ID])
			fp.write("{}\t{}\t{}\n".format(gname, v.ID, ",".join([str(x) for x in new_edges])))
		for nl in all_paths:
			assert(nl[0] != 0 and nl[-1] != len(g.vNodes) - 1)
			old_nodes = list(nl)
			i = 0
			while i < len(old_nodes) - 1:
				new_edges = new_g.match(old_nodes[i:])
				if new_edges != []:
					break
				i += 1
			if len(new_edges) > 0:
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
		write_mapping_simpleedge_prefixedge(Graphs, prefix_graphs, outprefix_trie)