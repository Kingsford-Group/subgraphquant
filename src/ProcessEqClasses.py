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