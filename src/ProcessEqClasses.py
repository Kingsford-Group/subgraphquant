#!/bin/python

import sys
import numpy as np
import pysam
import tqdm
from TranscriptClass import *
from TranscriptAlignmentClass import *
from GeneGraph import *
from trie_conversion import *


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python ProcessEquivalentClass.py <gtffile> <graphfile> <trans_bamfile> <out_eqclass> <outprefix_trie>")
	else:
		gtffile = sys.argv[1]
		graphfile = sys.argv[2]
		trans_bamfile = sys.argv[3]
		out_eqclass = sys.argv[4]
		outprefix_trie = sys.argv[5]

		Transcripts = ReadGTF(gtffile)
		Graphs = ReadGraphFile(graphfile)
		GraphNameIndex = {Graphs[i].GeneID:i for i in range(len(Graphs))}

		eqclasses = ProcessEquivalentClass(trans_bamfile, Transcripts, Graphs, GraphNameIndex)
		WriteEquivalentClass(out_eqclass, eqclasses)

		work(graphfile, outprefix_trie, outprefix_trie, debug = False)
