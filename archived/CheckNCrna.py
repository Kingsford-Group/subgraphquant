#!/bin/python

import sys
import numpy as np
from TranscriptClass import *


def GetNC_genes_trans(gtffile):
	NCGenes = []
	NCTrans = []
	fp = open(gtffile, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[2] == "gene":
			biotype = GetFeature(line, "gene_type")
			gene_id = GetFeature(line, "gene_id")
			if "ncRNA" in biotype:
				NCGenes.append(gene_id)
		elif strs[2] == "transcript":
			biotype = GetFeature(line, "transcript_type")
			trans_id = GetFeature(line, "transcript_id")
			if "ncRNA" in biotype:
				NCTrans.append(trans_id)
	fp.close()
	return NCGenes, NCTrans


Transcripts = ReadGTF("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf")
[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)
NCGenes, NCTrans = GetNC_genes_trans("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf")

print("total NCGenes = {}\ttotal NCTrans = {}".format(len(NCGenes), len(NCTrans)))

# count NCGenes that contain 1, 2, 3, >4 transcripts
count_trans = {1:0, 2:0, 3:0, 4:0}
for g in NCGenes:
	ntrans = len(GeneTransMap[g])
	if ntrans <= 4:
		count_trans[ntrans] += 1
	else:
		count_trans[4] += 1
print(count_trans)

# count NCTrans that are outside NCGenes
NCGenes_set = set(NCGenes)
count_outside = 0
for t in NCTrans:
	if not (TransGeneMap[t] in NCGenes_set):
		count_outside += 1
print(count_outside)