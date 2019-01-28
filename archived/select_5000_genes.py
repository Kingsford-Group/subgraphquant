#!/bin/python

import sys
import numpy as np
from trie_conversion import *
from TranscriptClass import *
from utils import *


def Genes_with_multialign(eq_classes):
	genes_multi = set()
	for eq in eq_classes:
		if len(eq) > 1:
			gids = set([ali.gid for ali in eq])
			if len(gids) > 1:
				for gid in gids:
					genes_multi.add(gid)
	print("There are {} genes that contain multi-mapping.".format(len(genes_multi)))
	return genes_multi


def Select_uniqalign_genes(graphs, genes_multi, TransNodes, num_selection = 5000, path_thresh = 5):
	graph_path_summary = {gname:0 for gname,g in graphs.items()}
	for info in TransNodes:
		(gname,tname,nodes) = info
		assert(gname in info)
		graph_path_summary[gname] += 1
	# selection
	selection = []
	for gname,g in graphs.items():
		if gname in genes_multi or graph_path_summary[gname] < path_thresh:
			continue
		selection.append(gname)
		if len(selection) > num_selection:
			break
	print("Final number of selection = " + str(len(selection)))
	return set(selection)


def UpdateGTFfile(ingtffile, outgtffile, selection):
	fpin = open(ingtffile, 'r')
	fpout = open(outgtffile, 'w')
	for line in fpin:
		if line[0] == '#':
			continue
		gene_id = GetFeature(line, "gene_id")
		if gene_id in selection:
			fpout.write(line)
	fpin.close()
	fpout.close()
	print("Finish updating GTF file.")


def UpdateFastafile(infasta, outfasta, selection, TransGeneMap):
	fpin = open(infasta, 'r')
	fpout = open(outfasta, 'w')
	to_write = False
	for line in fpin:
		if line[0] == '>':
			tname = line.split()[0][1:]
			gname = ""
			if tname in TransGeneMap:
				gname = TransGeneMap[tname]
			else:
				if "ENSG" not in tname:
					print(tname)
				assert("ENSG" in tname)
				gname = tname.split("_")[0]
			if gname in selection:
				to_write = True
		if to_write:
			fpout.write(line)
	fpin.close()
	fpout.close()
	print("Finish updating fasta file.")


def UpdateSequenceHeader(inheader, outheader, selection, TransGeneMap):
	fpin = open(inheader, 'r')
	fpout = open(outheader, 'w')
	for line in fpin:
		if line[0] == '>':
			tname = line.split()[0][1:]
			gname = ""
			if tname in TransGeneMap:
				gname = TransGeneMap[tname]
			else:
				assert("ENSG" in tname)
				gname = tname.split("_")[0]
			if gname in selection:
				fpout.write(line)
	fpin.close()
	fpout.close()
	print("Finish updating sequence header.")


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
	print("logExp_mean = {}\tlogExp_std = {}".format(logExp_mean, logExp_std))
	return SalmonCopyNumber, logExp_mean, logExp_std


def UpdateExpression(inexpfile, outexpfile, selection, TransGeneMap, logExp_mean):
	# read length and copy number of related transcripts
	copy_number = {}
	length = {}
	fpin = open(inexpfile, 'r')
	for line in fpin:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		tname = strs[0]
		gname = ""
		if tname in TransGeneMap:
			gname = TransGeneMap[tname]
		else:
			assert("ENSG" in tname)
			gname = tname.split("_")[0]
		if gname in selection:
			copy_number[tname] = float(strs[2])
			length[tname] = float(strs[1])
	fpin.close()
	# scale copy number and reads based on the mean
	for t,v in copy_number.items():
		if v == 0:
			copy_number[t] = 0.01
	cur_log_mean = np.mean( np.log( np.array(list(copy_number.values())) ) )
	multiplier = 1
	multiplier = max(1, logExp_mean / cur_log_mean)
	for t,v in copy_number.items():
		copy_number[t] *= multiplier
	numreads = {t:v*length[t] for t,v in copy_number.items()}
	# write new copy number ground truth
	fpout = open(outexpfile, 'w')
	fpout.write("# Name\tLength\tCopyNumber\tNumReads\n")
	for t,v in copy_number.items():
		fpout.write("{}\t{}\t{}\t{}\n".format(t, length[t], v, numreads[t]))
	fpout.close()
	print("Finish updating expression.")


if __name__=="__main__":
	# read GTF mainly to get TransGeneMap
	Transcripts = ReadGTF("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf")
	[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)
	# read graphs and eq_classes
	graphs, eq_classes = load_data("/home/hongyuz/src/gsalmon_fork/sim2/test")
	# read simulated TransNodes
	TransNodes = ReadSimulatedTransNodes("/home/congm1/savanna/savannacong33/GraphSalmonSimu/simu_Full_GEU_new_sequence_header.fa", TransGeneMap)
	# make selection
	genes_multi = Genes_with_multialign(eq_classes)
	selection = Select_uniqalign_genes(graphs, genes_multi, TransNodes)
	# estimate the mean and std of log expression
	_,logExp_mean, logExp_std = EstimatingExpDistribution("/home/congm1/savanna/savannacong33/SADrealdata/HumanBodyMap/salmon_Full_ERR030872/quant.sf")
	# update GTF, truth fasta, reference fasta, and truth expression file
	UpdateGTFfile("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf", "/mnt/disk33/user/congm1/small_graphsalmon/simu01/ref.annotation.gtf", selection)
	UpdateFastafile("/home/congm1/savanna/savannacong33/GraphSalmonSimu/simu_Full_GEU_new_sequence.fa", "/mnt/disk33/user/congm1/small_graphsalmon/simu01/truth.sequence.fa", selection, TransGeneMap)
	UpdateFastafile("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.full.transcripts.fa", "/mnt/disk33/user/congm1/small_graphsalmon/simu01/ref.sequence.fa", selection, TransGeneMap)
	UpdateExpression("/home/congm1/savanna/savannacong33/GraphSalmonSimu/simu_Full_GEU_new_expression.txt", "/mnt/disk33/user/congm1/small_graphsalmon/simu01/truth.expression.txt", selection, TransGeneMap, logExp_mean)
