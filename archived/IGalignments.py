#!/bin/python

import sys
import numpy as np
import pysam
import collections
from TranscriptClass import *


Interval = collections.namedtuple("Interval", ['chr','startpos','endpos','gid'])

def ReadIGPositions(gtffile, gene_type = "IG_V_gene"):
	intervals = []
	fp = open(gtffile, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[2] == "exon":
			if GetFeature(line, "gene_type") == gene_type:
				gid= GetFeature(line, "gene_id")
				intervals.append( Interval(chr=strs[0], startpos=int(strs[3])-1, endpos=int(strs[4]), gid=gid) )
	fp.close()
	intervals.sort(key = lambda x:(x.chr, x.startpos))
	return intervals


def FilterIntervalChr(intervals, chrname):
	intervals = [x for x in intervals if x.chr == chrname]
	return intervals


def WithinRange(genes_range, read):
	if read.reference_name != genes_range.chr:
		return False
	elif read.reference_start > genes_range.endpos:
		return False
	elif read.reference_end < genes_range.startpos:
		return False
	return True


def ChechAlignment(bamfile, IGVintervals, IGDintervals, IGJintervals):
	assert(len(set([x.chr for x in IGVintervals])) == 1)
	assert(len(set([x.chr for x in IGDintervals])) == 1)
	assert(len(set([x.chr for x in IGJintervals])) == 1)
	IGV_range = Interval(chr=IGVintervals[0].chr, startpos=np.min([x.startpos for x in IGVintervals]), endpos=np.max([x.endpos for x in IGVintervals]), gid="IGV")
	IGD_range = Interval(chr=IGDintervals[0].chr, startpos=np.min([x.startpos for x in IGDintervals]), endpos=np.max([x.endpos for x in IGDintervals]), gid="IGD")
	IGJ_range = Interval(chr=IGJintervals[0].chr, startpos=np.min([x.startpos for x in IGJintervals]), endpos=np.max([x.endpos for x in IGJintervals]), gid="IGJ")
	fp = pysam.AlignmentFile(bamfile)
	potential_alignments = []
	for read in fp:
		within_V = WithinRange(IGV_range, read)
		within_D = WithinRange(IGD_range, read)
		within_J = WithinRange(IGJ_range, read)
		if not (within_V or within_D or within_J):
			continue
		else:
			potential_alignments.append(read)
	# retrieve the VDJ genes of each paired-end reads
	Overlapping_Reads = {}
	potential_alignments.sort(key = lambda x:x.query_name)
	align_blocks = []
	readname = potential_alignments[0].query_name
	for ali in potential_alignments:
		if ali.query_name != readname:
			# check which VDJ genes it overlaps
			V_genes = []
			D_genes = []
			J_genes = []
			# assertion for the align_blocks should be sorted by coordinate
			assert(len(align_blocks) > 0)
			align_blocks.sort(key = lambda x:x[0])
			# IGV genes
			idx_search = 0
			for b in align_blocks:
				while IGVintervals[idx_search].endpos < b[0]:
					idx_search += 1
				for itvl in IGVintervals[idx_search:]:
					if itvl.startpos > b[1]:
						break
					elif max(itvl.startpos, b[0]) < min(itvl.endpos, b[1]):
						V_genes.append(itvl.gid)
			# IGD genes
			idx_search = 0
			for b in align_blocks:
				while IGDintervals[idx_search].endpos < b[0]:
					idx_search += 1
				for itvl in IGDintervals[idx_search:]:
					if itvl.startpos > b[1]:
						break
					elif max(itvl.startpos, b[0]) < min(itvl.endpos, b[1]):
						D_genes.append(itvl.gid)
			# IGJ genes
			idx_search = 0
			for b in align_blocks:
				while IGJintervals[idx_search].endpos < b[0]:
					idx_search += 1
				for itvl in IGDintervals[idx_search:]:
					if itvl.startpos > b[1]:
						break
					elif max(itvl.startpos, b[0]) < min(itvl.endpos, b[1]):
						J_genes.append(itvl.gid)
			if len(V_genes) + len(D_genes) + len(J_genes) > 0:
				Overlapping_Reads[readname] = [V_genes, D_genes, J_genes]
			# add new read info
			readname = ali.query_name
			align_blocks = ali.get_blocks()
		else:
			align_blocks += ali.get_blocks()
	return potential_alignments, Overlapping_Reads


def WriteResult(templatebam, outprefix, potential_alignments, Overlapping_Reads):
	fpin = pysam.AlignmentFile(templatebam)
	fpout = pysam.AlignmentFile(outprefix+"_align.bam", 'wb', template=fpin)
	for read in potential_alignments:
		if read.query_name in Overlapping_Reads:
			fpout.write(read)
	fpin.close()
	fpout.close()
	fp = open(outprefix+"_gid.txt", 'w')
	fp.write("# query_name\tV_genes\tD_genes\tJ_genes\n")
	for readname,genes in Overlapping_Reads.items():
		V_genes, D_genes, J_genes = genes
		fp.write("{}\t{}\t{}\t{}\n".format(readname, V_genes, D_genes, J_genes))
	fp.close()


IGVintervals = ReadIGPositions("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf", "IG_V_gene")
IGVintervals = FilterIntervalChr(IGVintervals, "chr14")
IGDintervals = ReadIGPositions("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf", "IG_D_gene")
IGDintervals = FilterIntervalChr(IGDintervals, "chr14")
IGJintervals = ReadIGPositions("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf", "IG_J_gene")
IGJintervals = FilterIntervalChr(IGJintervals, "chr14")

# potential_alignments, Overlapping_Reads = ChechAlignment("/home/congm1/savanna/savannacong33/SADrealdata/HumanBodyMap/star_Full_ERR030878/Aligned.sortedByCoord.out.bam", IGVintervals, IGDintervals, IGJintervals)
# WriteResult("/home/congm1/savanna/savannacong33/SADrealdata/HumanBodyMap/star_Full_ERR030878/Aligned.sortedByCoord.out.bam", "/home/congm1/savanna/savannacong33/SADrealdata/HumanBodyMap/star_Full_ERR030878/VDJ", potential_alignments, Overlapping_Reads)