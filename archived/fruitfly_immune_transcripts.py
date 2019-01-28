#!/bin/python

import sys
import numpy as np
import subprocess
from GeneGraph import *
from utils import *
from TranscriptClass import *


def UpdateGTF(ingtf, outgtf, involved_genes, graphs, GraphNameIndex, transpaths):
	fpin = open(ingtf, 'r')
	fpout = open(outgtf, 'w')
	havewritten = {g:False for g in involved_genes}
	count = 0
	for line in fpin:
		if line[0] == '#':
			continue
		gid = GetFeature(line, "gene_id")
		if not (gid in involved_genes):
			fpout.write(line)
		elif not havewritten[gid]:
			# write the gene and expanded transcripts
			strs = line.strip().split("\t")
			assert(strs[2] == "gene")
			fpout.write(line)
			# select relevant transcripts to write
			for info in transpaths:
				(gname,tname,edges) = info
				if gname == gid:
					g = graphs[GraphNameIndex[gid]]
					nodes = [g.vEdges[idx_e].Node_1 for idx_e in edges]
					assert(nodes[0] == 0)
					nodes = nodes[1:]
					chr = g.vNodes[nodes[0]].Chr
					startpos = np.min([g.vNodes[idx_v].StartPos for idx_v in nodes])
					endpos = np.max([g.vNodes[idx_v].EndPos for idx_v in nodes])
					strand = "-"
					if g.vNodes[nodes[0]].Strand:
						strand = "+"
					fpout.write("{}\tExpand\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\";\n".format(chr, startpos, endpos, strand, gname, tname))
					for idx_v in nodes:
						fpout.write("{}\tExpand\texon\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\";\n".format(chr, g.vNodes[idx_v].StartPos, g.vNodes[idx_v].EndPos, strand, gname, tname))
					count += 1
			havewritten[gid] = True
	fpin.close()
	fpout.close()
	print("Finish writing annotation GTF. Additional {} are written.".format(count))


def UpdateSequenceHeader(outheader, involved_genes, graphs, GraphNameIndex, transpaths):
	fpout = open(outheader, 'w')
	for g in graphs:
		if not (g.GeneID in involved_genes):
			trans_existing = sum([e.IncidentTranscripts for e in g.vEdges], [])
			trans_existing = list(set(trans_existing))
			node_existing = [tuple([e.Node_1 for e in g.vEdges if tname in e.IncidentTranscripts] + [g.vNodes[-1].ID]) for tname in trans_existing]
			for i in range(len(trans_existing)):
				tname = trans_existing[i]
				nodes = node_existing[i]
				node_str = ""
				for idx_v in nodes:
					node_str += "(" + str(idx_v) +","+ str(g.vNodes[idx_v].Chr) +","+ str(g.vNodes[idx_v].StartPos) +","+ str(g.vNodes[idx_v].EndPos) +","+ str(g.vNodes[idx_v].Strand) +")"
				fpout.write(">" + tname + " " + node_str + "\n")
		else:
			for info in transpaths:
				(gname, tname, edges) = info
				if gname == g.GeneID:
					nodes = tuple([g.vEdges[idx_e].Node_1 for idx_e in edges] + [g.vNodes[-1].ID])
					node_str = ""
					for idx_v in nodes:
						node_str += "(" + str(idx_v) +","+ str(g.vNodes[idx_v].Chr) +","+ str(g.vNodes[idx_v].StartPos) +","+ str(g.vNodes[idx_v].EndPos) +","+ str(g.vNodes[idx_v].Strand) +")"
					fpout.write(">" + tname + " " + node_str + "\n")
	fpout.close()
	print("Finish writing header.")


if __name__=="__main__":
	graphs = ReadGraphFile("/home/congm1/savanna/savannacong33/fruitfly/Drosophila_melanogaster.BDGP6.95.chr.graph_fragstart_beforebias.txt")
	GraphNameIndex = {graphs[i].GeneID:i for i in range(len(graphs))}
	involved_genes = set(["FBgn0033159"])
	transpaths = ExpandReference( [graphs[GraphNameIndex[gname]] for gname in involved_genes], max_num_paths_pergene=1000000)
	UpdateGTF("/home/congm1/savanna/savannacong33/fruitfly/Drosophila_melanogaster.BDGP6.95.chr.gtf", \
		"/home/congm1/savanna/savannacong33/fruitfly/2Drosophila_melanogaster.immune.expand.gtf", \
		involved_genes, graphs, GraphNameIndex, transpaths)
	UpdateSequenceHeader("/home/congm1/savanna/savannacong33/fruitfly/Drosophila_melanogaster.immune.expand.transcript.header.fa",\
		involved_genes, graphs, GraphNameIndex, transpaths)
	# update fasta
	cmd = "/home/congm1/ocean/oceancong02/Code/SalmonInput/GtftoFa "
	cmd += "/home/congm1/savanna/savannacong33/fruitfly/Drosophila_melanogaster.immune.expand.gtf "
	cmd += "/home/congm1/savanna/savannacong33/fruitfly/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa "
	cmd += "/home/congm1/savanna/savannacong33/fruitfly/Drosophila_melanogaster.immune.expand.transcript.fa"
	p = subprocess.Popen(cmd, shell=True)
	out,err = p.communicate()
