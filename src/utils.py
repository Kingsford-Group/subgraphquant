#!/bin/python

import sys
import numpy as np
import struct
import tqdm
import collections
import copy


Nucleotide={'A':'T', 'C':'G', 'G':'C', 'T':'A', 'R':'Y', 'Y':'R', 'S':'W', 'W':'S', 'K':'M', 'M':'K', 'B':'V', 'V':'B', 'D':'H', 'H':'D', 'N':'N', '.':'.', '-':'-'}

def ReverseComplement(seq):
	rcseq = [Nucleotide[o] for o in seq.upper()]
	rcseq = "".join(rcseq)
	rcseq = rcseq[::-1]
	return rcseq


def ReadTruth(truthfile, TransLength):
	TransFlow = {}
	fp = open(truthfile, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		if strs[0] in TransLength:
			TransFlow[strs[0]] = float(strs[1]) / TransLength[strs[0]]
	fp.close()
	return TransFlow


def ReadSalmon(salmonquant):
	TransFlow = {}
	fp = open(salmonquant, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		TransFlow[strs[0]] = float(strs[4]) / float(strs[2])
	fp.close()
	return TransFlow


def ReadTransEdges_fromgraph(Graphs):
	TransEdges = []
	for g in Graphs:
		tnames = sum([e.IncidentTranscripts for e in g.vEdges], [])
		tnames = list(set(tnames))
		for tname in tnames:
			edges = [e.ID for e in g.vEdges if tname in e.IncidentTranscripts]
			TransEdges.append( [g.GeneID, tname, edges] )
	return TransEdges


def ReadTransEdges(seqfile, Graphs, GraphNameIndex, TransGeneMap):
	TransEdges = []
	fp = open(seqfile, 'r')
	for line in fp:
		if line[0] != '>':
			continue
		strs = line.strip().split()
		# retrieve transcript name and gene name
		tname = strs[0][1:]
		gname = ""
		if tname in TransGeneMap:
			gname = TransGeneMap[tname]
		else:
			assert(tname[:4] == "ENSG")
			gname = "_".join(tname.split("_")[:-2])
		# retrieve node and path info
		nodes = []
		if not (gname in GraphNameIndex):
			continue
		g = Graphs[GraphNameIndex[gname]]
		for s in strs[1][1:-1].split(")("):
			info = s.split(",")
			nodes.append(int(info[0]))
			assert(g.vNodes[nodes[-1]].Chr == info[1] and g.vNodes[nodes[-1]].StartPos == int(info[2]) and g.vNodes[nodes[-1]].EndPos == int(info[3]))
		edges = []
		for i in range(len(nodes)-1):
			idx_v1 = nodes[i]
			idx_v2 = nodes[i+1]
			this_edge = [idx_e for idx_e in g.vNodes[idx_v1].OutEdges if g.vEdges[idx_e].Node_2 == idx_v2]
			assert(len(this_edge) == 1)
			edges.append(this_edge[0])
		TransEdges.append( (gname,tname,edges) )
	fp.close()
	return TransEdges


def ConvertEdges2Nodes(TransEdges, Graphs, GraphNameIndex):
	TransNodes = []
	for info in TransEdges:
		(gname,tname,edges) = info
		assert(gname in GraphNameIndex)
		g = Graphs[GraphNameIndex[gname]]
		nodes = [g.vEdges[idx_e].Node_1 for idx_e in edges] + [g.vEdges[edges[-1]].Node_2]
		assert(nodes[0] == 0 and nodes[-1] == len(g.vNodes)-1)
		TransNodes.append( (gname,tname,nodes) )
	return TransNodes


def ReadGenome(fafile):
	Genome = {}
	fp = open(fafile, 'r')
	chrname = ""
	chrseq = ""
	for line in fp:
		if line[0] == '>':
			if chrname != "":
				Genome[chrname] = chrseq
			chrname = line.strip().split()[0][1:]
			chrseq = ""
		else:
			chrseq += line.strip()
	fp.close()
	if chrname != "":
		Genome[chrname] = chrseq
	return Genome


def WriteTopExpressedPaths(outputfile, Graphs, Genome):
	fp = open(outputfile, 'w')
	for g in Graphs:
		total_flow = np.sum([e.Flow for e in g.vEdges if e.Node_1 == 0])
		# find top k paths that can explain at least 90% of total flow
		k = 2
		paths, flows = TopKPaths(g, k)
		while np.sum(flows) < 0.9 * total_flow:
			k += 1
			paths, flows = TopKPaths(g, k)
		print([g.GeneID, k])
		# find corresponding nodes and sequence of each path
		for i in range(len(paths)):
			p = paths[i]
			f = flows[i]
			nodes = [g.vEdges[i].Node_1 for i in p] + [g.vEdges[p[-1]].Node_2]
			assert(nodes[0] == 0 and nodes[-1] == len(g.vNodes)-1)
			nodes = nodes[1:-1]
			seq = ""
			for idx_v in nodes:
				if g.vNodes[idx_v].Strand:
					seq += Genome[g.vNodes[idx_v].Chr][g.vNodes[idx_v].StartPos:g.vNodes[idx_v].EndPos]
				else:
					seq += ReverseComplement( Genome[g.vNodes[idx_v].Chr][g.vNodes[idx_v].StartPos:g.vNodes[idx_v].EndPos] )
			assert(len(seq) == np.sum([g.vNodes[idx_v].Length for idx_v in nodes]))
			# write sequence to file
			fp.write(">" + g.GeneID +"_"+ str(i) +" "+ str(f) +" "+ ",".join([str(x) for x in nodes]) +"\n")
			count = 0
			while count < len(seq):
				fp.write(seq[count:min(count+70, len(seq))] + "\n")
				count += 70
	fp.close()


def ReadCorrection(correct_file):
	corrections = {}
	fp = open(correct_file, 'rb')
	numtrans = struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen = struct.unpack('i', fp.read(4))[0]
		seqlen = struct.unpack('i', fp.read(4))[0]
		name = ""
		bias = np.zeros(seqlen)
		for j in range(namelen):
			name += struct.unpack('c', fp.read(1))[0].decode('utf-8')
		for j in range(seqlen):
			bias[j] = struct.unpack('d', fp.read(8))[0]
		corrections[name] = bias
	fp.close()
	return corrections


def ReadCorrection_matrix(correct_file):
	corrections = {}
	fp = open(correct_file, 'rb')
	numtrans = struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen = struct.unpack('i', fp.read(4))[0]
		seqlen = struct.unpack('i', fp.read(4))[0]
		name = ""
		bias = np.zeros(seqlen)
		for j in range(namelen):
			name += struct.unpack('c', fp.read(1))[0].decode('utf-8')
		for j in range(seqlen):
			bias[j] = struct.unpack('d', fp.read(8))[0]
		corrections[name] = bias
	fp.close()
	for k,v in corrections.items():
		l = int(np.sqrt(len(v)))
		corrections[k] = v.reshape( (l,l) )
	return corrections


def ReadCorrection_sparsematrix(correct_file):
	corrections = {}
	fp = open(correct_file, 'rb')
	numtrans = struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen = struct.unpack('i', fp.read(4))[0]
		seqlen = struct.unpack('i', fp.read(4))[0]
		name = ""
		bias = np.zeros((seqlen, 3))
		for j in range(namelen):
			name += struct.unpack('c', fp.read(1))[0].decode('utf-8')
		for j in range(seqlen):
			bias[j][0] = struct.unpack('i', fp.read(4))[0]
			bias[j][1] = struct.unpack('i', fp.read(4))[0]
			bias[j][2] = struct.unpack('d', fp.read(8))[0]
		corrections[name] = bias
	fp.close()
	return corrections


def UpdateNodeBiasMultiplier(sequence_file, corrections, Graphs, GraphNameIndex, damping = 0.1):
	# read paths and flows from sequence_file
	trans_flow = {}
	trans_nodes = {}
	fp = open(sequence_file, 'r')
	for line in fp:
		if line[0] == '>':
			strs = line.strip().split(" ")
			trans_flow[strs[0][1:]] = float(strs[1])
			trans_nodes[strs[0][1:]] = [int(x) for x in strs[2].split(",")]
	fp.close()
	# calculate the new node bias multiplier from each path
	# map from GeneID to a list of (node ID, flow, multiplier)
	NewMultiplier = {}
	for tname, bias in corrections.items():
		gname = tname.split("_")[0]
		g = Graphs[GraphNameIndex[gname]]
		f = trans_flow[tname]
		nodes = trans_nodes[tname]
		covered_length = 0
		for idx_v in nodes:
			tmp = (idx_v, f, np.sum(bias[covered_length:(covered_length + g.vNodes[idx_v].EndPos - g.vNodes[idx_v].StartPos)]) / g.vNodes[idx_v].Length)
			covered_length += g.vNodes[idx_v].Length
			assert(covered_length <= len(bias))
			if gname in NewMultiplier:
				NewMultiplier[gname].append( tmp )
			else:
				NewMultiplier[gname] = [tmp]
	# for each graph, calculate new bias multiplier, new multiplier = (1 - damping) * weighted sum of new + damping * old
	for gname, newlist in NewMultiplier.items():
		newlist.sort(key = lambda x:x[0])
		# get all bias correction for one node
		idx_s = 0
		idx_t = 1
		while idx_s < len(newlist):
			while idx_t < len(newlist) and newlist[idx_t][0] == newlist[idx_s][0]:
				idx_t += 1
			sum_new_flow = np.sum([x[1] for x in newlist[idx_s:idx_t]])
			new_bias = np.sum([x[2] * x[1] / sum_new_flow for x in newlist[idx_s:idx_t]])
			old_bias = Graphs[GraphNameIndex[gname]].vNodes[newlist[idx_s][0]].BiasMultiplier
			Graphs[GraphNameIndex[gname]].vNodes[newlist[idx_s][0]].BiasMultiplier = (1 - damping) * new_bias + damping * old_bias
			idx_s = idx_t
	return Graphs


def ExpandReference(Graphs, max_num_paths_pergene = 500):
	# enumerate all s-t paths for each gene, use as salmon reference
	# the following variable records the s-t paths using edge index, in a tuple (gname, tname, path in edge index)
	TransPaths = []
	with tqdm.tqdm(total=len(Graphs)) as pbar:
		for g in Graphs:
			gname = g.GeneID
			paths = g.enumerate_paths(max_num_paths_pergene)
			paths = [tuple(p) for p  in paths]
			# separate the paths into known transcripts and unknown transcripts
			trans_existing = sum([e.IncidentTranscripts for e in g.vEdges], [])
			trans_existing = list(set(trans_existing))
			path_existing = [tuple([e.ID for e in g.vEdges if tname in e.IncidentTranscripts]) for tname in trans_existing]
			path_new = []
			for p in paths:
				if not (p in path_existing):
					path_new.append(p)
			# record the paths
			for i in range(len(path_existing)):
				TransPaths.append( (gname, trans_existing[i], path_existing[i]) )
			for i in range(len(path_new)):
				TransPaths.append( (gname, gname+"_ext_"+format(i, '04d'), path_new[i]) )
			pbar.update(1)
	print("There are {} paths in total after expansion.".format(len(TransPaths)))
	return TransPaths


def TransPath2Sequence(outputfile, TransPaths, Graphs, GraphNameIndex, Genome):
	fp = open(outputfile, 'w')
	for i in range(len(TransPaths)):
		(gname, tname, path) = TransPaths[i]
		g = Graphs[GraphNameIndex[gname]]
		nodes = [g.vEdges[e].Node_1 for e in path] + [g.vEdges[path[-1]].Node_2]
		node_str = ""
		seq = ""
		for idx_v in nodes:
			node_str += "(" + str(idx_v) +","+ str(g.vNodes[idx_v].Chr) +","+ str(g.vNodes[idx_v].StartPos) +","+ str(g.vNodes[idx_v].EndPos) +","+ str(g.vNodes[idx_v].Strand) +")"
			if g.vNodes[idx_v].Strand:
				seq += Genome[g.vNodes[idx_v].Chr][g.vNodes[idx_v].StartPos:g.vNodes[idx_v].EndPos]
			else:
				seq += ReverseComplement( Genome[g.vNodes[idx_v].Chr][g.vNodes[idx_v].StartPos:g.vNodes[idx_v].EndPos] )
		assert(len(seq) == np.sum([g.vNodes[idx_v].Length for idx_v in nodes]) - 2)
		# write transcript name and sequence
		fp.write(">" + tname + "\n")
		count = 0
		while count < len(seq):
			fp.write(seq[count:min(count+70, len(seq))] + "\n")
			count += 70
	fp.close()


def TransPath2Header(outputfile, TransPaths, Graphs, GraphNameIndex):
	fp = open(outputfile, 'w')
	for i in range(len(TransPaths)):
		(gname, tname, path) = TransPaths[i]
		g = Graphs[GraphNameIndex[gname]]
		nodes = [g.vEdges[e].Node_1 for e in path] + [g.vEdges[path[-1]].Node_2]
		node_str = ""
		for idx_v in nodes:
			node_str += "(" + str(idx_v) +","+ str(g.vNodes[idx_v].Chr) +","+ str(g.vNodes[idx_v].StartPos) +","+ str(g.vNodes[idx_v].EndPos) +","+ str(g.vNodes[idx_v].Strand) +")"
		# write transcript name and sequence
		fp.write(">" + tname + " " + node_str + "\n")
	fp.close()


def TransPath_summary(outputfile, TransPaths, Graphs):
	gene_path_map = {}
	for i in range(len(TransPaths)):
		gname = TransPaths[i][0]
		if gname in gene_path_map:
			gene_path_map[gname].append(i)
		else:
			gene_path_map[gname] = [i]
	# summarize for each gene, the number of nodes, edges, existing transcripts, novel transcripts
	fp = open(outputfile, 'w')
	fp.write("# GeneID\tNumNodes\tNumEdges\tNumExisting\tNumNovel\n")
	for g in Graphs:
		gname = g.GeneID
		assert(gname in gene_path_map)
		this_trans_paths = [TransPaths[i] for i in gene_path_map[gname]]
		this_trans_names = [info[1] for info in this_trans_paths]
		num_existing = len([tname for tname in this_trans_names if tname[:4] == "ENST"])
		num_new = len([tname for tname in this_trans_names if tname[:4] == "ENSG"])
		assert(num_existing + num_new == len(this_trans_names))
		fp.write("{}\t{}\t{}\t{}\t{}\n".format(gname, len(g.vNodes), len(g.vEdges), num_existing, num_new))
	fp.close()


def ReadTransNodes(seqfile, TransGeneMap = None):
	TransNodes = []
	fp = open(seqfile, 'r')
	for line in fp:
		if line[0] != '>':
			continue
		strs = line.strip().split()
		# retrieve transcript name and gene name
		tname = strs[0][1:]
		gname = ""
		if not (TransGeneMap is None) and tname in TransGeneMap:
			gname = TransGeneMap[tname]
		else:
			assert(tname[:4] == "ENSG")
			gname = "_".join(tname.split("_")[:-2])
		# retrieve node and path info
		nodes = []
		for s in strs[-1][1:-1].split(")("):
			info = s.split(",")
			nodes.append(int(info[0]))
		TransNodes.append( (gname,tname,nodes) )
	fp.close()
	return TransNodes


def ReadTransLogWeights(seqfile):
	log_weights = {}
	fp = open(seqfile, 'r')
	for line in fp:
		if line[0] == '>':
			strs = line.strip().split()
			tname = strs[0][1:]
			log_w = float(strs[1])
			log_weights[tname] = log_w
	fp.close()
	return log_weights


Interval = collections.namedtuple("Interval", ['chr','startpos','endpos','gid'])
def Graph2Interval(graphs):
	intervals = []
	for g in graphs:
		gid = g.GeneID
		chr = g.vNodes[1].Chr
		s = np.min([v.StartPos for v in g.vNodes[1:-1]])
		t = np.max([v.EndPos for v in g.vNodes[1:-1]])
		intervals.append( Interval(chr=chr, startpos=s, endpos=t, gid=gid) )
	intervals.sort(key = lambda x:(x.chr, x.startpos))
	return intervals


def BinarySearchPoint(intervals, chr, pos):
	indexes = []
	gids = []
	lo = 0
	hi = len(intervals)
	while lo < hi:
		mid = int((lo + hi) / 2)
		if intervals[mid].chr < chr or (intervals[mid].chr == chr and intervals[mid].endpos <= pos):
			lo = mid + 1
		elif intervals[mid].chr == chr and intervals[mid].startpos <= pos and intervals[mid].endpos > pos:
			lo = mid
			hi = mid
			break
		elif intervals[mid].chr > chr or (intervals[mid].chr == chr and intervals[mid].startpos > pos):
			hi = mid - 1
		else:
			print("error: chr = {}, pos = {}, lo = {}, mid = {}, hi = {}".format(chr, pos, intervals[lo], intervals[mid], intervals[hi]))
	# because of interval searching, search nearby regions around lo and hi
	assert(lo >= hi)
	i = lo
	while i >= 0 and i >= lo - 50:
		if intervals[i].chr == chr and intervals[i].startpos <= pos and intervals[i].endpos > pos:
			indexes.append(i)
			gids.append(intervals[i].gid)
		i -= 1
	i = lo + 1
	while i <= lo + 50 and i < len(intervals):
		if intervals[i].chr == chr and intervals[i].startpos <= pos and intervals[i].endpos > pos:
			indexes.append(i)
			gids.append(intervals[i].gid)
		i += 1
	assert( len(indexes) == len(gids) )
	return indexes, gids


def BinarySearchRegion(intervals, chr, region):
	indexes = []
	gids = []
	lo = 0
	hi = len(intervals)
	while lo < hi:
		mid = int((lo + hi) / 2)
		if intervals[mid].chr < chr or (intervals[mid].chr == chr and intervals[mid].endpos <= region[0]):
			lo = mid + 1
		elif intervals[mid].chr == chr and intervals[mid].startpos <= region[1] and intervals[mid].endpos > region[0]:
			lo = mid
			hi = mid
			break
		elif intervals[mid].chr > chr or (intervals[mid].startpos > region[1]):
			hi = mid - 1
		else:
			print("error: chr = {}, region = {}, lo = {}, mid = {}, hi = {}".format(chr, region, intervals[lo], intervals[mid], intervals[hi]))
	# because of interval searching, search nearby regions around lo and hi
	assert(lo >= hi)
	i = lo
	while i >= 0 and i >= lo - 50:
		if intervals[i].chr == chr and intervals[i].startpos <= region[1] and intervals[i].endpos > region[0]:
			indexes.append(i)
			gids.append(intervals[i].gid)
		i -= 1
	i = lo + 1
	while i <= lo + 50 and i < len(intervals):
		if intervals[i].chr == chr and intervals[i].startpos <= region[1] and intervals[i].endpos > region[0]:
			indexes.append(i)
			gids.append(intervals[i].gid)
		i += 1
	assert( len(indexes) == len(gids) )
	return indexes, gids


def FindPaths4Transcript(graphs, NameIndex, intervals, t):
	TransPaths = []
	indexes, gids = BinarySearchRegion(intervals, t.Chr, [t.StartPos, t.EndPos])
	for i in range(len(indexes)):
		g = graphs[NameIndex[gids[i]]]
		exons = copy.deepcopy(t.Exons)
		if not g.vNodes[1].Strand:
			exons = exons[::-1]
		# loop over all exons and all graph nodes to find index of nodes that overlap with the exons of transcript
		overlap_nodes = [0]
		idx_exon = 0
		idx_node = 1
		while idx_exon < len(exons) and idx_node + 1 < len(g.vNodes):
			if g.vNodes[1].Strand:
				# exon to the left of nodes
				if exons[idx_exon][1] <= g.vNodes[idx_node].StartPos:
					idx_exon += 1
				# overlap
				elif max(exons[idx_exon][0], g.vNodes[idx_node].StartPos) < min(exons[idx_exon][1], g.vNodes[idx_node].EndPos):
					overlap_nodes.append(idx_node)
					idx_node += 1
				# node to the left of exon
				else:
					idx_node += 1
			else:
				# exon to the right of nodes
				if exons[idx_exon][0] >= g.vNodes[idx_node].EndPos:
					idx_exon += 1
				# overlap
				elif max(exons[idx_exon][0], g.vNodes[idx_node].StartPos) < min(exons[idx_exon][1], g.vNodes[idx_node].EndPos):
					overlap_nodes.append(idx_node)
					idx_node += 1
				# node to the right of exon
				else:
					idx_node += 1
		overlap_nodes.append(g.vNodes[-1].ID)
		overlap_nodes = list(set(overlap_nodes))
		overlap_nodes.sort()
		# convert node list to edge list
		edges = []
		all_edges = {(e.Node_1, e.Node_2):e.ID for e in g.vEdges}
		for j in range(1, len(overlap_nodes)):
			if (overlap_nodes[j-1], overlap_nodes[j]) in all_edges:
				edges.append( all_edges[ (overlap_nodes[j-1], overlap_nodes[j]) ] )
		if len(edges) != 0:
			TransPaths.append( (g.GeneID, t.TransID, edges) )
	return TransPaths


def PathsFromAssembly(asm_transcripts, graphs):
	intervals = Graph2Interval(graphs)
	NameIndex = {graphs[i].GeneID:i for i in range(len(graphs))}
	TransPaths = []
	for tname, t in asm_transcripts.items():
		paths = FindPaths4Transcript(graphs, NameIndex, intervals, t)
		TransPaths += paths
	return TransPaths
