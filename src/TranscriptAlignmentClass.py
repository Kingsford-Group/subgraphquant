#!/bin/python

import sys
import copy
import numpy as np
import pysam
import tqdm
from pathlib import Path


def ReadSalmonTPM(quantfile):
	TPM = {}
	fp = open(quantfile, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		TPM[strs[0]] = float(strs[3])
	fp.close()
	return TPM


def GetSalmonAssignWeight(salmon_eq_file, TPM):
	fp = open(salmon_eq_file, 'r')
	n_trans = int(fp.readline().strip())
	n_eq = int(fp.readline().strip())
	# read the ordered transcript list
	trans_list = []
	for i in range(n_trans):
		trans_list.append( fp.readline().strip() )
	# read the aligned set of transcripts and part1 weight of each eq_classes
	# process the final weight by multiplying the TPM to the part1 weight and normalize to 1
	map_transtuple_weight = {}
	for i in range(n_eq):
		strs = fp.readline().strip().split("\t")
		m = int(strs[0])
		if len(strs) != 2 + 2 * m:
			print(strs)
		assert( len(strs) == 2 + 2 * m)
		trans = np.array([trans_list[int(x)] for x in strs[1:(1+m)]])
		part1weight = np.array([float(x) for x in strs[(1+m):-1]])
		assert( len(trans) == len(part1weight) )
		# sort the transtuple and corresponding part1 weight by the transname
		index = np.argsort(trans)
		trans = trans[index]
		part1weight = part1weight[index]
		# process the final weight
		assert( np.all([x in TPM for x in trans]) )
		weight = np.array([part1weight[j] * TPM[trans[j]] for j in range(m)])
		if np.sum(weight) == 0:
			weight = np.array([part1weight[j] for j in range(m)])
		assert( np.sum(weight) > 0 )
		weight /= np.sum(weight)
		map_transtuple_weight[tuple(list(trans))] = weight
	return map_transtuple_weight


def MergeBlocks(align_blocks):
	new_align_blocks = []
	align_blocks.sort(key = lambda x:x[0])
	for e in align_blocks:
		if len(new_align_blocks) == 0 or e[0] > new_align_blocks[-1][1]:
			new_align_blocks.append( [e[0], e[1]] )
		else:
			new_align_blocks[-1][1] = max(new_align_blocks[-1][1], e[1])
	new_align_blocks = [(x[0], x[1]) for x in new_align_blocks]
	return new_align_blocks


def ConvertCoord_trans2genome(t, blocks):
	assert(len(blocks) == 2 and blocks[0] <= blocks[1])
	newblocks = []
	covered_length = 0
	for e in t.Exons:
		if covered_length + e[1] - e[0] < blocks[0]:
			covered_length += e[1] - e[0]
		elif covered_length >= blocks[1]:
			break
		else:
			# the exon overlaps with the blocks
			overlap_s = max(blocks[0], covered_length)
			overlap_t = min(blocks[1], covered_length + e[1] - e[0])
			if t.Strand:
				newblocks.append( (e[0] + overlap_s - covered_length, e[0] + overlap_t - covered_length) )
			else:
				newblocks.append( (e[1] - (overlap_t - covered_length), e[1] - (overlap_s - covered_length)) )
			covered_length += e[1] - e[0]
	return newblocks


def ShortestPath(g, n1, n2):
	# starting from n1, breast first search
	visitied = [False] * len(g.vNodes)
	tobecheckes = [n1] # the middle elements that may lead to n2
	inedges = np.ones(len(g.vNodes), dtype=np.int) * (-1) # the edge that lead to the corresponding node with the smallest distance from n1
	distances = np.ones(len(g.vNodes)) * np.inf
	distances[n1] = 0
	while len(tobecheckes) > 0:
		last_node =  tobecheckes[0]
		tobecheckes = tobecheckes[1:]
		if visitied[last_node]:
			continue
		for e in [g.vEdges[i] for i in g.vNodes[last_node].OutEdges]:
			if e.Node_2 == n2 and distances[n2] > distances[last_node] + 1:
				inedges[n2] = e.ID
				distances[n2] = distances[last_node] + 1
				break
			elif distances[e.Node_2] > distances[last_node] + 1:
				inedges[e.Node_2] = e.ID
				distances[e.Node_2] = distances[last_node] + 1
				tobecheckes.append(e.Node_2)
		visitied[last_node] = True
	# now, the vertices are indicated with the edge with the shortest distance from n1
	# back trace from n2, follow the inedges, and construct the whole path from n1 to n2
	assert(inedges[n2] != -1)
	last_node = n2
	backtrace_edges = []
	while last_node != n1:
		backtrace_edges.append(inedges[last_node])
		last_node = g.vEdges[inedges[last_node]].Node_1
	return backtrace_edges[::-1]


class ReadAlignment_t(object):
	def __init__(self, alignments, Transcripts, transtuple = None, weight = None):
		self.ReadName = alignments[0].query_name
		self.GeneIDs = []
		self.Aligned_blocks = []
		self.Aligned_Nodes = []
		self.Aligned_Edges = []
		self.Aligned_Paths = []
		self.Weights = []
		self.IsMultiMapped = False
		# map transcripts to blocks
		TransIDs = []
		trans_blocks = {}
		for a in alignments:
			if a.reference_name in trans_blocks:
				trans_blocks[a.reference_name] += a.get_blocks()
			else:
				trans_blocks[a.reference_name] = a.get_blocks()
		gene_blocks = {}
		gene_weight = {}
		for tname,blocks in trans_blocks.items():
			g = Transcripts[tname].GeneID
			# merge the blocks into non-overlapping regions
			blocks = MergeBlocks(blocks)
			# convert to genomic coordinate
			blocks = sum([ConvertCoord_trans2genome(Transcripts[tname], b) for b in blocks], [])
			# print([tname, g, blocks])
			if g in gene_blocks:
				is_redundant = False
				for o in range(len(gene_blocks[g])):
					oldblocks = gene_blocks[g][o]
					is_equal = True
					if len(oldblocks) != len(blocks):
						is_equal = False
					else:
						for i in range(len(oldblocks)):
							if oldblocks[i] != blocks[i]:
								is_equal = False
					if is_equal:
						is_redundant = True
						if not (transtuple is None):
							assert(tname in transtuple)
							idx_tuple = transtuple.index(tname)
							gene_weight[g][o] += weight[idx_tuple]
						break
				if not is_redundant:
					gene_blocks[g].append(blocks)
					if not (transtuple is None):
						assert(tname in transtuple)
						idx_tuple = transtuple.index(tname)
						gene_weight[g].append( weight[idx_tuple] )
			else:
				gene_blocks[g] = [blocks]
				if not (transtuple is None):
					assert(tname in transtuple)
					idx_tuple = transtuple.index(tname)
					gene_weight[g] = [weight[idx_tuple]]
			assert( len(gene_weight) == len(gene_blocks) )
		list_gene_blocks = list(gene_blocks.items())
		list_gene_blocks.sort(key = lambda x:x[0])
		for (gname, multi_blocks) in list_gene_blocks:
			self.GeneIDs += [gname] * len(multi_blocks)
			self.Aligned_blocks += multi_blocks
			if not (transtuple is None):
				self.Weights += gene_weight[gname]
		if (transtuple is None):
			self.Weights = [1.0 / len(self.GeneIDs)] * len(self.GeneIDs)
		assert(len(self.GeneIDs) == len(self.Aligned_blocks))
		assert(len(self.GeneIDs) == len(self.Weights))
		assert( np.abs(np.sum(self.Weights) - 1) < 1e-4)
		# check whether is multi-mapped
		if len(self.GeneIDs) > 1:
			self.IsMultiMapped = True

	def find_overlapping_nodes(self, Graphs, GraphNameIndex):
		shrink_GeneIDs = []
		shrink_Aligned_blocks = []
		shrink_Aligned_Nodes = []
		shrink_Aligned_Edges = []
		for i in range(len(self.GeneIDs)):
			gname = self.GeneIDs[i]
			blocks = self.Aligned_blocks[i]
			assert(gname in GraphNameIndex)
			g = Graphs[GraphNameIndex[gname]]
			overlapping_nodes = []
			for v in g.vNodes:
				for b in blocks:
					if max(v.StartPos, b[0]) < min(v.EndPos, b[1]):
						overlapping_nodes.append(v.ID)
			overlapping_nodes = list(set(overlapping_nodes))
			overlapping_nodes.sort()
			set_overlapping_nodes = set(overlapping_nodes)
			# find the edges between aligned nodes
			# overlapping_edges = []
			# for j in range(1, len(overlapping_nodes)):
			# 	overlapping_edges += ShortestPath(g, overlapping_nodes[j-1], overlapping_nodes[j])
			# find overlapping edges by check whether each transcript contain all the nodes
			trans_nodes_map = {}
			for e in g.vEdges:
				for t in e.IncidentTranscripts:
					if t in trans_nodes_map:
						trans_nodes_map[t].append( e.Node_1 )
					else:
						trans_nodes_map[t] = [e.Node_1]
			overlapping_edges = []
			for t,nl in trans_nodes_map.items():
				nl.append( g.vNodes[-1].ID )
				setnl = set(nl)
				if len(setnl & set_overlapping_nodes) == len(overlapping_nodes):
					tmpoverlapping_edges = [e.ID for e in g.vEdges if t in e.IncidentTranscripts and e.Node_1 >= overlapping_nodes[0] and e.Node_2 <= overlapping_nodes[-1]]
					if len(overlapping_edges) == 0 or len(tmpoverlapping_edges) < len(overlapping_edges):
						overlapping_edges = tmpoverlapping_edges
			overlapping_edges = list(set(overlapping_edges))
			overlapping_edges.sort()
			# sanity check: make sure that the edges are connected
			for j in range(1, len(overlapping_edges)):
				assert(g.vEdges[overlapping_edges[j-1]].Node_2 == g.vEdges[overlapping_edges[j]].Node_1)
			# add Aligned_Edges and Aligned_Nodes to attributes
			self.Aligned_Edges.append(overlapping_edges)
			if len(overlapping_edges) > 0:
				self.Aligned_Nodes.append([g.vEdges[e].Node_1 for e in overlapping_edges] + [g.vEdges[overlapping_edges[-1]].Node_2])
			else:
				self.Aligned_Nodes.append(overlapping_nodes)
			assert( (len(self.Aligned_Nodes[-1]) == 0 and len(self.Aligned_Edges[-1]) == 0) or (len(self.Aligned_Nodes[-1]) == len(self.Aligned_Edges[-1]) + 1))
			# handle the case where the corresponding node is deleted by the 0 bias-corrected expected coverage
			if len(self.Aligned_Nodes[0]) != 0:
				shrink_GeneIDs.append( self.GeneIDs[i] )
				shrink_Aligned_blocks.append( self.Aligned_blocks[i] )
				shrink_Aligned_Nodes.append( self.Aligned_Nodes[i] )
				shrink_Aligned_Edges.append( self.Aligned_Edges[i] )
				assert(len(overlapping_nodes) == 1 or len(overlapping_edges) > 0)
		self.GeneIDs = shrink_GeneIDs
		self.Aligned_blocks = shrink_Aligned_blocks
		self.Aligned_Nodes = shrink_Aligned_Nodes
		self.Aligned_Edges = shrink_Aligned_Edges
		assert(len(self.GeneIDs) == len(self.Aligned_Nodes))
		assert(len(self.GeneIDs) == len(self.Aligned_Edges))
		# make nodes and edges interleaving to generate the path
		for i in range(len(self.GeneIDs)):
			path = sum([[self.Aligned_Nodes[i][j], self.Aligned_Edges[i][j]] for j in range(len(self.Aligned_Edges[i]))], []) + [self.Aligned_Nodes[i][-1]]
			self.Aligned_Paths.append(path)

	def remove_duplicated_alignment(self):
		# if a read is aligned to the same gene, with same set of paths
		unique_align_index = []
		new_GeneIDs = []
		new_Aligned_blocks = []
		new_Aligned_Nodes = []
		new_Aligned_Edges = []
		new_Aligned_Paths = []
		new_Weights = []
		for i in range(len(self.GeneIDs)):
			is_duplicates = False
			for j in range(len(new_GeneIDs)):
				if self.GeneIDs[i] == new_GeneIDs[j] and self.Aligned_Paths[i] == new_Aligned_Paths[j]:
					is_duplicates = True
					new_Weights[j] += self.Weights[i]
					break
			if not is_duplicates:
				new_GeneIDs.append( self.GeneIDs[i] )
				new_Aligned_blocks.append( self.Aligned_blocks[i] )
				new_Aligned_Nodes.append( self.Aligned_Nodes[i] )
				new_Aligned_Edges.append( self.Aligned_Edges[i] )
				new_Aligned_Paths.append( self.Aligned_Paths[i] )
				new_Weights.append( self.Weights[i] )
		self.GeneIDs = new_GeneIDs
		self.Aligned_blocks = new_Aligned_blocks
		self.Aligned_Nodes = new_Aligned_Nodes
		self.Aligned_Edges = new_Aligned_Edges
		self.Aligned_Paths = new_Aligned_Paths
		self.Weights = new_Weights

	def sort_alignments(self):
		index = list(np.arange(len(self.GeneIDs)))
		index.sort(key = lambda o : ([self.GeneIDs[o]] + self.Aligned_Paths[o]))
		self.GeneIDs = [self.GeneIDs[o] for o in index]
		self.Aligned_blocks = [self.Aligned_blocks[o] for o in index]
		self.Aligned_Nodes = [self.Aligned_Nodes[o] for o in index]
		self.Aligned_Edges = [self.Aligned_Edges[o] for o in index]
		self.Aligned_Paths = [self.Aligned_Paths[o] for o in index]

	def convert2tuple(self):
		self.GeneIDs = tuple(self.GeneIDs)
		self.Aligned_Nodes = [tuple(x) for x in self.Aligned_Nodes]
		self.Aligned_Edges = [tuple(x) for x in self.Aligned_Edges]
		self.Aligned_Paths = [tuple(x) for x in self.Aligned_Paths]

	def convert2list(self):
		self.GeneIDs = list(self.GeneIDs)
		self.Aligned_Nodes = [list(x) for x in self.Aligned_Nodes]
		self.Aligned_Edges = [list(x) for x in self.Aligned_Edges]
		self.Aligned_Paths = [list(x) for x in self.Aligned_Paths]


class EquivalentClass_t(object):
	def __init__(self, ra):
		if not (ra is None):
			self.Count = 1
			self.GeneIDs = copy.copy(ra.GeneIDs)
			self.Aligned_Nodes = copy.copy(ra.Aligned_Nodes)
			self.Aligned_Edges = copy.copy(ra.Aligned_Edges)
			self.Aligned_Paths = copy.copy(ra.Aligned_Paths)
			self.Weights = np.array(ra.Weights)
		else:
			self.Count = 0
			self.GeneIDs = []
			self.Aligned_Nodes = []
			self.Aligned_Edges = []
			self.Aligned_Paths = []
			self.Weights = np.array([])

	def IsEqual(self, ra):
		if self.GeneIDs != ra.GeneIDs:
			return False
		for i in range(len(self.Aligned_Paths)):
			# compare paths
			if self.Aligned_Paths[i] != ra.Aligned_Paths[i]:
				return False
		return True

	def addread(self, ra):
		self.Count += 1
		self.Weights +=  np.array(ra.Weights)

	def initializeweight(self):
		self.Weights = 1.0 * self.Count * np.ones(len(self.GeneIDs)) / len(self.GeneIDs)

	def assignweights(self, Graphs, GraphNameIndex, Opts):
		flows = []
		for i in range(len(self.GeneIDs)):
			gname = self.GeneIDs[i]
			path = self.Aligned_Paths[i]
			g = Graphs[GraphNameIndex[gname]]
			opt = Opts[GraphNameIndex[gname]]
			# retrieve optimized flow of the graph
			opt_results = opt.Results
			if opt_results is None:
				opt_results = opt.initialflow()
				assert(np.sum([opt_results[e.ID] * g.vNodes[e.Node_1].Length * g.vNodes[e.Node_1].BiasMultiplier for e in g.vEdges]) > 0)
				s = np.sum(opt.Weights)  / np.sum([opt_results[e.ID] * g.vNodes[e.Node_1].Length * g.vNodes[e.Node_1].BiasMultiplier for e in g.vEdges])
				opt_results *= s
			# find the corresponding edge / hyper edge
			assert(g.GeneID == opt.graph.GeneID and len(g.vNodes) == len(opt.graph.vNodes) and len(g.vEdges) == len(opt.graph.vEdges))
			assert(path in opt.PathGroups)
			idx_p = opt.PathGroups.index(path)
			flow = 0
			# aligned to single node or single junction
			if idx_p < opt.offset_hyper:
				 assert(len(path) == 1 or len(path) == 3)
				 if len(path) == 1:
				 	flow = np.sum([opt_results[idx_e] for idx_e in g.vNodes[path[0]].OutEdges])
				 else:
				 	flow = opt_results[path[1]]
			# alignment is a hyper edge
			else:
				flow = opt_results[idx_p + len(g.vEdges) - opt.offset_hyper]
			flows.append(flow)
		# re-assign weight according to the flows
		flows = np.array(flows)
		if np.sum(flows) != 0:
			flows /= np.sum(flows)
		self.Weights = self.Count * flows

	def Print(self):
		s = "("+ ",".join(self.GeneIDs) + ")\t" 
		s += "["
		for fs in self.Aligned_Nodes:
			s += "("+ ",".join(str(x) for x in fs) + "),"
		s = s[:-1] + "]\t"
		s += "["
		for fs in self.Aligned_Edges:
			s += "(" + ",".join(str(x) for x in fs) + "),"
		s = s[:-1] + "]\t"
		s += "["
		for p in self.Aligned_Paths:
			s += "(" + ",".join(str(x) for x in p) + "),"
		s = s[:-1] + "]\t"
		s += str(self.Count) + "\t"
		s += "["+ ",".join(str(x) for x in self.Weights) +"]"
		return s


def ProcessEquivalentClass(salmon_folder, Transcripts, Graphs, GraphNameIndex, out_eqclass):
	TPM = ReadSalmonTPM(salmon_folder + "/quant.sf")
	map_transtuple_weight = GetSalmonAssignWeight(salmon_folder + "/aux_info/eq_classes.txt", TPM)
	ras = []
	tmpreads = []
	# reading sam file
	if Path(salmon_folder + "/mapping.bam").exists():
		fp = pysam.AlignmentFile(salmon_folder + "/mapping.bam")
	else:
		fp = pysam.AlignmentFile(salmon_folder + "/mapping.sam")
	for read in fp:
		if len(tmpreads) == 0 or read.query_name == tmpreads[-1].query_name:
			tmpreads.append(read)
		elif read.query_name != tmpreads[-1].query_name:
			assert(len(tmpreads) > 0)
			assert(len(set([r.query_name for r in tmpreads])) == 1)
			# retrieve the salmon weight
			transtuple = list(set([r.reference_name for r in tmpreads]))
			transtuple.sort()
			transtuple = tuple(transtuple)
			assert(transtuple in map_transtuple_weight)
			weight = map_transtuple_weight[transtuple]
			# Finding the aligned nodes and edges in the graph
			print(tmpreads[0].query_name)
			tmp = ReadAlignment_t(tmpreads, Transcripts, transtuple, weight)
			oldlen = len(tmp.GeneIDs)
			tmp.find_overlapping_nodes(Graphs, GraphNameIndex)
			assert( oldlen == len(tmp.GeneIDs) )
			if len(tmp.GeneIDs) != 0:
				tmp.remove_duplicated_alignment()
				tmp.sort_alignments()
				if len(tmp.GeneIDs) != len(tmp.Weights):
					print(tmp.GeneIDs)
				assert( len(tmp.GeneIDs) == len(tmp.Weights) )
				ras.append(tmp)
			tmpreads = [read]
	# the last read group
	if len(tmpreads) > 0:
		assert(len(set([r.query_name for r in tmpreads])) == 1)
		# retrieve the salmon weight
		transtuple = list(set([r.reference_name for r in tmpreads]))
		transtuple.sort()
		transtuple = tuple(transtuple)
		assert(transtuple in map_transtuple_weight)
		weight = map_transtuple_weight[transtuple]
		# Finding the aligned nodes and edges in the graph
		tmp = ReadAlignment_t(tmpreads, Transcripts, transtuple, weight)
		tmp.find_overlapping_nodes(Graphs, GraphNameIndex)
		if len(tmp.GeneIDs) != 0:
			tmp.remove_duplicated_alignment()
			tmp.sort_alignments()
			assert( len(tmp.GeneIDs) == len(tmp.Weights) )
			ras.append(tmp)
	print("Finish finding graph nodes and edges of the reads.")
	fp.close()
	# # count number of reads with 0 or 1 junction
	# count_1_junction = 0
	# count_2_junction = 0
	# count_1_degree = 0
	# for ra in ras:
	# 	is_1_junction = 0
	# 	is_2_junction = 0
	# 	is_1_degree = 0
	# 	for i in range(len(ra.GeneIDs)):
	# 		if len(ra.Aligned_Edges) <= 1:
	# 			is_1_junction += 1
	# 		elif np.all(np.array([len(Graphs[GraphNameIndex[ra.GeneIDs[i]]].vNodes[v].InEdges) for v in ra.Aligned_Nodes[i]]) == 1) \
	# 			or np.all(np.array([len(Graphs[GraphNameIndex[ra.GeneIDs[i]]].vNodes[v].OutEdges) for v in ra.Aligned_Nodes[i]]) == 1):
	# 			is_1_degree += 1
	# 		if len(ra.Aligned_Edges) <= 2:
	# 			is_2_junction += 1
	# 	if is_1_junction == len(ra.GeneIDs):
	# 		count_1_junction += 1
	# 	elif is_1_junction + is_1_degree == len(ra.GeneIDs):
	# 		count_1_degree += 1
	# 	if is_2_junction == len(ra.GeneIDs):
	# 		count_2_junction += 1
	# print("There are {}% reads with less than 1 junction ({} out of {})".format(100.0*count_1_junction / len(ras), count_1_junction, len(ras)))
	# # There are 91.72455901708847% reads with less than 1 junction (21027214 out of 22924301)
	# merge to equivalent class
	eqclasses = []
	for i in range(len(ras)):
		oldlen = len(ras[i].GeneIDs)
		ras[i].convert2tuple()
		assert(oldlen == len(ras[i].GeneIDs) and len(ras[i].GeneIDs) == len(ras[i].Weights))
	ras.sort(key = lambda x : (x.GeneIDs, x.Aligned_Paths[0]))
	print("Finish sorting reads according to the graph coordinates.")
	with tqdm.tqdm(total=len(ras)) as pbar:
		for i in range(len(ras)):
			if len(eqclasses) == 0 or not eqclasses[-1].IsEqual(ras[i]):
				neweq = EquivalentClass_t(ras[i])
				eqclasses.append(neweq)
			else:
				eqclasses[-1].addread(ras[i])
			if len(eqclasses[-1].GeneIDs) != len(eqclasses[-1].Weights):
				print(i)
				print(eqclasses[-1].GeneIDs)
				print(eqclasses[-1].Weights)
				eqclasses[-1].addread(ras[i])
				print(eqclasses[-1].Weights)
			assert( len(eqclasses[-1].GeneIDs) == len(eqclasses[-1].Weights) )
			pbar.update(1)
	print("Finish merging equivalent classes.")
	for eq in eqclasses:
		eq.Weights = np.array(eq.Weights)
	print("There are {} equivalent classes.".format(len(eqclasses)))
	WriteEquivalentClass(out_eqclass, eqclasses)
	print("Finish writing.")
	return eqclasses


def WriteEquivalentClass(outputfile, eqclasses):
	fp = open(outputfile, 'w')
	fp.write("# GeneIDs\tAligned_Nodes\tAligned_Edges\tReadCount\tWeights\n")
	for eq in eqclasses:
		s = eq.Print()
		fp.write(s + "\n")
	fp.close()


def ReadEquivalentClass(inputfile):
	eqclasses = []
	fp = open(inputfile, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		GeneIDs = [x for x in strs[0][1:-1].split(",") if x != ""]
		Aligned_Nodes = []
		for s in strs[1][2:-2].split("),("):
			tmp_nodes = [int(x) for x in s.split(",") if x != ""]
			tmp_nodes = tuple(tmp_nodes)
			Aligned_Nodes.append(tmp_nodes)
		Aligned_Edges = []
		for s in strs[2][2:-2].split("),("):
			tmp_edges = [int(x) for x in s.split(",") if x != ""]
			tmp_edges = tuple(tmp_edges)
			Aligned_Edges.append(tmp_edges)
		Aligned_Paths = []
		for s in strs[3][2:-2].split("),("):
			tmp_paths = [int(x) for x in s.split(",") if x != ""]
			tmp_paths = tuple(tmp_paths)
			Aligned_Paths.append(tmp_paths)
		ReadCount = int(strs[4])
		Weights = [float(x) for x in strs[5][1:-1].split(",")]
		if np.abs(np.sum(Weights) - ReadCount) >= 1e-4:
			print("Error of Weights: " + str(Weights) +"\t"+ str(ReadCount) +"\t"+ line)
		assert(np.abs(np.sum(Weights) - ReadCount) < 1e-4)
		assert(len(GeneIDs) == len(Aligned_Nodes))
		assert(len(GeneIDs) == len(Aligned_Edges))
		assert(len(GeneIDs) == len(Weights))
		tmp = EquivalentClass_t(None)
		tmp.Count = ReadCount
		tmp.GeneIDs = GeneIDs
		tmp.Aligned_Nodes = Aligned_Nodes
		tmp.Aligned_Edges = Aligned_Edges
		tmp.Aligned_Paths = Aligned_Paths
		tmp.Weights = Weights
		eqclasses.append(tmp)
	fp.close()
	return eqclasses