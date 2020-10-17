#!/bin/python

import sys
import numpy as np


class Node_t(object):
	def __init__(self, _ID, _Chr, _StartPos, _EndPos, _Strand, _IsSegment, _InEdges, _OutEdges, _BiasMultiplier):
		self.ID = _ID
		self.Chr = _Chr
		self.StartPos = _StartPos
		self.EndPos = _EndPos
		self.Strand = _Strand
		self.IsSegment = _IsSegment
		self.InEdges = _InEdges
		self.OutEdges = _OutEdges
		self.BiasMultiplier = _BiasMultiplier
		# initialize or compute other attributes
		self.Length = self.EndPos - self.StartPos
		if self.Length <= 0:
			self.Length = 1
		self.Coverage = 0


class Edge_t(object):
	def __init__(self, _ID, _Node_1, _Node_2, _IncidentTranscripts):
		self.ID = _ID
		self.Node_1 = _Node_1
		self.Node_2 = _Node_2
		self.IncidentTranscripts = _IncidentTranscripts
		self.Flow = 0


class GeneGraph_t(object):
	def __init__(self, _GeneID, _vNodes, _vEdges):
		self.GeneID = _GeneID
		self.vNodes = _vNodes
		self.vEdges = _vEdges

	def find_containing_nodes(self, regions):
		containing_nodes = []
		for r in regions:
			for v in self.vNodes:
				if r[0] <= v.StartPos and r[1] >= v.EndPos:
					containing_nodes.append(v.ID)
		containing_nodes = list(set(containing_nodes))
		containing_nodes.sort()
		return containing_nodes

	def get_node_flow(self):
		node_flows = []
		for v in self.vNodes:
			node_flows.append( np.sum(self.vEdges[i].Flow for i in v.OutEdges) )
		node_flows = np.array(node_flows)
		node_flows[0] = 0
		node_flows.append(0)
		return np.array(node_flows)

	def WriteGraph(self):
		s = ""
		# string for nodes
		for v in self.vNodes:
			s += "node\t" + str(v.ID) +"\t"+ v.Chr +"\t"+ str(v.StartPos) +"\t"+ str(v.EndPos) +"\t"
			if v.Strand:
				s += "1\t"
			else:
				s += "0\t"
			s += (str(v.Coverage)) +"\t"+ str(v.BiasMultiplier) +"\t("+ ",".join([str(e) for e in v.InEdges]) +")\t(" + ",".join([str(e) for e in v.OutEdges]) +")\n"
		# string for edges
		for e in self.vEdges:
			s += "edge\t"+ str(e.ID) +"\t"+ str(e.Node_1) +"\t"+ str(e.Node_2) +"\t("+ ",".join(e.IncidentTranscripts) +")\n"
		return s

	def enumerate_paths(self, max_num_paths = 1000):
		# paths encoded by the edge index
		all_paths = [[e.ID] for e in self.vEdges if e.Node_1 == 0]
		# keep record of how many paths have ended with node t
		num_ending_properly = np.sum([self.vEdges[x[-1]].Node_2 == self.vNodes[-1].ID for x in all_paths])
		while num_ending_properly < len(all_paths) and num_ending_properly < max_num_paths:
			extended_all_paths = []
			for i in range(len(all_paths)):
				if self.vEdges[all_paths[i][-1]].Node_2 != self.vNodes[-1].ID:
					this_path = all_paths[i]
					# extending this path
					v = self.vEdges[this_path[-1]].Node_2
					for e in self.vNodes[v].OutEdges:
						extended_all_paths.append( this_path + [e] )
				else:
					extended_all_paths.append( all_paths[i] )
			all_paths = extended_all_paths
			num_ending_properly = np.sum([self.vEdges[x[-1]].Node_2 == self.vNodes[-1].ID for x in all_paths])
		index_ending_properly = [i for i in range(len(all_paths)) if self.vEdges[all_paths[i][-1]].Node_2 == self.vNodes[-1].ID]
		assert(len(index_ending_properly) == len(all_paths) or len(index_ending_properly) >= max_num_paths)
		return [all_paths[i] for i in index_ending_properly]


def ReadGraphFile(filename):
	Graphs = []
	fp = open(filename, 'r')
	_GeneID = ""
	_vNodes = []
	_vEdges = []
	numnodes = 0
	numedges = 0
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[0] == "graph":
			if _GeneID != "":
				g = GeneGraph_t(_GeneID, _vNodes, _vEdges)
				Graphs.append(g)
			_GeneID = strs[1]
			numnodes = int(strs[2])
			numedges = int(strs[3])
			_vNodes = []
			_vEdges = []
		elif strs[0] == "node":
			_ID = int(strs[1])
			_Chr = strs[2]
			_StartPos = int(strs[3])
			_EndPos = int(strs[4])
			_Strand = (strs[5] == "1")
			_BiasMultiplier = float(strs[7])
			_InEdges = [int(x) for x in strs[8][1:-1].split(",") if x != ""]
			_OutEdges = [int(x) for x in strs[9][1:-1].split(",") if x != ""]
			n = Node_t(_ID, _Chr, _StartPos, _EndPos, _Strand, (_StartPos != _EndPos), _InEdges, _OutEdges, _BiasMultiplier)
			_vNodes.append(n)
		elif strs[0] == "edge":
			_ID = int(strs[1])
			_Node_1 = int(strs[2])
			_Node_2 = int(strs[3])
			_IncidentTranscripts = strs[4][1:-1].split(",")
			e = Edge_t(_ID, _Node_1, _Node_2, _IncidentTranscripts)
			_vEdges.append(e)
		else:
			print("wrong key")
	# append the last graph
	if _GeneID != "":
		g = GeneGraph_t(_GeneID, _vNodes, _vEdges)
		Graphs.append(g)
	fp.close()
	return Graphs


def WriteGraphFile(outputfile, Graphs):
	fp = open(outputfile, 'w')
	fp = open(outputfile, 'w')
	fp.write("# graph\tGeneID\tNumNodes\tNumEdges\n")
	fp.write("# node\tID\tChr\tStartPos\tEndPos\tStrand\tCoverage\tBiasMultiplier\tInEdges\tOutEdges\n")
	fp.write("# edge\tID\tNode_1\tNode_2\tIncidentTranscripts\n")
	for g in Graphs:
		fp.write("graph\t{}\t{}\t{}\n".format(g.GeneID, len(g.vNodes), len(g.vEdges)))
		s = g.WriteGraph()
		fp.write(s)
	fp.close()
