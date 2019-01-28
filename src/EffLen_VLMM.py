#!/bin/python

import numpy as np
import struct
from GeneGraph import *


def ReadCorrectedLength(filename):
	corrections={}
	fp=open(filename, 'rb')
	numtrans=struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen=struct.unpack('i', fp.read(4))[0]
		seqlen=struct.unpack('i', fp.read(4))[0]
		name=""
		correction=np.zeros(seqlen)
		for j in range(namelen):
			name+=struct.unpack('c', fp.read(1))[0].decode('utf-8')
		for j in range(seqlen):
			correction[j]=struct.unpack('d', fp.read(8))[0]
		corrections[name]=correction
	fp.close()
	print("Finish reading theoretical distributions for {} transcripts.".format(len(corrections)))
	return corrections


def ReadSalmonQuant(filename):
	Exp = {}
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		Exp[strs[0]] = float(strs[3])
	fp.close()
	return Exp


# deprecated
def TranscriptPaths(oldgraphfile, new_graphs):
	TransPaths = [] #(gname, tname, new edges, length of nodes before each edge)
	old_graphs = ReadGraphFile(oldgraphfile)
	for g in old_graphs:
		gname = g.GeneID
		trans = sum([e.IncidentTranscripts for e in g.vEdges], [])
		trans = list(set(trans))
		for tname in trans:
			old_edges = [e.ID for e in g.vEdges if tname in e.IncidentTranscripts]
			old_nodes = [g.vNodes[g.vEdges[idx_e].Node_1].ID for idx_e in old_edges] + [g.vNodes[-1].ID]
			node_lengths = [g.vNodes[idx_v].EndPos - g.vNodes[idx_v].StartPos for idx_v in old_nodes]
			new_nodes, new_edges = new_graphs[gname].walk(old_nodes)
			assert(len(old_nodes) == len(new_nodes))
			assert(len(node_lengths) == len(new_edges) + 1)
			TransPaths.append( (gname, tname, new_edges, node_lengths[:-1]) )
	return TransPaths


def CalculateOldNodeEffLen(oldgraphfile, Exp, corrections):
	old_graphs = ReadGraphFile(oldgraphfile)
	old_nodes_efflen = {g.GeneID:np.zeros(len(g.vNodes)) for g in old_graphs}
	for g in old_graphs:
		gname = g.GeneID
		# sum of correct length * TPM per node
		numer = np.zeros(len(g.vNodes))
		# sum of TPM per node
		denom = np.zeros(len(g.vNodes)) 
		# loop over transcripts
		trans = sum([e.IncidentTranscripts for e in g.vEdges], [])
		trans = list(set(trans))
		for tname in trans:
			old_edges = [e.ID for e in g.vEdges if tname in e.IncidentTranscripts]
			old_nodes = [g.vNodes[g.vEdges[idx_e].Node_1].ID for idx_e in old_edges] + [g.vNodes[-1].ID]
			node_lengths = [g.vNodes[idx_v].EndPos - g.vNodes[idx_v].StartPos for idx_v in old_nodes]
			assert(tname in corrections)
			correction = corrections[tname]
			this_TPM = 1e-8
			if tname in Exp:
				this_TPM = max(Exp[tname], 1e-8)
			cur_length = 0
			for i in range(len(old_nodes)):
				numer[old_nodes[i]] += this_TPM * np.sum(correction[cur_length:(cur_length + node_lengths[i])])
				denom[old_nodes[i]] += this_TPM
				cur_length += node_lengths[i]
		if not np.all(denom > 0):
			print(g)
			print(denom)
			print(trans)
		assert(np.all(denom > 0))
		old_nodes_efflen[gname] = numer / denom
	return old_nodes_efflen


# def CalculateEdgewiseEffLen(new_graphs, TransPaths, Exp, corrections):
# 	efflen = {}
# 	efflen_numer = {gname:np.zeros(len(g.edges)) for gname,g in new_graphs.items()}
# 	efflen_denom = {gname:np.zeros(len(g.edges)) for gname,g in new_graphs.items()}
# 	for info in TransPaths:
# 		(gname, tname, new_edges, node_lengths) = info
# 		assert(len(new_edges) == len(node_lengths))
# 		TPM = 1e-8
# 		if tname in Exp:
# 			TPM = max(Exp[tname], 1e-8)
# 		assert(tname in corrections)
# 		correction = corrections[tname]
# 		assert(np.sum(node_lengths) == len(correction))
# 		cur_length = 0
# 		for i in range(len(new_edges)):
# 			efflen_numer[gname][new_edges[i]] += np.sum(correction[cur_length:(cur_length + node_lengths[i])]) * TPM
# 			efflen_denom[gname][new_edges[i]] += TPM
# 			cur_length += node_lengths[i]
# 	for gname, numer in efflen_numer.items():
# 		denom = efflen_denom[gname]
# 		# avoid 0 denominator
# 		denom[np.where(denom == 0)[0]] = 1
# 		value = numer / denom
# 		efflen[gname] = value
# 	return efflen


def CalculateEdgewiseEffLen(new_graphs, old_nodes_efflen):
	efflen = {}
	for gname, g in new_graphs.items():
		assert(gname in old_nodes_efflen)
		this_old_nodes_efflen = old_nodes_efflen[gname]
		corres_nodes = [g.nodes[e[0]][-1] for e in g.edges]
		efflen[gname] = np.array([this_old_nodes_efflen[idx_v] for idx_v in corres_nodes])
	return efflen


def UpdateEdgewiseEfflen(new_graphs, old_graphs, corrections, transnodes, log_weights, efflen, first_node=True, damping_factor = 0.1):
	old_NameIndex = {old_graphs[i].GeneID:i for i in range(len(old_graphs))}
	new_weighted_efflen = {}
	for info in transnodes:
		(gname,tname,old_nodes) = info
		old_g = old_graphs[old_NameIndex[gname]]
		g = new_graphs[gname]
		if not (tname in corrections):
			continue
		bias = corrections[tname]
		w = log_weights[tname]
		# the efflen of each old nodes
		node_efflen = np.zeros(len(old_nodes))
		covered_length = 0
		for idx_v in range(len(node_efflen)):
			assert(old_nodes[idx_v] <= len(old_g.vNodes))
			v = old_g.vNodes[old_nodes[idx_v]]
			assert(covered_length + v.EndPos - v.StartPos <= len(bias))
			node_efflen[idx_v] += np.sum(bias[covered_length:(covered_length + v.EndPos - v.StartPos)])
			covered_length += v.EndPos - v.StartPos
		assert(covered_length >= len(bias))
		new_nodes, new_edges = g.walk(old_nodes)
		assert(len(new_edges) == len(old_nodes) - 1)
		for idx_e in range(len(new_edges)):
			# XXX: the following can be changed when edge efflen is the second nodes's efflen
			if (gname, new_edges[idx_e]) in new_weighted_efflen:
				if first_node:
					new_weighted_efflen[ (gname, new_edges[idx_e]) ].append( (node_efflen[idx_e], w) )
				else:
					new_weighted_efflen[ (gname, new_edges[idx_e]) ].append( (node_efflen[idx_e+1], w) )
			else:
				if first_node:
					new_weighted_efflen[ (gname, new_edges[idx_e]) ] = [ (node_efflen[idx_e], w) ]
				else:
					new_weighted_efflen[ (gname, new_edges[idx_e]) ] = [ (node_efflen[idx_e+1], w) ]
	# update the existing efflen
	new_efflen = {}
	for gname,v in efflen.items():
		new_efflen[gname] = copy.copy(v)
	for key, value in new_weighted_efflen.items():
		gname, idx_e = key
		updates = np.array([x[0] for x in value])
		log_w = np.array([x[1] for x in value])
		w = np.exp(log_w)
		old_l = efflen[gname][idx_e]
		new_l = np.dot(updates, w) / np.sum(w)
		new_efflen[gname][idx_e] = damping_factor * old_l + (1 - damping_factor) * new_l
	return new_efflen

