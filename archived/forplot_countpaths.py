#!/bin/python
# For plotting purpose. This script shows the num splicing edges, num annotated transcripts, and num total paths for human genes with gencode v26 annotation. 

import sys
import numpy as np
import tqdm
from GeneGraph import *
from utils import *
from TranscriptClass import *

# for human
old_graphs = ReadGraphFile("/home/congm1/savanna/savannacong33/GraphSalmonSimu/salmon/nb_Full_GEU_new/gs_graph_fragstart.txt")
old_NameIndex = {old_graphs[i].GeneID:i for i in range(len(old_graphs))}

num_edges = [len(g.vEdges) for g in old_graphs]
num_annotated = [0] * len(old_graphs)
num_enumerate = [0] * len(old_graphs)

# annotated
for i in range(len(old_graphs)):
	g = old_graphs[i]
	transnames = sum([e.IncidentTranscripts for e in g.vEdges], [])
	num_annotated[i] = len(set(transnames))

# enumeration
batch_size = 5000
num_batches = int(np.ceil(1.0 * len(old_graphs) / batch_size))
for b in range(num_batches):
	s = batch_size*b
	t = batch_size*(b+1)
	transpaths = ExpandReference( old_graphs[s:t], max_num_paths_pergene=1e7)
	for info in transpaths:
		idx = old_NameIndex[info[0]]
		num_enumerate[idx] += 1


assert(np.all(np.array(num_annotated) > 0))
assert(np.all(np.array(num_enumerate) > 0))


fp = open("/home/congm1/savanna/savannacong33/GraphSalmonSimu/result_general/human_num_edge_path.txt", 'w')
fp.write("# GeneID\tSplicingEdges\tAnnotated\tTotalPaths\n")
for i in range(len(old_graphs)):
	fp.write("{}\t{}\t{}\t{}\n".format(old_graphs[i].GeneID, num_edges[i], num_annotated[i], num_enumerate[i]))
fp.close()


# for fruit fly
old_graphs = ReadGraphFile("/home/congm1/savanna/savannacong33/fruitfly/Drosophila_melanogaster.BDGP6.95.chr.graph_fragstart_beforebias.txt")
old_NameIndex = {old_graphs[i].GeneID:i for i in range(len(old_graphs))}

num_edges = [len(g.vEdges) for g in old_graphs]
num_annotated = [0] * len(old_graphs)
num_enumerate = [0] * len(old_graphs)

# annotated
for i in range(len(old_graphs)):
	g = old_graphs[i]
	transnames = sum([e.IncidentTranscripts for e in g.vEdges], [])
	num_annotated[i] = len(set(transnames))

# enumeration
batch_size = 5000
num_batches = int(np.ceil(1.0 * len(old_graphs) / batch_size))
for b in range(num_batches):
	s = batch_size*b
	t = batch_size*(b+1)
	transpaths = ExpandReference( old_graphs[s:t], max_num_paths_pergene=1e7)
	for info in transpaths:
		idx = old_NameIndex[info[0]]
		num_enumerate[idx] += 1


assert(np.all(np.array(num_annotated) > 0))
assert(np.all(np.array(num_enumerate) > 0))


fp = open("/home/congm1/savanna/savannacong33/GraphSalmonSimu/result_general/fruitfly_num_edge_path.txt", 'w')
fp.write("# GeneID\tSplicingEdges\tAnnotated\tTotalPaths\n")
for i in range(len(old_graphs)):
	fp.write("{}\t{}\t{}\t{}\n".format(old_graphs[i].GeneID, num_edges[i], num_annotated[i], num_enumerate[i]))
fp.close()
