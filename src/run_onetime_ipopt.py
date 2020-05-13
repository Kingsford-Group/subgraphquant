#!/bin/python

import sys
import numpy as np
from pathlib import Path
from TranscriptClass import *
from TranscriptAlignmentClass import *
from IpoptOptimizeClass import *
from utils import *
from trie_conversion import *
from EffLen_VLMM import *
from flow_graph import FlowGraph
from test_barrier import *
import pyipopt
import tqdm
import copy
import pickle
import concurrent.futures
import time
from func_timeout import func_timeout


def SamplePaths(g, opt):
	fg = FlowGraph(g.edges, opt.results, 0, len(g.nodes)-1)
	tmp_edge_map = {}
	existing_edges = set()
	for e in range(len(g.edges)):
		if e in existing_edges:
			continue
		tmp_edge_list, tmp_log_weights = fg.sample_centered_path(1, e, weighted=True, seed=0)
		tmp_edge_list = tuple(tmp_edge_list[0])
		if not tmp_edge_list in tmp_edge_map:
			tmp_edge_map[tmp_edge_list] = tmp_log_weights[0]
			existing_edges = existing_edges | set(list(tmp_edge_list))
	edge_list = []
	log_weights = []
	for k,v in tmp_edge_map.items():
		edge_list.append(k)
		log_weights.append(v)
	return edge_list, log_weights


def WriteSampleSTPaths(graphs, Opts, old_graphs, old_NameIndex, genomefasta, outputfile):
	genome = ReadGenome(genomefasta)
	fp = open(outputfile, 'w')
	total_paths = 0
	with tqdm.tqdm(total=len(Opts)) as pbar:
		for i in range(len(Opts)):
			opt = Opts[i]
			if np.sum(opt.weights) != 0:
				g = graphs[opt.gid]
				old_g = old_graphs[old_NameIndex[opt.gid]]
				edge_list, log_weights = SamplePaths(g, opt)
				total_paths += len(edge_list)
				assert(len(edge_list) == len(log_weights))
				for j in range(len(edge_list)):
					edges = edge_list[j]
					w = log_weights[j]
					nodes = [g.nodes[g.edges[e][0]][-1] for e in edges] + [g.nodes[g.edges[edges[-1]][1]][-1]]
					node_str =  ""
					seq = ""
					for idx_v in nodes:
						if old_g.vNodes[idx_v].Strand:
							seq += genome[old_g.vNodes[idx_v].Chr][old_g.vNodes[idx_v].StartPos:old_g.vNodes[idx_v].EndPos]
						else:
							seq += ReverseComplement( genome[old_g.vNodes[idx_v].Chr][old_g.vNodes[idx_v].StartPos:old_g.vNodes[idx_v].EndPos] )
						node_str += "(" + str(idx_v) +","+ str(old_g.vNodes[idx_v].Chr) +","+ str(old_g.vNodes[idx_v].StartPos) +","+ str(old_g.vNodes[idx_v].EndPos) +","+ str(old_g.vNodes[idx_v].Strand) +")"
					fp.write(">" + opt.gid + "_sample_" + format(j, '04d') +" "+ str(w) +" "+ node_str + "\n")
					count = 0
					while count < len(seq):
						fp.write(seq[count:min(count+70, len(seq))] + "\n")
						count += 70
			pbar.update(1)
	print("Total paths = {}".format(total_paths))
	fp.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python run_onetime_ipopt.py <salmon_dir> <prefix_graph> <prefix_output> (<threads>)")
	else:
		salmon_dir = sys.argv[1]
		prefix_graph = sys.argv[2]
		prefix_output = sys.argv[3]
		if len(sys.argv) > 4:
			nthreads = int(sys.argv[4])
		else:
			nthreads = 8

		# check whether the first round optimization has already been done
		if Path(prefix_output + "_opt_ipopt_round0.pkl").exists():
			print("No more work need to do.")
			sys.exit()

		graphs, eq_classes = load_data(prefix_graph)
		# efflens = ReadPathEfflen(prefix_graph + "_path_efflen.txt", graphs)
		# assert( sanity_check_pathefflen(efflens, graphs) )
		efflen_hyper = ReadEfflen_hyper(prefix_graph)
		efflen_simple = ReadEfflen_simple(prefix_graph)
		# forbidden_edges = ReadForbiddenEdge_pseudo_eqclass(pseudoeq_filename, graphs)
		forbidden_edges = {}

		Gene_Eq_Index = {gname:[] for gname in graphs.keys()}
		for i in range(len(eq_classes)):
			for g in [ali.gid for ali in eq_classes[i]]:
				Gene_Eq_Index[g].append(i)
		Gene_Eq_Index = {g:list(set(v)) for g,v in Gene_Eq_Index.items() if len(v) > 0}

		Names = list( set(Gene_Eq_Index.keys()) & set(efflen_hyper.keys()) )
		Names.sort()
		Opts = []
		for gname in Names:
			relevant_eq = [eq_classes[j] for j in Gene_Eq_Index[gname]]
			relevant_eq = sum([[z for z in x if z.gid == gname] for x in relevant_eq], [])
			if np.max([x.weights for x in relevant_eq]) < 1e-6:
				continue
			if gname in forbidden_edges:
				print(gname + " contain forbidden edges")
				Opts.append( IpoptObject_split(graphs[gname], efflen_simple[gname], efflen_hyper[gname], [eq_classes[j] for j in Gene_Eq_Index[gname]], forbidden_edges[gname]) )
			else:
				Opts.append( IpoptObject_split(graphs[gname], efflen_simple[gname], efflen_hyper[gname], [eq_classes[j] for j in Gene_Eq_Index[gname]]) )

		# NameIndex = {Opts[i].gid:i for i in range(len(Opts))}
		# pool = concurrent.futures.ProcessPoolExecutor(nthreads)
		# index = [i for i in range(len(Opts)) if np.sum(Opts[i].weights) > 0]
		# future_to_opt = [pool.submit(optimizegraph, Opts[i], max_iter = 500, max_cpu_time = 200) for i in index]
		# with tqdm.tqdm(total=len(index)) as pbar:
		# 	for future in concurrent.futures.as_completed(future_to_opt):
		# 		opt, status = future.result()
		# 		idx = NameIndex[opt.gid]
		# 		Opts[idx].results = opt.results
		# 		pbar.set_postfix(gname=opt.gid, n_var = opt.n_var)
		# 		pbar.update(1)

		with tqdm.tqdm(total = len(Opts)) as pbar:
			for i in range(len(Opts)):
				x0 = initialize_flow(Opts[i], np.sum(Opts[i].weights))
				try:
					x, _ = func_timeout(10, barrier_method, kwargs={"opt":Opts[i], "x0":x0, "stop_criteria":1e-14})
					Opts[i].results = x
				except:
					opt,_ = optimizegraph(Opts[i], max_iter = 3000, max_cpu_time = 100)
					Opts[i] = opt
				pbar.update(1)

		# add back the forbidden edges
		Opts = [opt for opt in Opts if (not (opt.results is None)) and (not np.any(np.isnan(opt.results))) and (not np.any(np.isinf(opt.results)))]
		collected_results = {opt.gid : opt.results for opt in Opts}

		pickle.dump( Opts, open(prefix_output + "_opt_ipopt_round0.pkl", 'wb') )
		pickle.dump( collected_results, open(prefix_output + "_result_ipopt_round0.pkl", 'wb') )

