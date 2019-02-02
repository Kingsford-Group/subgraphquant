#!/bin/python

import sys
import subprocess
from pathlib import Path
from trie_conversion import *

if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python process.py <readprefix> <salmonindex> <outdir_salmon> <outdir_flow> <gtffile> <genomefasta>")
	else:
		codedir = "/".join(sys.argv[0].split("/")[:-2])

		readprefix = sys.argv[1]
		salmonindex = sys.argv[2]
		outdir_salmon = sys.argv[3]
		outdir_gs = sys.argv[4]
		gtffile = sys.argv[5]
		genomefasta = sys.argv[6]

		# salmon quantify
		if not Path(outdir_salmon + "/quant.sf").exists():
			print("RUNNING SALMON...")
			p = subprocess.Popen("mkdir -p "+outdir_salmon, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			read1 = readprefix + "_1.fastq.gz"
			if not Path(readprefix + "_1.fastq.gz").exists() and Path(readprefix + "_1.fasta.gz").exists():
				read1 = readprefix + "_1.fasta.gz"
			read2 = readprefix + "_2.fastq.gz"
			if not Path(readprefix + "_2.fastq.gz").exists() and Path(readprefix + "_2.fasta.gz").exists():
				read2 = readprefix + "_2.fasta.gz"
			p = subprocess.Popen("salmon quant -l A -p 4 -i {} -1 {} -2 {} --gcBias --seqBias --posBias --dumpEqWeights -o {} --writeMappings={}".format(salmonindex, read1, read2, outdir_salmon, outdir_salmon+"/mapping.sam"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if b"ERROR" in err:
				print(err)
				sys.exit()

		# convert sam to bam
		if not Path(outdir_salmon + "/mapping.bam").exists():
			print("CONVERTING SAM TO BAM...")
			print("samtools view -Shb {} -o {}".format(outdir_salmon+"/mapping.sam", outdir_salmon+"/mapping.bam"))
			p = subprocess.Popen("samtools view -Shb {} -o {}".format(outdir_salmon+"/mapping.sam", outdir_salmon+"/mapping.bam"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if err != b'':
				print(err)
				sys.exit()

		# generate graph
		if not Path(outdir_gs + "/gs_graph_fragstart.txt").exists() or not Path(outdir_gs + "/gs_corrections_fragstart.dat").exists():
			print("GENERATING GRAPH...")
			p = subprocess.Popen("mkdir -p "+outdir_gs, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			salmonauxdir = outdir_salmon + "/aux_info"
			salmonquant = outdir_salmon + "/quant.sf"
			outprefix = outdir_gs + "/gs"
			p = subprocess.Popen("{}/bin/testgraph 1 {} {} {} {} {}".format(codedir, gtffile, genomefasta, salmonauxdir, salmonquant, outprefix), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if err != b'':
				print(err)
				sys.exit()

		# generate equivalent class and prefix trie
		if not Path(outdir_gs + "/eq_classes.txt").exists():
			print("GENERATING EQ CLASSES AND PREFIX TRIE...")
			graphfile = outdir_gs + "/gs_graph_fragstart.txt"
			trans_bamfile = outdir_salmon + "/mapping.bam"
			p = subprocess.Popen("python {}/src/ProcessEqClasses.py {} {} {} {} {}".format(codedir, gtffile, graphfile, trans_bamfile, outdir_gs+"/eq_classes.txt", outdir_gs+"/gs"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if (b'ERROR' in err) or (b'Error' in err):
				print(err)
				sys.exit()

		# deal with the bug of old version of ProcessEqClasses.py
		#if not Path(outdir_gs + "/gs_prefix_tries.dat").exists():
		#	print("GENERATING PREFIX TRIE...")
		#	work(outdir_gs + "/gs_graph_fragstart.txt", outdir_gs + "/eq_classes.txt", outdir_gs + "/gs", debug=False)

		# estimating flow
		if not Path(outdir_gs + "/gs_result_ipopt_round0.pkl").exists():
			print("ESTIMATING PREFIX GRAPH EDGE FLOW...")
			p = subprocess.Popen("python {}/src/run_onetime_ipopt.py {} {} {}".format(codedir, outdir_salmon, outdir_gs + "/gs", outdir_gs + "/gs"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if (b'ERROR' in err) or (b'Error' in err):
				print(err)
				sys.exit()
