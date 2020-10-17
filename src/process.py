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
		if codedir == "":
			p = subprocess.Popen("pwd", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			codedir = out.decode('utf-8').strip()
			if codedir[-4:] == "/src":
				codedir = codedir[:-4]
		assert(codedir != "")

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
			if Path(readprefix + "_1.fastq.gz").exists():
				read1 = readprefix + "_1.fastq.gz"
			elif Path(readprefix + "_1.fastq").exists():
				read1 = readprefix + "_1.fastq"
			elif Path(readprefix + "_1.fasta.gz").exists():
				read1 = readprefix + "_1.fasta.gz"
			elif Path(readprefix + "_1.fasta").exists():
				read1 = readprefix + "_1.fasta"
			else:
				print("Error: No (gzipped) fasta or fastq files with such prefix.")
				sys.exit()
			if Path(readprefix + "_2.fastq.gz").exists():
				read2 = readprefix + "_2.fastq.gz"
			elif Path(readprefix + "_2.fastq").exists():
				read2 = readprefix + "_2.fastq"
			elif Path(readprefix + "_2.fasta.gz").exists():
				read2 = readprefix + "_2.fasta.gz"
			elif Path(readprefix + "_2.fasta").exists():
				read2 = readprefix + "_2.fasta"
			else:
				print("Error: No (gzipped) fasta or fastq files with such prefix.")
				sys.exit()
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
		if not Path(outdir_gs + "/gs_graph_fragstart_beforebias.txt").exists():
			print("GENERATING GRAPH...")
			p = subprocess.Popen("mkdir -p "+outdir_gs, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			salmonauxdir = outdir_salmon + "/aux_info"
			salmonquant = outdir_salmon + "/quant.sf"
			outprefix = outdir_gs + "/gs"
			print( "{}/bin/testgraph {} {} {} {} {}".format(codedir, gtffile, genomefasta, salmonauxdir, salmonquant, outprefix) )
			p = subprocess.Popen("{}/bin/testgraph {} {} {} {} {}".format(codedir, gtffile, genomefasta, salmonauxdir, salmonquant, outprefix), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if err != b'':
				print(err)
				sys.exit()

		# generate equivalent class and prefix trie
		if not Path(outdir_gs + "/eq_classes.txt").exists():
			print("GENERATING EQ CLASSES AND PREFIX TRIE...")
			graphfile = outdir_gs + "/gs_graph_fragstart_beforebias.txt"
			# trans_bamfile = outdir_salmon + "/mapping.bam"
			print( "python {}/src/ProcessEqClasses.py {} {} {} {} {}".format(codedir, gtffile, graphfile, outdir_salmon, outdir_gs+"/eq_classes.txt", outdir_gs+"/gs") )
			p = subprocess.Popen("python {}/src/ProcessEqClasses.py {} {} {} {} {}".format(codedir, gtffile, graphfile, outdir_salmon, outdir_gs+"/eq_classes.txt", outdir_gs+"/gs"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if (b'ERROR' in err) or (b'Error' in err):
				print(err)
				sys.exit()

		# path effective length
		if not Path(outdir_gs + "/gs_simpleedge_efflen.txt").exists():
			print("CALCULATING PATH EFFECTIVE LENGTH...")
			salmonquant = outdir_salmon + "/quant.sf"
			outprefix = outdir_gs + "/gs"
			print( "{}/bin/pathbias {} {} {} {}".format(codedir, gtffile, genomefasta, salmonquant, outprefix) )
			p = subprocess.Popen("{}/bin/pathbias {} {} {} {}".format(codedir, gtffile, genomefasta, salmonquant, outprefix), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if (b'ERROR' in err) or (b'Error' in err):
				print(err)
				sys.exit()

		# estimating flow
		if not Path(outdir_gs + "/gs_result_ipopt_round0.pkl").exists():
			print("ESTIMATING PREFIX GRAPH EDGE FLOW...")
			print("python {}/src/run_onetime_ipopt.py {} {} {}".format(codedir, outdir_salmon, outdir_gs + "/gs", outdir_gs + "/gs"))
			p = subprocess.Popen("python {}/src/run_onetime_ipopt.py {} {} {}".format(codedir, outdir_salmon, outdir_gs + "/gs", outdir_gs + "/gs"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if (b'ERROR' in err) or (b'Error' in err):
				print(err)
				sys.exit()

		# lb and ub under the assumption that all S-T paths freely express
		if not Path(outdir_gs + "/gs_maxflow_bound.txt").exists():
			print("BOUNDING TRANSCRIPTS: IPOPT + ALL PATH")
			p = subprocess.Popen("python {}/src/BoundingTranscriptFlows.py {} {} {}".format(codedir, outdir_gs + "/gs", outdir_gs + "/gs_result_ipopt_round0.pkl", outdir_gs + "/gs_maxflow_bound.txt"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if (b'ERROR' in err) or (b'Error' in err):
				print(err)
				sys.exit()

		# projecting Salmon quantification result on the edge flow for running LP
		if not Path(outdir_gs + "/salmon_result.pkl").exists():
			print("PROJECTING SALMON QUANTIFICATION ON GRAPH")
			p = subprocess.Popen("python {}/src/GetFlow_quantifier.py {} {} {}".format(codedir, outdir_gs + "/gs", outdir_salmon + "/quant.sf", outdir_gs + "/salmon_result.pkl"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if (b'ERROR' in err) or (b'Error' in err):
				print(err)
				sys.exit()

		# lb and ub under the assumption that reference is complete
		if not Path(outdir_gs + "/salmon_lp_bound.txt").exists():
			print("BOUNDING TRANSCRIPTS: COMPLETE REFERENCE ASSUMPTION")
			p = subprocess.Popen("python {}/src/lp_fixed_transcripts.py {} {} {}".format(codedir, outdir_gs + "/gs", outdir_gs + "/salmon_result.pkl", outdir_gs + "/salmon_lp_bound.txt"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if (b'ERROR' in err) or (b'Error' in err):
				print(err)
				sys.exit()

