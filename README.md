# Overview
This repo contain the scripts to quantify the uncertainty of expression estimation for subgraphs. Due to the non-identifiability of RNA-seq expression quantification, there can be multiple configurations of transcript abundances equally optimal in terms of the likelihood. Lower bound and upper bound of abundances of transcripts or a general subgraph can be calculated using this repo.

Three major modules are included:
+ Uncertainty quantification module.
+ Flow estimation module, for preparing the graph edge abundances required in the subgraph uncertainty quantification.
+ Analysis module for the 30 GEUVADIS samples.


# Uncertainty quantification module
The degree of uncertainty is represented by the lower bound and upper bound of the set of optimal abundances of the queried subgraph.

Subgraph abundance query can be in either disjunctive (OR) or conjunctive (AND) form. 
+ Disjunctive: Abundance of a disjunctive subgraph is the sum of abundances of RNA molecules that contain at least one edge in the subgraph. 
+ Conjunctive: Abundance of a conjunctive subgraph is the sum of abundances of RNA molecules that contain all the edges in the subgraph. We extend of the usage of the conjective quantification, so that given a list of edge lists, an edge must be contained for each edge list.

## Prerequisite
+ python3 (>=3.4.3)
+ pickle
+ [networkx](https://networkx.github.io/documentation/stable/)
+ [numpy](http://www.numpy.org/)
+ [tqdm](https://tqdm.github.io/)

## Usage
src/flow_graph.py contains the main functions to calculate the lower bound and upper bound of subgraph abundances. Other packages and class definitions need to be imported for preparing the input of the bound calculation.
```
from GeneGraph import *
from flow_graph import FlowGraph
from trie_conversion import *
```

The conjunctive and disjunctive subgraph abundance bounds can be calculated using the following python statements
```
# g_prefix: an instance of PrefixTrie class
# g_splice: the instance of GeneGraph class of the same gene as g_prefix
# flows: the estimated edge abundance of g_prefix
# edge_list: the list of edges in OR (disjunctive) quantification. Edges are indexed in g_prefix
# edge_lists: the list of edge lists in AND (conjunctive) quantification. Edges are indexed in g_prefix, and are required to be in topological order.
# return: a scalar denote the lower / upper bound of the query

# constructing a FlowGraph object, this object will be transformed for calculating max-flow
fg = FlowGraph(g_prefix.edges, flows, 0, len(g_prefix.nodes)-1)

# Lower bound for OR-quantification
f = fg.diff_flow(edge_list)
# Upper bound for OR-quantification
f = fg.split_flow(edge_list)

# Lower bound for AND-quantification
f = split_block_flow(edge_lists)
# Upper bound for AND-Quantifications
f = fg.block_flow(edge_lists)
```
Note that bounds can have floating point number inaccuracies.


# Flow estimation module
To estimate the flow of an prefix graph, the following steps are calculated sequentially: generating the bias correction model using Salmon, construct a splice graph from annotation GTF file, estimating the node biases using Salmon bias correction model, projecting the read alignment onto the graph, estimating flow and assigning multi-mapped reads using EM algorithm. These steps are runned using the script src/process.py

## Prerequisite
+ [Salmon](https://salmon.readthedocs.io/en/latest/)
+ [Boost](https://www.boost.org/)
+ [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
+ [Jellyfish](https://github.com/gmarcais/Jellyfish)
+ [Spline](https://kluge.in-chemnitz.de/opensource/spline/)
+ python3 (>=3.4.3)
+ pickle
+ [numpy](http://www.numpy.org/)
+ [tqdm](https://tqdm.github.io/)
+ [pysam](https://pysam.readthedocs.io/en/latest/)
+ [pyipopt](https://github.com/xuy/pyipopt)
+ [Ipopt](https://projects.coin-or.org/Ipopt)

## Compile
The code for processing Salmon bias correction model is written in c++, which depends on Boost, Eigen, Jellyfish, and Spline library. Download the library from the website indicated above. Among these libraries, only Boost need to be compiled, while the others are header-only.

After obtaining the libraries and compiling Boost, compile the bias correction processing code using
```
cd subgraphquant
make (BOOST=) (EIGEN=) (JELLYFISH=) (SPLINE=)
```

## Usage
```
python process.py <readprefix> <salmonindex> <outdir_salmon> <outdir_flow> <gtffile> <genomefasta>
```


# Analysis module for the 30 GEUVADIS samples

## Data accession
The SRA accession of the 30 GEUVADIS RNA-seq samples can be found in the last column of file data/Metadata.txt. The following command will download the gziped fastq files in the current directory.
```
./data/getdata.sh
```
The reference genome, transcriptome, and annotation can be downloaded from [Gencode website](https://www.gencodegenes.org/). Version 26 of human is used in the manuscript. Salmon is used for learning the bias correction model. Using the following command to download the required reference and generate Salmon index
```
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz | gunzip > gencode.v26.annotation.gtf
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.p10.genome.fa.gz | gunzip > GRCh38.p10.genome.fa
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz | gunzip > gencode.v26.transcripts.fa
salmon index --gencode -t gencode.v26.transcripts.fa -i gencode.v26.full
```

## quantifying the uncertainty for annotated transcripts
After the flow is estimated by src/process.py, src/BoundingTranscriptFlows.py will calculate the lower and upper bounds of abundance estimation of the annotated transcripts.
```
python BoundingTranscriptFlows.py <graph_prefix> <result_pickle> <outputfile>
```

## Meta-script to process from the beginning to the end
The script scripts/metarun.sh will start from downloading the data, and process to the end of quantifying the abundance uncertainty. A folder named GEUVADIS will be created under the current directory, which stores the output.
```
./scripts/metarun.sh
```

## plotting the results
Using the results generated by scripts.metarun.sh script, the following R script will read from the uncertainty bounding outputs to plot the figures that are used in the manuscript. The figures are located in the directory where metarun.sh is called.
```
Rscript scripts/draw_figure.R <folder where metarun.sh is called>
```
