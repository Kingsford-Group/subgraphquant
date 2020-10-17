# Overview
This repo contain the scripts to perform non-identifiability aware transcript expression quantification based on [Salmon](https://salmon.readthedocs.io/en/latest/). The optimization problem for expression quantification may have multiple optimal solutions, which is known as non-identifiability problem. The multiple optimal abundances for each transcript is a continuous interval due to the convexity of the optimization problem. This repo provides the codes to calculate the range of optimal transcript abundances under any of the assumptions (1) the expression all comes from the reference transcriptome; (2) the expression comes from any of the splice graph paths; (3) a certain proportion of expression comes from the reference transcriptome and the rest comes from any of the splice graph paths, where the proportion is a user-defined argument.

Three major modules are included:
+ Non-identifiability aware transcript quantification module.
+ Analysis module for 16 Human Body Map samples.
+ Analysis module for 6 MCF10 cell line samples.


# Non-identifiability aware transcript quantification module

## Prerequisite
+ python3 (>=3.4.3)
+ pickle
+ [networkx](https://networkx.github.io/documentation/stable/)
+ [numpy](http://www.numpy.org/)
+ [tqdm](https://tqdm.github.io/)
+ [pysam](https://pysam.readthedocs.io/en/latest/)
+ [pyipopt](https://github.com/xuy/pyipopt)
+ [cvxopt](http://cvxopt.org/userguide/intro.html)
+ [Salmon](https://salmon.readthedocs.io/en/latest/)
+ [Boost](https://www.boost.org/)
+ [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
+ [Jellyfish](https://github.com/gmarcais/Jellyfish)
+ [Spline](https://kluge.in-chemnitz.de/opensource/spline/)
+ [Ipopt](https://projects.coin-or.org/Ipopt)
+ [GLPK](https://www.gnu.org/software/glpk/)
+ [progress-cpp](https://github.com/prakhar1989/progress-cpp)

## Compile
The code for processing Salmon bias correction model is written in c++, which depends on Boost, Eigen, Jellyfish, and Spline library. Download the library from the website indicated above. Among these libraries, only Boost need to be compiled, while the others are header-only.

After obtaining the libraries and compiling Boost, compile the bias correction processing code using
```
cd subgraphquant
make (BOOST=) (EIGEN=) (JELLYFISH=) (SPLINE=) (PBAR=)
```

## Usage
Given the (gzipped) fastq/fasta files of an RNA-seq sample, the following script will 
+ perform Salmon quantification; 
+ calculate the range of optimal abundances under (1) complete reference assumption and (2) incomplete reference but correct splice graph assumption.
```
python process.py <readprefix> <salmonindex> <outdir_salmon> <outdir_flow> <gtffile> <genomefasta>
```
Input specification:

Input name    | Type   | Description
---           | :---:  | ---
readprefix    | string | Prefix of the paired-end RNA-seq fastq or fastq.gz files, including path. For example, the argument for paired-end files `/home/sample_1.fastq.gz` and `/home/sample_2.fastq.gz` should be `/home/sample`.
salmonindex   | string | Path to the directory of Salmon indexes. The indexes should be pre-computed. See [Salmon](https://salmon.readthedocs.io/en/latest/) documentation for performing the indexing.
outdir_salmon | string | Path to the directory to store Salmon quantification result. Salmon quantification will be performed by the script of `process.py`.
outdir_flow   | string | Path to the directory to store the computed range of optima for all transcripts.
gtffile       | string | Path to the GTF annotation file. The annotation file should have the same set of transcripts (with the same id) as in Salmon indexes.
genomefasta   | string | Path to the genome fasta file.

(TODO, add a script to take in the proportion as argument and interpolate the two ranges for the third assumption)

## Output specification
The main output are the follows:
+ `<outdir_flow>/salmon_lp_bound.txt`: range of optimal abundances under the assumption of complete reference.
+ `<outdir_flow>/gs_maxflow_bound.txt`: range of optimal abundances under the assumption of incomplete reference but correct splice graph.

Both files have the following 3 columns:

Column name | Description
---         | ---
Name        | Name of the transcript
lb          | Lower bound of the range of optimal abundances, normalized. Under the complete reference assumption, it has the same unit as TPM; under the incomplete reference assumption, the normalization is the defined as sum (flow * corresponding effective length) = number of reads.
ub          | Upper bound of the range of optimal abundances, normalized in the same way as lb.

Note that bounds can have floating point number inaccuracies.


# Analysis module for 16 Human Body Map samples

## Data accession
The SRA accession of the 16 Human Body Map RNA-seq samples can be found in `data/Metadata_hubmap.txt`. The following command will download all required reference files, Human Body Map samples, and MCF10 samples, and perform non-identifiability aware quantification for both datasets.
```
./script/metarun.sh <output folder>
```

## reproduce the plots
Using the results generated by scripts.metarun.sh script, the following R script will read the range of optima for 16 Human Body Map samples and plot the figures that are used in the manuscript. The figures are located in the directory where metarun.sh is called. The plot file will be located at `<output folder>/DeviationfromBound.pdf`.
```
Rscript scripts/draw_figure_hubmap.R <output folder as specified in metarun.sh>
```


# Analysis module for 6 MCF10 cell line samples

## Data accession
The SRA accession of the 6 MCF10 cell line samples can be found in `data/Metadata_mcf10.txt`. The following command will download all required reference files, Human Body Map samples, and MCF10 samples, and perform non-identifiability aware quantification for both datasets.
```
./script/metarun.sh <output folder>
```

## DE analysis and reproduce the plots
Using the results generated by scripts.metarun.sh script, the following R script will read in the range of optima for MCF10 cell line samples, perform DE detection using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html), and plot the histogram and example range of optima with varying reference proportion parameter. The plot file will be located at `<output folder>/IValue.pdf`.
```
Rscript scripts/draw_figure_mcf10.R <output folder as specified in metarun.sh>
```
