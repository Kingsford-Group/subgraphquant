#!/bin/bash

# download reference and build index
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz | gunzip > gencode.v26.annotation.gtf
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.p10.genome.fa.gz | gunzip > GRCh38.p10.genome.fa
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz | gunzip > gencode.v26.transcripts.fa
salmon index --gencode -t gencode.v26.transcripts.fa -i gencode.v26.full

# download data
./../data/getdata.sh


GTFfile="gencode.v26.annotation.gtf"
GenomeFasta="GRCh38.p10.genome.fa"
OutPrefix="GEUVADIS/salmon_Full"

# estimate flow
while read -r line; do
	read -ra x <<< ${line}
	ID=${x[${#x[@]}-1]}

	python ../src/process.py ${ID}/${ID} gencode.v26.full ${OutPrefix}_${ID} ${OutPrefix}_${ID}/prefixgraph ${GTFfile} ${GenomeFasta}

	# get salmon edge flow
	python ../src/GetFlow_NodeEdge.py 1 ${GTFfile} ${OutPrefix}_${ID}/prefixgraph/gs_graph_fragstart.txt ${OutPrefix}_${ID}/prefixgraph/gs ${OutPrefix}_${ID}/quant.sf ${OutPrefix}_${ID}/prefixgraph/salmon

	# bounding uncertainty of annotated transcripts
	python ../src/BoundingTranscriptFlows.py ${OutPrefix}_${ID}/prefixgraph/gs ${OutPrefix}_${ID}/prefixgraph/gs_opt_ipopt_round0.pkl ${OutPrefix}_${ID}/prefixgraph/gs_bound_round0.txt
done < ../data/Metadata.txt
