#!/bin/bash

folder=$0
folder=${folder%/*}

# download reference and build index
#wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz | gunzip > gencode.v26.annotation.gtf
#wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.p10.genome.fa.gz | gunzip > GRCh38.p10.genome.fa
#wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz | gunzip > gencode.v26.transcripts.fa
#salmon index --gencode -t gencode.v26.transcripts.fa -i gencode.v26.full
#awk '{if($3=="transcript") print substr($10,2,length($10)-3)"\t"substr($12,2,length($12)-3)}' gencode.v26.annotation.gtf > gencode.v26.Gene_Trans_Map.txt

# download data
#${folder}/../data/getdata.sh


GTFfile="gencode.v26.annotation.gtf"
GenomeFasta="GRCh38.p10.genome.fa"
OutPrefix="GEUVADIS/salmon_Full"

# estimate flow
while read -r line; do
	read -ra x <<< ${line}
	ID=${x[${#x[@]}-1]}

	python ${folder}/../src/process.py ${ID}/${ID} gencode.v26.full ${OutPrefix}_${ID} ${OutPrefix}_${ID}/prefixgraph ${GTFfile} ${GenomeFasta}

	# get salmon edge flow
	python ${folder}/../src/GetFlow_NodeEdge.py 1 ${GTFfile} ${OutPrefix}_${ID}/prefixgraph/gs_graph_fragstart.txt ${OutPrefix}_${ID}/prefixgraph/gs ${OutPrefix}_${ID}/quant.sf ${OutPrefix}_${ID}/prefixgraph/salmon

	# bounding uncertainty of annotated transcripts
	python ${folder}/../src/BoundingTranscriptFlows.py ${OutPrefix}_${ID}/prefixgraph/gs ${OutPrefix}_${ID}/prefixgraph/gs_result_ipopt_round0.pkl ${OutPrefix}_${ID}/prefixgraph/gs_bound_round0.txt
done < ${folder}/../data/Metadata.txt

# calculate weighted Jaccard similarity and output examples
python ${folder}/../src/GEUVADIS_bound_example.py
