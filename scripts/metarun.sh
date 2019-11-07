#!/bin/bash

script_folder=$0
script_folder=${script_folder%/*}
out_folder=$1

# download reference and build index
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz | gunzip > ${out_folder}/gencode.v26.annotation.gtf
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.p10.genome.fa.gz | gunzip > ${out_folder}/GRCh38.p10.genome.fa
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz | gunzip > ${out_folder}/gencode.v26.transcripts.fa
salmon index --gencode -t ${out_folder}/gencode.v26.transcripts.fa -i ${out_folder}/gencode.v26.full
awk '{if($3=="transcript") print substr($10,2,length($10)-3)"\t"substr($12,2,length($12)-3)}' ${out_folder}/gencode.v26.annotation.gtf > ${out_folder}/gencode.v26.Gene_Trans_Map.txt

# download data
${script_folder}/../data/getdata.sh ${out_folder}/Reads/


GTFfile="${out_folder}/gencode.v26.annotation.gtf"
GenomeFasta="${out_folder}/GRCh38.p10.genome.fa"
SalmonIndex="${out_folder}/gencode.v26.full"
OutPrefix="${out_folder}/HumanBodyMap/salmon_Full"

# estimate flow
while read -r line; do

	ID=${line}

	python ${script_folder}/../src/process.py ${out_folder}/Reads/${ID} ${SalmonIndex} ${OutPrefix}_${ID} ${OutPrefix}_${ID}/prefixgraph ${GTFfile} ${GenomeFasta}

	# get salmon edge flow
	python ${script_folder}/../src/GetFlow_NodeEdge.py 1 ${GTFfile} ${OutPrefix}_${ID}/prefixgraph/gs_graph_fragstart.txt ${OutPrefix}_${ID}/prefixgraph/gs ${OutPrefix}_${ID}/quant.sf ${OutPrefix}_${ID}/prefixgraph/salmon

	# bounding uncertainty of annotated transcripts
	python ${script_folder}/../src/BoundingTranscriptFlows.py ${OutPrefix}_${ID}/prefixgraph/gs ${OutPrefix}_${ID}/prefixgraph/gs_result_ipopt_round0.pkl ${OutPrefix}_${ID}/prefixgraph/gs_maxflow_bound.txt
done < ${script_folder}/../data/Metadata.txt

