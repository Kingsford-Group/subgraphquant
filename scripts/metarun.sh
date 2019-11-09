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

#################### Human Body Map dataset ####################
# download data
${script_folder}/../data/getdata.sh ${out_folder}/Reads_hubmap/

GTFfile="${out_folder}/gencode.v26.annotation.gtf"
GenomeFasta="${out_folder}/GRCh38.p10.genome.fa"
SalmonIndex="${out_folder}/gencode.v26.full"
OutPrefix="${out_folder}/HumanBodyMap"

# estimate flow for Human Body Map dataset
while read -r line; do

	ID=${line}

	python ${script_folder}/../src/process.py ${out_folder}/Reads_hubmap/${ID} ${SalmonIndex} ${OutPrefix}/${ID} ${OutPrefix}/${ID}/prefixgraph ${GTFfile} ${GenomeFasta}

done < ${script_folder}/../data/Metadata_hubmap.txt


#################### MCF10 dataset ####################
# download data
${script_folder}/../data/getdata.sh ${out_folder}/Reads_mcf10/

GTFfile="${out_folder}/gencode.v26.annotation.gtf"
GenomeFasta="${out_folder}/GRCh38.p10.genome.fa"
SalmonIndex="${out_folder}/gencode.v26.full"
OutPrefix="${out_folder}/MCF10"

# estimate flow for Human Body Map dataset
while read -r line; do

	IFS=$'\t' read -a x <<< ${line}
	id=${x[4]}
	echo ${id}

	python ${script_folder}/../src/process.py ${out_folder}/Reads_mcf10/${id} ${SalmonIndex} ${OutPrefix}/${id}/ ${OutPrefix}/${id}/prefixgraph ${GTFfile} ${GenomeFasta}

done < ${script_folder}/../data/Metadata_mcf10.txt

# calculate I value
python ${script_folder}/ComputeIValue.py ${out_folder}
