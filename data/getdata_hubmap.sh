#!/bin/bash

script_folder=$0
folder=$1

while read -r line; do

	if [[ ! -e ${folder}/${line}_1.fastq.gz ]]; then
		wget -P ${folder} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${line:0:6}/${line}/${line}_1.fastq.gz
	fi

	if [[ ! -e ${folder}/${line}_2.fastq.gz ]]; then
		wget -P ${folder} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${line:0:6}/${line}/${line}_2.fastq.gz
	fi
done < ${script_folder%/*}/Metadata_hubmap.txt
