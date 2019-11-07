#!/bin/bash

script_folder=$0
folder=$1

while read -r line; do
	IFS=$'\t' read -ra x <<< "${line}"

	echo "fastq-dump --split-files --gzip ${x[4]}"
	fastq-dump --split-files --gzip -O ${folder} ${x[4]}
done < ${script_folder%/*}/Metadata_mcf10.txt
