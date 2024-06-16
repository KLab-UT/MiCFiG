#!/usr/bin/env bash

input_fasta_dir="$1"
input_fasta="$2"
input_gff_dir="$3"
nfout="gff_not_found.txt"

cd "$input_fasta_dir"

python3 extract_gene_name.py -i "$input_fasta" -o "extract_output.txt"

> "$nfout"

cd "$input_gff_dir"

while read -r file_id gene_name; do
	gff_file="${input_gff_dir}/${file_id}.gff"
	if [ -f "$gff_file" ]; then
		sed -i "s#blast_id=#mito_carta_id=${gene_name};blast_id=#" "$gff_file"
		echo "Added gene_name to $gff_file"
	else
		echo "$gff_file does not exist." >> "$nfout"
	fi
done < extract_output.txt
