#!/usr/bin/env/ bash

input_dir="$1"
input_fasta="$2"

python3 extract_gene_name.py -i "$input_fasta" -o "extract_output.txt"

while read file_id gene_name; do
	gff_file="${input_dir}/${file_id}.gff"
	if [ -f "$gff_file" ]; then
		sed -i "s#blast_id=#mito_carta_id=${gene_name};blast_id=#" "$gff_file"
		echo "Added gene_name to $gff_file"
	else
		echo "$gff_file does not exist."
	fi
done < extract_output.txt
