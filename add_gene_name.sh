#!/usr/bin/env bash

input_fasta="$1"
input_gff_dir="$2"
nfout="gff_not_found.txt"

python3 extract_gene_name.py -i "$input_fasta" -o "extract_output.txt"

> "$nfout"

while read -r file_id gene_name; do
	gff_file="${input_gff_dir}/${file_id}.gff"
	if [ -f "$gff_file" ]; then
		if grep -q "mito_carta_id=${gene_name}" "$gff_file"; then
  			echo "Already added."
			exit 0
                else	
			sed -i "s#blast_id=#mito_carta_id=${gene_name};blast_id=#" "$gff_file"
			echo "Added gene_name to $gff_file"
		fi
	else
		echo "$file_id" >> "$nfout"
		echo "$gene_name not found"
	fi
done < extract_output.txt

echo "Done updating MitoCarta names!"
