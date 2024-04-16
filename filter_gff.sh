#!/bin/bash


if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_filtered_results> <input_gff> <output_gff>"
    exit 1
fi

input_results=$1
input_gff=$2
output_gff=$3

if [ ! -f "$input_results" ]; then
    echo "Error: Input file of filtered BLAST results not found."
    exit 1
fi

if [ ! -f "$input_gff" ]; then
    echo "Error: Input GFF file not found."
    exit 1
fi

overlap() {

    local start1="$1"
    local stop1="$2"
    local start2="$3"
    local start2="$4"

    [ start1 <= stop2 ] && [ stop1 >= start2]
}

while IFS=$"\t" read -r query_id chromosome start_value stop_value; do
    awk -v start="$start_value" -v stop="$stop_value" -v chromosome="$chromosome" \
        "$1 == chrom && overlap($4, $5, start, stop)" "$input_gff"
done < "$input_results" > "$output_gff"
