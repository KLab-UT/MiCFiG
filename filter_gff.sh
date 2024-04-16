#!/bin/bash


if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_filtered_results> <input_gff> <output_gff>"
    exit 1
fi

input_bed=$1
input_gff=$2
output_gff=$3

if [ ! -f "$input_bed" ]; then
    echo "Error: Input BED file of filtered BLAST results not found."
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
    local stop2="$4"

    [ "$start1" -le "$stop2" ] && [ "$stop1" -ge "$start2" ]
}

while IFS=$"\t" read -r chromosome start_value stop_value; do
    awk -v start="$start_value" -v stop="$stop_value" -v chromosome="$chromosome" \
        {if ($1 == chromosome && overlap($4, $5, start, stop)) print} "$input_gff"
done < "$input_bed" > "$output_gff"

# I'm assuming that the BED file is simply chromosome number, start coordinate,
# and stop coordinate.
