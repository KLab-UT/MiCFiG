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

    local blast_start="$1"
    local blast_stop="$2"
    local gff_start="$3"
    local gff_stop="$4"

    [ "$blast_start" -le "$gff_stop" ] && [ "$blast_stop" -ge "$gff_start" ]
}

<<<<<<< HEAD
while IFS=$"\t" read -r b_chrom b_start b_stop; do
    while IFS=$"\t" read -r g_chrom g_source g_feature g_start g_stop g_rest; do
        if [[ "$b_chrom" == "$g_chrom" && $(overlap "$b_start" "$b_stop" "$g_start" "$g_stop") ]]; then
            echo -e "$g_chrom\t$g_source\t$g_featutre\t$g_start\t$g_stop\t$g_rest" >> "$3"
        fi
    done < "$2"
done < "$1"

#while IFS=$"\t" read -r chromosome start_value stop_value; do
#    echo($chromsome)
#    echo($start_value)
#    echo($stop_value)
=======
cat $input_bed

awk '{overlap($1, $2, $3)}' $input_bed

#awk -v blast_start="$start_value" -v blast_stop="$stop_value" -v blast_chromo="$chromosome" \
#    {if ($1 == blast_chromo && overlap(blast_start, blast_stop, $4, $5)) print} "$input_gff"
#
#while IFS=$"\t" read -r chromosome start_value stop_value; do
#    echo $chromsome
#    echo $start_value
#    echo $stop_value
>>>>>>> 4af798a6d0a086475a2a22cffa7675d1388b5144
##    awk -v blast_start="$start_value" -v blast_stop="$stop_value" -v blast_chromo="$chromosome" \
##        {if ($1 == blast_chromo && overlap(blast_start, blast_stop, $4, $5)) print} "$input_gff"
#done < "$input_bed" > "$output_gff"
#
# I'm assuming that the BED file is simply chromosome number, start coordinate,
# and stop coordinate.
