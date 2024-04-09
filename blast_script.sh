### Test usage: bash blast_script.sh -d fastas/test_fastas/ -b /scratch/general/nfs1/utu_4310/whiptail_shared_data/references/blast_dbs/a_marmoratus/a_marmoratus_AspMarm2.0.fasta.db -o test_blast_script

#!/bin/bash

# Default values
output_dir="output"

# Function to display script usage
usage() {
    echo "Usage: $0 -d <cds_fastas_dir> -b <db> [-o <output_dir>]"
    echo "Options:"
    echo "  -d <cds_fastas_dir> : Path to the directory containing CDS FASTA files."
    echo "  -b <db>             : BLAST database name."
    echo "  -o <output_dir>     : Output directory (default: output)."
    exit 1
}

# Parse arguments
while getopts ":d:b:o:" opt; do
    case $opt in
        d) cds_fastas_dir=$OPTARG ;;
        b) db=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check if mandatory arguments are provided
if [ -z "$cds_fastas_dir" ] || [ -z "$db" ]; then
    echo "Mandatory arguments missing!"
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Load blast 
module load blast

# Loop through cds fastas
for fasta in "$cds_fastas_dir"/*fasta; do
    # Get the filename without extension
    filename=$(basename "$fasta" .fasta)

    # Define output file path
    output_file="$output_dir/$filename.csv"
    
    # Perform BLASTN and save output to the defined file
    blastn -query "$fasta" -db "$db" -outfmt "6 qseqid sseqid sstart send" -out "$output_file"
done

