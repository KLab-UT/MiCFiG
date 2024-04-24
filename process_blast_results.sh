#script for processing blast results
# For each chromosome hit, grab the minimum subject start and maximum subject
# end


module load blast

results_dir=$1 # Directory to blast results in XML format
output_dir=$2 # Directory to save BED outputs
gene_ids=$3 # CSV file with genes and associated ids

cd /uufs/chpc.utah.edu/common/home/u6052680/MiCFiG

# Check if the input directories exist
if [ ! -d "$results_dir" ]; then
    echo "Error: Directory '$results_dir' does not exist."
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Parse the results files in the results directory
for result_file in "$results_dir"*.xml; do
    if [ -f "$result_file" ]; then
        result_filename=$(basename "$result_file")
        bed_out="$output_dir${result_filename%.xml}.bed"

        python3 process_blast.py "$result_file" "$gene_ids" --bed_output_path "$bed_out"
    fi
done
