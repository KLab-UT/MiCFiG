#script for processing blast results
# For each chromosome hit, grab the minimum subject start and maximum subject
# end
# Creates output file for grabbing annotations with range overlap from gff
#   Ouput = Query ID / Chromosome / Start / Stop
results_dir=$1
output_file=$2

# Parse the results files in results dir
