#script for processing blast results
# For each chromosome hit, grab the minimum subject start and maximum subject
# end
# Creates output file for grabbing annotations with range overlap from gff
#   Ouput = Query ID / Chromosome / Start / Stop


#Results directory: /scratch/general/nfs1/utu_4310/whiptail_shared_data/blastn_results

# Results dir =
results_dir=./test_results_dir
output_file=$2


ls $results_dir

# Parse the results files in results dir
