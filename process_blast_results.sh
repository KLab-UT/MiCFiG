#script for processing blast results
# For each chromosome hit, grab the minimum subject start and maximum subject
# end
# Creates output file for grabbing annotations with range overlap from gff
#   Ouput = Query ID , Chromosome , Start , Stop

# Can start with the test files that regan made

#Results file
results_dir=$1 # sratch blast_results
output_file=$2 # scratch processed_ ... file...

# Parse the results files in results dir (blast_results in scratch)
