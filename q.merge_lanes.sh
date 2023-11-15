#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N


# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_nmt_variation_data/mapped_reads
home=~/Biol_4310/whiptail_nmt_variation
# run bash script for counting unmerged reads
cd $wd

# */ is a wildcard pattern that matches all directories in the current directory
for reference in */; do
    echo $reference
	bash $home/merge_lanes.sh -i $wd/$reference -g ${reference///} -o $wd/$reference
    bash $home/merge_unmerged_and_merged.sh -i $wd/$reference -g ${reference///} -o $wd/$reference
done
