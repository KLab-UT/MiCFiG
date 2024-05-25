#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=24
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N

# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_nmt_variation_data
mkdir -p $wd/VCF_files

echo "Beggining variant calling"
call_variants() {
	wd=/scratch/general/nfs1/utu_4310/whiptail_nmt_variation_data
    species=$1
    genome=$2
	echo "############################"
	echo -e "Species: ${species}\nGenome: ${genome}"
######
	bash get_VCF.sh -i $wd/mapped_reads/${species} -g $wd/references/${genome} -o $wd/VCF_files/${species}
	echo "${species} variant calling done"
}
export -f call_variants

# use references in ref_genome in the MapReads function and run each line in parallel
# {1} is the directory name and {2} is the file name of the genome reference mathcing that directory
# grep -v filter out lines starting with #

#grep -v '^#' ref_genomes.txt | cut -d " " -f 1,2 | parallel MapReads {1} {2}
grep -v '^#' ref_genomes.txt | parallel --colsep ' ' call_variants {1} {2}

echo "variant calling complete"

