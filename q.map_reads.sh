#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=24
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N

# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_nmt_variation_data/
mkdir -p $wd/mapped_reads/Aspidoscelis_arizonae_fasta
mkdir -p $wd/mapped_reads/Aspidoscelis_marmoratus_fasta
#mkdir -p $wd/mapped_reads/Aspidoscelis_arizonae_gff
#mkdir -p $wd/mapped_reads/Aspidoscelis_marmoratus_gff

echo "Beggining mapping"
MapReads() {
	wd=/scratch/general/nfs1/utu_4310/whiptail_nmt_variation_data/
	echo "############################"
	echo -e "Species: ${1}\nGenome: ${2}"
	# you're passing in n threads where n is number of reads (6) multiplied by number of threads used by functions in map_reads.sh (2).
	# This is done for every genome you want to map (27 genomes listed in ref_genomes.txt).
	bash map_reads.sh -i $wd/trimmed_reads/paired_reads -g $wd/references/${2} -o $wd/mapped_reads/${1} -t 12
	echo "merched read mapped"
	bash map_reads_paired.sh -i $wd/trimmed_reads/unpaired_reads -g $wd/references/${2} -o $wd/mapped_reads/${1} -t 12
	echo "unmerged read mapped"
}
export -f MapReads

# use references in ref_genome in the MapReads function and run each line in parallel
# {1} is the directory name and {2} is the file name of the genome reference mathcing that directory
# grep -v filter out lines starting with #

#grep -v '^#' ref_genomes.txt | cut -d " " -f 1,2 | parallel MapReads {1} {2}
grep -v '^#' ref_genomes.txt | parallel --colsep ' ' MapReads {1} {2}

echo "mapping done"

# Keep in mind that if your MapReads function itself uses multiple threads (e.g., with the -t option set to 4 as in your example), the total number of threads used will be a combination of those used by MapReads and those used by parallel.

