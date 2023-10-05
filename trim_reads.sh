#! /bin/bash
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N

# PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...

# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_nmt_variation_data
mkdir -p $wd/trimmed_reads
mkdir -p $wd/trimmed_reads/paired_reads
mkdir -p $wd/trimmed_reads/unpaired_reads

module load trimmomatic

cd $wd/raw_reads
ls *1.fq.gz | cut -d "_" -f 1,2,3,4 | sort | uniq > fq_list.txt
while read sample; do
	read1=$wd/raw_reads/${sample}_1.fq.gz 
	read2=$wd/raw_reads/${sample}_2.fq.gz
	output_p_1=$wd/trimmed_reads/paired_reads/${sample}_paired_1.fq.gz
	output_u_1=$wd/trimmed_reads/unpaired_reads/${sample}_unpaired_1.fq.gz
	output_p_2=$wd/trimmed_reads/paired_reads/${sample}_paired_2.fq.gz
	output_u_2=$wd/trimmed_reads/unpaired_reads/${sample}_unpaired_2.fq.gz
	trimmomatic PE $read1 $read2 $output_p_1 $output_u_1 $output_p_2 $output_u_2 
	#ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True
	#trimmomatic PE -threads 8 -trimlog trimlog.txt -summary summary.txt -basein $read1 $read2 -baseout $output_paired_1 $output_unpaired_1 $output_paired_2 $output_unpaired_2

done<fq_list.txt

module unload trimmomatic
