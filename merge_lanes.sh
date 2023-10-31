#!/bin/bash

{
usage="$(basename "$0") [-h] [-i <working_directory>] [-g reference genome directory] [-o <output_directory>] This program merges bam files that were sequenced on lanes and change the names of bam files sequences on single lanes
    -h  show this help text
    -i  Path to the working directory (the main directory for the raw reads)    
    -o  Path to the output directory (the main directory for the clean reads)   
    -g is the name of the reference directory"
options=':h:i:g:o:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    i) i=$OPTARG;;
    g) g=$OPTARG;;
    o) o=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done


echo ""
echo "Working directory: " $i
echo "Output directory: " $o
echo "Reference directory: " $g
echo ""


# mandatory arguments
if [ ! "$i" ] || [ ! "$g" ] || [ ! "$o" ]; then
  echo "arguments -o, and -i  must be provided"
  echo "$usage" >&2; exit 1
fi

#module load bwa/2020_03_19
module load samtools/1.16

### files look like this...
# RLK034_USD16091387L_HJNHCDSXX_L2_unpaired_2.bam 

cd $i
#merge (paired.bam) files sequenced across multiple lanes
# I use "[ -f <file_name> ] && do thing" to check if file exists before doing tasks
# merge/combine paired_1 and paired_2 and merge unpaired_1 and unpaired_2
# KLC098
[[ ! -f "KLC098_${g}_L4_paired.bam" ]] && samtools merge KLC098_${g}_L4_paired.bam KLC098_*_L4_paired_1.bam KLC098_*_L4_paired_2.bam
[[ ! -f "KLC098_${g}_L1_paired.bam" ]] && samtools merge KLC098_${g}_L1_paired.bam KLC098_*_L1_paired_1.bam KLC098_*_L1_paired_2.bam

[[ ! -f "KLC098_${g}_L4_unpaired.bam" ]] && samtools merge KLC098_${g}_L4_unpaired.bam KLC098_*_L4_unpaired.bam KLC098_*_L4_unpaired.bam
[[ ! -f "KLC098_${g}_L1_unpaired.bam" ]] && samtools merge KLC098_${g}_L1_unpaired.bam KLC098_*_L1_unpaired.bam KLC098_*_L1_unpaired.bam

# RLK019
[[ ! -f "RLK019_${g}_L3_paired.bam" ]] && samtools merge RLK019_${g}_L3_paired.bam RLK019_*_L3_paired_1.bam RLK019_*_L1_paired_2.bam
[[ ! -f "RLK019_${g}_L1_paired.bam" ]] && samtools merge RLK019_${g}_L1_paired.bam RLK019_*_L1_paired_1.bam RLK019_*_L1_paired_2.bam


[[ ! -f "RLK019_${g}_L3_unpaired.bam" ]] && samtools merge RLK019_${g}_L3_unpaired.bam RLK019_*_L3_unpaired.bam RLK019_*_L1_unpaired.bam
[[ ! -f "RLK019_${g}_L1_unpaired.bam" ]] && samtools merge RLK019_${g}_L1_unpaired.bam RLK019_*_L1_unpaired.bam RLK019_*_L1_unpaired.bam

# RLC0004
[[ ! -f "RLK004_${g}_paired.bam" ]] && samtools merge RLK004_${g}_paired.bam RLK004_*_paired_1.bam RLK004_*_paired_2.bam

[[ ! -f "RLK004_${g}_unpaired.bam" ]] && samtools merge RLK004_${g}_unpaired.bam RLK004_*_unpaired.bam RLK004_*_unpaired.bam

# RLK034
[[ ! -f "RLK034_${g}_paired.bam" ]] && samtools merge RLK034_${g}_paired.bam RLK034_*_paired_1.bam RLK034_*_paired_2.bam

[[ ! -f "RLK034_${g}_unpaired.bam" ]] && samtools merge RLK034_${g}_unpaired.bam RLK034_*_unpaired.bam RLK034_*_unpaired.bam

# merge/combine lanes for KLC098 and RLK019
[[ ! -f "KLC098_${g}_paired.bam" ]] && samtools merge KLC098_${g}_paired.bam KLC098_${g}_L4_paired.bam KLC098_${g}_L1_paired.bam

[[ ! -f "RLK019_${g}_paired.bam" ]] && samtools merge RLK019_${g}_paired.bam RLK019_${g}_L3_paired.bam RLK019_${g}_L1_paired.bam

[[ ! -f "KLC098_${g}_unpaired.bam" ]] && samtools merge KLC098_${g}_unpaired.bam KLC098_${g}_L4_unpaired.bam KLC098_${g}_L1_unpaired.bam

[[ ! -f "RLK019_${g}_unpaired.bam" ]] && samtools merge RLK019_${g}_unpaired.bam RLK019_${g}_L3_unpaired.bam RLK019_${g}_L1_unpaired.bam

#files look like this...
# RLK034_USD16091387L_HJNHCDSXX_L2_unpaired_2.bam 

#module unload bwa/2020_03_19
module unload samtools/1.16
}
