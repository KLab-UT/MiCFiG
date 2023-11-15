#!/bin/bash

{
usage="$(basename "$0") [-h] [-i <working_directory>] [-g <mapped reads directory>] [-o <output_directory>] This program merges bam files that were sequenced on lanes and change the names of bam files sequences on single lanes
    -h  show this help text
    -i  Path to the working directory (the main directory for the raw reads)

    -o  Path to the output directory (the main directory for the clean reads)
    -g is the name of the mapped reads directory"
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
echo "Mapped reads directory: " $g
echo ""


# mandatory arguments
if [ ! "$i" ] || [ ! "$g" ] || [ ! "$o" ]; then
  echo "arguments -o, and -i  must be provided"
  echo "$usage" >&2; exit 1
fi

module load samtools/1.16

cd $i
#merge (merged.bam) files and (unmerged.bam) files
# KLC098, RLK019, RLK004, RLK034
# "[ -f <file_name> ] && do thing" checks if file exists before doing tasks
[[ ! -f "KLC098_${g}.bam" ]] && samtools merge KLC098_${g}.bam KLC098_*_paired.bam KLC098_*_unpaired.bam

[[ ! -f "RLK019_${g}.bam" ]] && samtools merge RLK019_${g}.bam RLK019_*_paired.bam RLK019_*_unpaired.bam

[[ ! -f "RLK004_${g}.bam" ]] && samtools merge RLK004_${g}.bam RLK004_*_paired.bam RLK004_*_unpaired.bam

[[ ! -f "RLK034_${g}.bam" ]] && samtools merge RLK034_${g}.bam KLC098_*_paired.bam RLK034_*_unpaired.bam

#files look like this...
#RLK004_am_unmerged.bam
#RLK034_am_merged.bam

module unload samtools/1.16
}
