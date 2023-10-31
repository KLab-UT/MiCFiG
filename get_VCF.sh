#! /bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input files] [-g mapped reads directory]
This program will map reads onto the reference directory
        -h show help text
        -i directory name where input files are located
        -g mapped reads directory"

options=':h:i:g:'
while getopts $options option; do
        case "$option" in
                h) echo "$usage"; exit;;
                i) i=$OPTARG;;
                g) g=$OPTARG;;
                :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
                \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
        esac
done

echo ""
echo "Working directory: ${i}"
echo "Mapped reads directory: ${g}"
echo "Output directory: ${o}"
echo ""

# mandatory arguments
if [ ! "$i" ] || [ ! "$g" ] || [ ! "$o" ]; then
  echo "arguments -i, -g, and -o  must be provided"
  echo "$usage" >&2; exit 1
fi

#==========================================================================================================
#Consensus/variant calling options:
#       -c  SNP calling (force -e)
#       -d FLOAT  skip loci where less than FLOAT fraction of samples covered [0]
#       -e        likelihood based analyses
#       -g        call genotypes at variant sites (force -c)
#       -i FLOAT  indel-to-substitution ratio [-1]
#       -I        skip indels
#       -m FLOAT  alternative model for multiallelic and rare-variant calling, include if P(chi^2)>=FLOAT
#       -p FLOAT  variant if P(ref|D)<FLOAT [0.5]
#       -t FLOAT  scaled substitution mutation rate [0.001]
#       -T STR    constrained calling; STR can be: pair, trioauto, trioxd and trioxs (see manual) [null]
#       -v        output potential variant sites only (force -c)
#==========================================================================================================

module load BCFtools/1.3.1
module load SAMtools/1.3.1

##adding -I skips Indels
#bcftools view -I -bvcg - > variants.bcf

#samtools mpileup -uf reference.fasta output.bam | bcftools view -bvcg - > variants.bcf
#bcftools view variants.bcf | vcfutils.pl varFilter > filtered_variants.vcf

# List all BAM files in current directory
bam_files=$(ls *.bam | tr '\n' ' ')

# Create a space-separated string of bam files
BAM_FILES=$(echo $bam_files | tr ' ' '\n' | paste -s -d ' ')

# Run mpileup on all bam files and pipe the output to BCFtools
samtools mpileup -g -B $BAM_FILES | bcftools view -Ou - > all_samples.bcf

module unload BCFtools/1.3.1
module unload SAMtools/1.3.1
