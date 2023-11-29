#! /bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input bam files] [-g reference genome] [-o ouput file directory]
This program will get Variants and make a VCF file
        -h show help text
        -i directory name where input bam files are located
        -g reference genome
		-o output file directory"

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

module load bcftools/1.16
module load samtools/1.16

##adding -I skips Indels
#bcftools view -I -bvcg - > variants.bcf

#samtools mpileup -uf reference.fasta output.bam | bcftools view -bvcg - > variants.bcf
#bcftools view variants.bcf | vcfutils.pl varFilter > filtered_variants.vcf

get_variants() {
	name=$1
	reference=$2
	output=$3
    echo "Name: ${name}\nReference: ${reference}\nOutput: ${output}"
# -Ob = output type is combressed or binary
# -o is output
# -f fasta reference file
	# bcftools mpileup -Ob -o <study.bcf> -f <ref.fa> <sample1.bam>
	bcftool mpileup -Ob -o ${output}/${name}.bcf > -f ${reference} ${name}_combined.bam
# -v = vatriants only
# -m = mark sites
# -Oz = compressed VCF
# -Ov = uncompressed VCF
	# bcftools call -vmO z -o <study.vcf.gz> <study.bcf>
	bcftools call -vmOv -o ${name}.vcf.gz ${name}.bcf
}
export -f get_variants

echo "Getting Variants."
cd $i
ls *_combined.bam | cut -d "." -f "1" | parallel fastqToBam {} $g $o
# bam file example KLC098_USD16091388L_HKFJFDSXX_L4_paired_1.bam
# bam file example KLC098_USD16091388L_HKFJFDSXX.bam

module unload bcftools/1.16
module unload samtools/1.16
}
