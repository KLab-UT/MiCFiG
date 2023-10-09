#! bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input files] [-g reference genome file] [-o output file where mapped reads go] [-t number of threds]
This program will map reads onto the reference directory
	-h show help text
	-i directory name where input files are located
	-g path and name of reference file
	-o output file where mapped reads go
	-t number of threads (WARNING: bwa is set to use 4 threads per sample)"
options=':h:i:g:o:t:'
while getopts $options option; do
	case "$option" in
		h) echo "$usage"; exit;;
		i) i=$OPTARG;;
		g) g=$OPTARG;;
		o) o=$OPTARG;;
		t) t=$OPTARG;;
		:) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
		\?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
	esac
done

echo ""
echo "working directory: " $i
echo "reference genome: " $g
echo "output file: " $o
echo "number of threads: " $t

# mandatory arguements

#if [ ! "$i" ] || [ ! "$g" ] || [ ! "$o" ] || [ ! "$t" ]; then
#	echo "arguments -i, -g, -o, and -t  must be provided"
#	echo "$usage" >&2; exit 1
#fi

# map reads --------------------------------------------------------------------------------

######################
# Indexing Reference #
######################


###########################################################################################
# Aligning datasets againts reference with bwa mem #
###########################################################################################

module load bwa/2020_03_19
module load samtools/1.16

echo "Indexing Reference: ${g}"
echo ""
#cd /scratch/general/nfs1/utu_4310/whiptail_shared_data/references
bwa index ${g}

# Create function that runs bwa and converts sam to bam
# Include below line in fastqToBam
# -t threads equal to number of reads. when I use ${t} in samtools it gives the following error:
# sort: option requires an argument -- '@'
#so instead of ${t} i put 2
echo "Beginning mapping"
echo ""
fastqToBam() {
    [[ ! -f "${3}/${1}_paired.sam" ]] && bwa mem -t 2 "$2" ${1}_paired_1.fq.gz ${1}_paired_2.fq.gz > ${3}/${1}_paired.sam
    echo "${3}/${1}_paired.sam completed."
    samtools sort ${3}/${1}_paired.sam > ${3}/${1}_paired_sorted.bam -@ 2
    echo "${3}/${1}_paired_sorted.bam completed."
}
export -f fastqToBam

echo "Aligning reads with reference with bwa mem."
cd $i
ls *fq.gz | cut -d "." -f "1,2" | parallel fastqToBam {} $g $o
	
module unload bwa/2020_03_19
module unload samtools/1.16
echo "Mapping complete"
}

