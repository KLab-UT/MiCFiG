module load blast


#Declare environment variables
cds_fastas_dir=$1
db=$2
output_dir=$3




# LOOP through cds fastas
for fasta in /uufs/chpc.utah.edu/common/home/u6052680/MiCFiG/fastas/cds_fastas/*fasta;
	output_file=basesname("$fasta")
	blastn 

