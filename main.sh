blast_db=/scratch/general/nfs1/utu_4310/whiptail_shared_data/references/blast_dbs/a_marmoratus/a_marmoratus_AspMarm2.0.fasta.db
outputDir=/scratch/general/nfs1/utu_4310/whiptail_shared_data/blast_results

cd /uufs/chpc.utah.edu/common/home/u6052680/MiCFiG
# Use blastn with cds_fastas to search for alignments in blast_db 
#bash blast_script.sh -d fastas/cds_fastas -b $blast_db -o $outputDir/blastn_output -t blastn


#Use tblastn with protein_fastas to search for alignments in blast_db
#bash blast_script.sh -d fastas/protein_fastas -b $blast_db -o $outputDir/tblastn_output -t tblastn

# Convert blastn results to bed files
#bash process_blast_results.sh /scratch/general/nfs1/utu_4310/whiptail_shared_data/blast_results/blastn_output/ /scratch/general/nfs1/utu_4310/MiCFiG_wd/blastn_bed_files/

# Convert tblast results to bed files
bash process_blast_results.sh /scratch/general/nfs1/utu_4310/whiptail_shared_data/blast_results/tblastn_output/ /scratch/general/nfs1/utu_4310/MiCFiG_wd/tblastn_bed_files/
#
# 
#
#
#
#

