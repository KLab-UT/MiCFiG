## ADJUST ENVIRONMENT VARIABLES AS NEEDED ##
blast_db=/scratch/general/nfs1/utu_4310/whiptail_shared_data/references/blast_dbs/a_marmoratus/a_marmoratus_AspMarm2.0.fasta.db
gff_path=/scratch/general/nfs1/utu_4310/whiptail_shared_data/references/a_marmoratus_AspMarm2.0_v1.gff
outputDir=/scratch/general/nfs1/utu_4310/MiCFiG_wd

cd /uufs/chpc.utah.edu/common/home/u6052680/MiCFiG

# Make individual fasta files for blast queries 
python3 split_fasta.py mitocarta_proteins.fasta fastas/proteins/


# Use tblastn with protein_fastas to search for alignments in blast_db
bash blast_script.sh -d fastas/proteins -b $blast_db -o $outputDir/blast_results/tblastn_output -t tblastn

# Convert tblastn results to bed files
bash process_blast_results.sh $outputDir/blast_results/tblastn_output/ $outputDir/tblastn_bed_files/ genes_ids.csv

# Filter gff with tblastn BED files
rm -rf /scratch/general/nfs1/utu_4310/MiCFiG_wd/tblastn_gff_filtered
python3 filter_gff.py $outputDir/tblastn_bed_files $gff_path $outputDir/tblastn_gff_filtered genes_ids.csv

#Process filtered gff files
python3 process_gffs.py genes_ids.csv $outputDir/tblastn_gff_filtered/ tblastn_log.csv tblastn_merged.gff $outputDir/tblastn_gff_filtered/not_found.txt


