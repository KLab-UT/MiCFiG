## ADJUST ENVIRONMENT VARIABLES AS NEEDED ##
blast_db=/scratch/general/nfs1/utu_4310/whiptail_shared_data/references/blast_dbs/a_marmoratus/a_marmoratus_AspMarm2.0.fasta.db
gff_path=/scratch/general/nfs1/utu_4310/whiptail_shared_data/references/a_marmoratus_AspMarm2.0_v1.gff
outputDir=/scratch/general/nfs1/utu_4310/MiCFiG_wd

script_dir="$(dirname "$(realpath "$0")")"

# Make individual fasta files for blast queries 
# python3 "$script_dir/split_fasta.py" "$script_dir/mitocarta_proteins.fasta" "$script_dir/fastas/proteins/"

# Use tblastn with protein_fastas to search for alignments in blast_db
# bash "$script_dir/blast_script.sh" -d "$script_dir/fastas/proteins" -b $blast_db -o $outputDir/blast_results/tblastn_output -t tblastn

# Convert tblastn results to bed files
# bash "$script_dir/process_blast_results.sh" $outputDir/blast_results/tblastn_output/ $outputDir/tblastn_bed_files/ "$script_dir/genes_ids.csv"

# Filter gff with tblastn BED files
#rm -rf /scratch/general/nfs1/utu_4310/MiCFiG_wd/tblastn_gff_filtered
#python3 "$script_dir/filter_gff.py" $outputDir/tblastn_bed_files $gff_path $outputDir/tblastn_gff_filtered "$script_dir/genes_ids.csv"

# Add gene names from MitoCarta database to individual GFF files of species

python3 "$script_dir/add_gene_name.sh" "$script_dir/mitocarta_proteins.fasta" "$ouput_dir/tblastn_gff_filtered"

#Process filtered gff files
python3 "$script_dir/process_gffs.py" "$script_dir/genes_ids.csv" $outputDir/tblastn_gff_filtered/ "$script_dir/tblastn_log.csv" "$script_dir/tblastn_merged.gff" $outputDir/tblastn_gff_filtered/not_found.txt


