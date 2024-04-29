import argparse
import os
import csv

def get_query_gene(query_id, gene_ids_csv):
    """
    Uses the query id to find gene name in csv with genes and ids. Gene needs to be first field in csv w/ ids
    """
    with open(gene_ids_csv, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if query_id in row:
                return row[0]  # Return the first field of the row where the ID is found
    return None  # Return None if the query ID is not found in the CSV

def main(gene_ids_csv, gff_directory, output_csv, merged_gff, not_found_file):
    # List all files in the directory
    files = os.listdir(gff_directory)
    
    # Filter out only the .gff files
    gff_files = [file for file in files if file.endswith('.gff')]
    
    # Create a set to keep track of processed GFF files
    processed_files = set()

    # get ids for query not hit by blast
    with open(not_found_file, 'r') as nf:
        not_found_list = nf.read().splitlines()

    # Open the output CSV file for writing
    with open(output_csv, 'w', newline='') as csvfile:
        # Define column headers
        fieldnames = ['ID', 'Gene', 'Annotation']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Write column headers
        writer.writeheader()
        
        # Open the merged GFF file for writing
        with open(merged_gff, 'w') as merged_file:
            # Loop through each .gff file
            for gff_file in gff_files:
                query_id = gff_file.split(".")[0]
                human_gene = get_query_gene(query_id, gene_ids_csv)
                gff_path = os.path.join(gff_directory, gff_file)
                
                # Check if the current GFF file has been processed
                if gff_file not in processed_files:
                    # Log message for each new GFF file
                    print(f"processing: {gff_path}")
                    processed_files.add(gff_file)  # Add the current GFF file to the set of processed files

                    # Read the first line of the GFF file and use it as annotation
                    with open(gff_path, 'r') as f:
                        gff_annotation = f.readline().strip()
                        if gff_annotation == "":
                            gff_annotation = "no overlap"

                # get annotations
                if query_id in not_found_list:
                    gff_annotation = "not found by blast"
                # Write overlaping annotations to a merged gff
                if gff_annotation != "not found by blast" and gff_annotation != "no overlap":
                    merged_file.write(gff_annotation)
                    merged_file.write("\n")
                
                # Write id, gene, and annotation to output csv with column headers
                writer.writerow({'ID': query_id, 'Gene': human_gene, 'Annotation': gff_annotation})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process .gff files and gene IDs CSV.')
    parser.add_argument('gene_ids_csv', type=str, help='Path to the CSV file containing gene IDs.')
    parser.add_argument('gff_directory', type=str, help='Path to the directory containing .gff files.')
    parser.add_argument('output_csv', type=str, help='Path to the output CSV file.')
    parser.add_argument('merged_gff', type=str, help='Path to the merged GFF file.')
    parser.add_argument('not_found', type=str, help='Path to not_found.txt')
    args = parser.parse_args()

    main(args.gene_ids_csv, args.gff_directory, args.output_csv, args.merged_gff, args.not_found)
