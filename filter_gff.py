#!/usr/bin/env python3
#
# Script to idnetify overlaps betwenn BED and GFF files for each BED file in the provided directory
#
# 

import sys
import os

class GFF_entry:
    def __init__(self, gff_line):
        gff_fields = gff_line.strip().split('\t')
        self.g_chrom = gff_fields[0]
        self.g_start = int(gff_fields[3])
        self.g_stop = int(gff_fields[4])
        self.g_line = gff_line

class BED_entry:
    def __init__(self, bed_line):
        b_fields = bed_line.strip().split('\t')
        try:
            self.b_chrom = b_fields[0]
            self.b_start = int(b_fields[1])
            self.b_stop = int(b_fields[2])
            self.b_line = bed_line
        except:
            print("error making BED_entry for: {bed_line}")

def overlap(bed_entry, gff_entry):
        # Check for overlap in the start stop values of GFF and BED entry
        if bed_entry.b_start <= gff_entry.g_stop and bed_entry.b_stop >= gff_entry.g_start:
            print("Overlap found!")
            return True  # There's an overlap
        else:
            return False  # There's no overlap
        
def make_gff_dict(input_gff):
    # Read a GFF file into a dictionary organized by chromosome, return the dictionary
    gff_dict_by_chrom = {}

    # Skip comment lines at the beginning of the GFF file
    with open(input_gff, 'r') as gff_file:
        for line in gff_file:
            if not line.startswith('##'):
                break

    # Process non-comment lines
    with open(input_gff, 'r') as gff_file:
        for line in gff_file:
            if not line.startswith('##'):
                # Initilize GFF_entry object
                g_entry = GFF_entry(line)

                # If the chromosome key doesn't exist, initialize it with an empty array
                if g_entry.g_chrom not in gff_dict_by_chrom:
                    gff_dict_by_chrom[g_entry.g_chrom] = []

                # Append the GFF entry to the list for the corresponding chromosome
                gff_dict_by_chrom[g_entry.g_chrom].append(g_entry)

    # Print all the chromosomes
    print("Chromosomes found in the GFF file:")
    for chrom in gff_dict_by_chrom:
        print(chrom)
    return gff_dict_by_chrom


def main():
    if len(sys.argv) != 4:
        print("Usage: {} <input_bed_directory> <input_gff> <output_gff_directory>".format(sys.argv[0]))
        sys.exit(1)

    input_bed_dir = sys.argv[1]
    input_gff = sys.argv[2]
    output_gff_dir = sys.argv[3]

    # Create output directory if it doesn't exist
    os.makedirs(output_gff_dir, exist_ok=True)

    if not os.path.isdir(input_bed_dir):
        print("Error: Input BED directory from blast results not found.")
        sys.exit(1)

    if not os.path.isfile(input_gff):
        print("Error: Input GFF file not found.")
        sys.exit(1)

    if not os.path.isdir(output_gff_dir):
        print("Error: Output GFF directory not found.")
        sys.exit(1)

    print(input_bed_dir)
    print(input_gff)
    print(output_gff_dir)

    gff_dict = make_gff_dict(input_gff)

    # loop through bed files and check for overlaps in gff
    for input_bed in [f for f in os.listdir(input_bed_dir) if f.endswith('.bed')]:
        input_bed_path = os.path.join(input_bed_dir, input_bed)
        print(input_bed_path)
        no_overlap_output = os.path.join(output_gff_dir, "no_overlap.bed") # Output for sequences hit by blast but with no overlap in gff
        not_found_output = os.path.join(output_gff_dir, "not_found.txt") # Output for sequences not hit by blast
        if os.path.getsize(input_bed_path) > 0:  # Check if the BED file is not empty
            print("checking for overlap in", input_bed_path)
            output_gff = os.path.join(output_gff_dir, os.path.basename(input_bed_path) + ".gff")

            print("output_gff:", output_gff)
            with open(input_bed_path, 'r') as bed_file:
                for b_line in bed_file:
                    b_entry = BED_entry(b_line)

                    # Check if the chromosome key exists in the dictionary
                    if b_entry.b_chrom in gff_dict:
                        # Loop through GFF entries for the same chromosome
                        for g_entry in gff_dict[b_entry.b_chrom]:
                            if overlap(b_entry, g_entry):
                                with open(output_gff, 'a') as out_gff_file:
                                    print(f"Writing overlapping gff line to {output_gff}\n{g_entry.g_line}")
                                    out_gff_file.write(g_entry.g_line + '\n')
                    else:
                        print(f"No overlap found in GFF, writing BED entry to: {no_overlap_output}")
                        with open(no_overlap_output, 'a') as no_overlap:
                            no_overlap.write(b_line) 
        else:
            print(f"No overlap, BED file is empty, writing seq_id to: {not_found_output}")
            with open(not_found_output, 'a') as not_found:
                seq_id = input_bed.replace(".bed", "")
                not_found.write(seq_id + '\n')

if __name__ == "__main__":
    main()
