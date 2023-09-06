#!/usr/bin/env python
'''
This script is to search a gff file (a genome annotation file) for nuclear-encoded, mitochondrial-targeting gens (Nmt genes)
sys allows the user to provide arguments from the command line
sys.argv is a list of the arguments.
sys.argv[0] is the name of the python script.
sys.argv[1] should be the name of the Human_MitoCarta3.0.csv file
sys.argv[2] should be the name of the gff file being read in
sys.argv[3] should be the name of the gff file being written out
'''

import sys

# Functions here

if __name__ == "__main__":
    with open(sys.argv[1], 'r') as mito_carta_file:
        with open(sys.argv[2], 'r') as gff_file_in:
            with open(sys.argv[3], 'w') as gff_file_out:
                # for each line in mito_carta_file, create search criteria
                    # for each line in gff_file_in, check search criteria created above
                        # if search criteria matches gff_file_in line, write this line to gff_file_out
