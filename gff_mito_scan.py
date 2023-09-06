#!/usr/bin/env python

import sys

# Functions here

if __name__ == "__main__":
    with open(sys.argv[1], 'r') as mito_carta_file:
        with open(sys.argv[2], 'r') as gff_file_in:
            with open(sys.argv[3], 'w') as gff_file_out:
                # for each line in mito_carta_file, create search criteria
                    # for each line in gff_file_in, check search criteria created above
                        # if search criteria matches gff_file_in line, write this line to gff_file_out
