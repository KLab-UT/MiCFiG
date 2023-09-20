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
    mito_carta_file = open(sys.argv[1], 'r')
    gff_file_in = open(sys.argv[2], 'r')
    gff_file_out = open(sys.argv[3], 'w')
    # for each line in mito_carta_file, create search criteria
    mito_carta_file.readline()
    for line in mito_carta_file:
        print(line)
        # split the line up, and only select the element that contains the gene name
        line = line.split(',')
        human_symbol = line[0]
        synonyms = line[1].split('|')
        description = line[2].strip()
        search_terms = synonyms
        search_terms = search_terms.append(human_symbol)
        print(search_terms)
        search_terms = search_terms.append(description)
        # for each line in gff_file_in, search for gene names
        lines_with_search_terms = [gff_line for gff_line in gff_file_in if any(term in gff_line for term in search_terms)]
        print(lines_with_search_terms)
        # if search criteria matches gff_file_in line, write this line to gff_file_out
    mito_carta_file.close()
    gff_file_in.close()
    gff_file_out.close()
