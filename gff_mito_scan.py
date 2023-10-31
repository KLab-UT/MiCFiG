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
    log = open("log.txt", 'w')
    # for each line in mito_carta_file, create search criteria
    mito_carta_file.readline()
    search_terms = []
    for line in mito_carta_file:
        #print(line)
        # split the line up, and only select the element that contains the gene name
        line = line.split(',')
        human_symbol = line[0]
        search_terms.append(human_symbol)
        synonyms = line[1].split('|')
        if synonyms[0] != '-':
            search_terms.extend(synonyms)
        new_terms = ['|' + word + '_HUMAN' for word in search_terms]

       # description = line[2].strip()
       #search_terms.append(description)
       #we no longer want to search for the description so I commented this out of the search terms list
        #^splits up all of the terms in the csv file and appends them to a list and then prints them
        #to a file, they are all separate lists and some contain nothing

        #.extend allows you to add the terms from one list to another without
        #adding the whole list

        # for each line in gff_file_in, search for gene names


    for gff_line in gff_file_in:
        for term in new_terms:
            if term in gff_line:
                # if search criteria matches gff_file_in line, write this line to gff_file_out
                if gff_line.split("\t")[2] == "gene":
                    print(term)
                    log.write("HIT:\n\tterm: " + term + "\n\tgff: " + gff_line)
                    print(gff_line.strip(), file=gff_file_out)
    mito_carta_file.close()
    gff_file_in.close()
    gff_file_out.close()
    log.close()
