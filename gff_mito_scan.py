#!/usr/bin/env python
'''
This script is to search a gff file (a genome annotation file) for nuclear-encoded,
mitochondrial-targeting gens (Nmt genes)

Argparse allows the user to provide arguments from the command line with tags.

Parser is an argument object that contains the arguments that are required of the user.

To run the code the user should input the name of the python script with the
following inputs:

The input with the tag -m should be the name of the Human_MitoCarta3.0.csv file

The input with the tag -g should be the name of the gff file being read in

If the user should want a specific name for the output gff file, the input with
the tag -o should be the name of the gff file being written out. There is a
default value for this of output.gff.

If the user should want a specific name for the output log file, the input with
the tag -l should be the name of the log file being written out. There is a
default value for this of log.txt

'''

import argparse

parser = argparse.ArgumentParser()
#creates the arguments to be passed through the functions
parser.add_argument("--mito_carta_file", '-m', help="Mitochondrial DNA CSV file")
parser.add_argument("--gff_file_in", "-g", help="Input gff file")
parser.add_argument("--gff_file_out", "-o", default="output.gff", help="Output gff file")
parser.add_argument("--log_file", "-l", default="log.txt", help="Log to keep track of the matches in the genes")

args = parser.parse_args()

#the names that the arguments will be refered to as throughout the script
mito = args.mito_carta_file
gff = args.gff_file_in
out = args.gff_file_out
log = args.log_file

def Create_Search_Terms(mito):
    '''
    This function takes an input of a file with a list of mitochondrial encoded
    genes (CSV file).
    It outputs a list of all of the search terms in the format "|term_HUMAN" without
    the phrase description of the gene and blank inputs.

    '''
    terms_file.readline()
    search_terms = []
    for line in mito:
        line = line.split(',')
        #divides up each column
        human_symbol = line[0]
        search_terms.append(human_symol)
        synonyms = line[1].split('|')
        #divides multiple inputs in column 2
        if synonyms[0] != '-':
            #ignores the blank inputs
            search_terms.extend(synonyms)
        terms_list = ['|' + word + '_HUMAN' for word in search_terms]
        #formats all of the terms correctly for better search acuracy
    return terms_list

def Find_Matches(terms_list, gff, out, log):
    '''
    This function takes an input of the terms_list (from Create_Search_Terms),
    the .gff file with the genes to be sorted through, the output file to write
    the sorted genes to, and the log file to write the log of the matched term with
    the gene to.

    It outputs the output file with the sorted genes and the log file detailing
    the matches corresponding to the gene found.

    '''
    for line in gff:
        for term in terms_list:
            if term in line:
                #checks if the term in terms_list is in the gff file line
                if line.split('\t')[2] == 'gene':
                    #if the term is there the line must be a gene not an exon or transcript
                    #this avoids duplicaions
                    log.write('HIT:\n\tterm: ' + term + '\n\tgff: ' + line)
                    #writes to the log file
                    print(line.strip(), file = out)
                    #writes to the output gff file
    return

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
