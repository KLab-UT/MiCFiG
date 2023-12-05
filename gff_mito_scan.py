#!/usr/bin/env python
'''
This script is to search a gff file (a genome annotation file) for nuclear-encoded,
mitochondrial-targeting gens (Nmt genes)

Argparse allows the user to provide arguments from the command line with tags.

Parser is an argument object that contains the arguments that are required of the user.

To run the code, the user should input the name of the python script with the
following inputs:

The input with the tag -m should be the name of the Human_MitoCarta3.0.csv file

The input with the tag -g should be the name of the gff file being read in

If the user should want a specific name for the output gff file, the input with
the tag -o should be the name of the gff file being written out. There is a
default value for this of output.gff.

If the user wants a specific name for the output log file, the input with
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

#the names that the arguments will be referred to as throughout the script
mito = args.mito_carta_file
gff = args.gff_file_in
out = args.gff_file_out
log = args.log_file

def Create_Search_Terms(mito):
    '''
    This function takes an input of a file with a list of mitochondrial-encoded
    genes (CSV file).

    It outputs a list of all of the search terms in the format "|term_HUMAN" without
    the phrase description of the gene and blank inputs.

    '''
    mito.readline()
    search_terms = []
    for line in mito:
        line = line.split(',')
        #divides up each column
        human_symbol = line[0]
        search_terms.append(human_symbol)
        synonyms = line[1].split('|')
        #divides multiple inputs in column 2
        if synonyms[0] != '-':
            #ignores the blank inputs
            search_terms.extend(synonyms)
        terms_list = ['|' + word + '_HUMAN' for word in search_terms]
        #formats all of the terms correctly for better search accuracy
    return terms_list

def Find_Matches(terms_list, gff, out):
    '''
    This function takes an input of the terms_list (from Create_Search_Terms),
    the .gff file with the genes to be sorted through, the output file to write
    the sorted genes to.

    It outputs the output file with the sorted genes.

    '''
    for line in gff:
        for term in terms_list:
            if term in line:
                #checks if the term in terms_list is in the gff file line
                if line.split('\t')[2] == 'gene':
                    #if the term is there the line must be a gene not an exon or transcript
                    #this avoids duplications
                    print(line.strip(), file = out)
                    #writes to the output gff file
    return

def Create_Log_File(terms_list, gff, log):
    '''
    This function takes an input of the terms_list (from Create_Search_Terms),
    the .gff file with the genes to be sorted through, the log file to print
    the matched genes to.

    It outputs the log file where the keyword that was matched to the gene
    precedes the gene that was matched.

    '''
    for line in gff:
        for term in terms_list:
            if term in line:
                #checks if the term in terms_list is in the gff file line
                if line.split('\t')[2] == 'gene':
                    #if the term is there the line must be a gene not an exon or transcript
                    #this avoids duplications
                    #print('HIT:\n\tterm:' + term + '\n\tgff:' + line, file=log)
                    log.write('HIT:\n\tterm:' + term + '\n\tgff:' + line)
                    #writes to the log file
    return



def main():
    '''
    This function takes no inputs but compiles all of the functions and runs
    them to first create the list of terms to search through the gff file and
    then create the output gff file and the log file.

    '''
    with open(mito, "r") as mito_terms:
        with open(gff, "r") as gff_in:
            with open(out, "w") as gff_out:
                with open(log, "w") as log_out:
                    terms = Create_Search_Terms(mito_terms)
                    Find_Matches(terms, gff_in, gff_out)
                    Create_Log_File(terms, gff_in, log_out)
    return

if __name__ == "__main__":
    main()


