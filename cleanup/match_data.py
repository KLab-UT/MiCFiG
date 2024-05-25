import sys
import numpy as np
import matplotlib.pyplot as plt

'''
This script functions to compile data from the mito carta file and the matched
genes and display them in a graph

It uses sys.argv for passing through arguments

sys.argv[0] should be the name of the function
sys.argv[1] should be the name of the mitocarta file
sys.argv[2] shoul be the name of the log file created in gff_mito_scan.py

'''

def Hits_List(log_file):
        '''
        creates a list of all of the terms that hit a gene in the gff file

        input: log file (txt file)

        returns: a list

        '''
        hits = []
        log_file.readline()
        for line in log_file:
            #only want to use the term that hit
            if "term: |" in line:
                #get rid of everything but the term
                new_line = line.split('|')
                hits.append(new_line[1].strip('_HUMAN\n'))
        return hits

def Mk_Dict(mito):
    '''
    creates a dictionary where the keys are the names from the symbol and
    synonym columns of the mitocarta.csv file and the values are the categories
    given to each gene from the mitocarta.csv file

    input: mitocarta file (CSV file)

    returns: dictionary

    '''
    mito.readline()
    new_dict = {}
    for line in mito:
        #splits up into the columns
        line = line.split(',')
        #name columns
        symbol = line[0]
        synonym = line[1]
        category = line[3].strip()
        #append each term to the dictionary with its corresponding category
        new_dict[symbol] = category
        each = synonym.split('|')
        for term in each:
            if term != '-':
                new_dict[term] = category
    return new_dict

def Mtch_Dict(hits, new_dict):
    '''
    uses the outputs from the previous two files to create a simplified
    dictionary of only the genes from the gff file that were a match in
    the mitocarta.csv file

    input: list of mitocarta terms that hit and a dictionary of the terms and
    categories

    output: dictionary of all of the matches and their category

    '''
    matches = {}
    for item in hits:
        for key, value in new_dict.items():
            if key in item:
                matches[key] = value
    return matches

def Category_Count(dictionary):
    '''
    counts how many genes in each category

    input: dictonary where the category is the value

    returns: ETC, ribosome, mitochondria, trna counts in this order

    '''
    #initialize counts
    etc_count = 0
    ribosome_count = 0
    mitochondria_count = 0
    trna_count = 0
    for value in dictionary.values():
        #determine which category the term is in and add to the count
        if value == 'ETC':
            etc_count += 1
        elif value == 'Ribosomal':
            ribosome_count += 1
        elif value == 'mitochondria':
            mitochondria_count += 1
        elif value == 'tRNA':
            trna_count += 1
    return etc_count, ribosome_count, mitochondria_count, trna_count

if  __name__ == "__main__":
    with open(sys.argv[1], 'r') as mito:
        with open(sys.argv[2], 'r') as log:
            hit = Hits_List(log)
            new = Mk_Dict(mito)
            matches = Mtch_Dict(hit, new)
            match_count = Category_Count(matches)
            overall_count = Category_Count(new)
            r = np.arange(4)
            width = 0.25
            plt.bar(r, overall_count, color = 'b', width = width, label = 'All NMT Genes')
            plt.bar(r + width, match_count, color = 'g', width = width, label = 'Matched NMT Genes')
            plt.xlabel('Category')
            plt.ylabel('Value')
            plt.title('NMT Gene Categories All vs. Matched')
            plt.xticks(r + width/2, ['ETC', 'Ribosome', 'Mitochondria','tRNA'])
            plt.legend()
            plt.show()
            plt.savefig('my_fig.png')





