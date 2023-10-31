#!/usr/bin/env python

import sys

'''
This script reads through the sorted nmt genes and removes any blank lines and
any duplicate lines
'''

def find_dup(nmt_genes, new_nmtgenes):
    nmt_genes_nodups = set()
    with open(nmt_genes, 'r') as nmt_in:
        for line in nmt_in:
            nmt_genes_nodups.add(line)
    print(nmt_genes_nodups)
    with open(new_nmtgenes, 'w') as nmt_out:
        nmt_out.write(''.join(nmt_genes_nodups))


if __name__ == '__main__':
    nmt_genes = sys.argv[1]
    new_nmtgenes = sys.argv[2]
    find_dup(nmt_genes, new_nmtgenes)




