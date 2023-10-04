#!/usr/bin/env python

import sys

'''
This script reads through the sorted nmt genes and removes any blank lines and
any duplicate lines
'''

def find_dup(nmt_genes, new_nmtgenes):
    with open(nmt_genes, 'r') as nmt_in:
        lines = nmt_genes.readline()
    clean_lines = set(line.strip() for line in lines if line.strip())
    with open(new_nmtgenes, 'w') as nmt_out:
        file.write('\n'.join(clean_lines))


if __name__ == '__main__':
    nmt_genes = home/ebunch/whiptail_nmt_variation/test.txt
    new_nmtgenes = home/ebunch/whiptail_nmt_variation/nmtgenes.txt
    find_dup(nmt_genes, new_nmtgenes)




