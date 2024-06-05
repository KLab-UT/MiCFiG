#!/usr/bin/env python3

from Bio import SeqIO
import argparse

def extract_gene_name(in_fasta, out_txt):
    in_fasta = open(in_fasta, "r")
    out_txt = open(out_txt, "w")
    for record in SeqIO.parse(in_fasta, "fasta"):
        header = record.description
        head_split = header.split("\t")
        if len(head_split) == 3:
            extract = head_split[0] + "\t" + head_split[2]
            out_txt.write(extract + "\n")
            print("Extract complete!")
    in_fasta.close()
    out_txt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, help = "Input FASTA file")
    parser.add_argument("-o", "--output", type = str, help = "Output text file for extracted gene names")
    args = parser.parse_args()
    extract_gene_name(args.input, args.output)

if __name__=="__main__":
    main()
    print("Done!")
