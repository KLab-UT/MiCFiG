#! /usr/bin/env python3

import sys
import argparse

def parse_filtered_blast(results):
    # Assuming the filtered BLAST results is in the format "Query / Chromosome / Start / Stop
    parsed_results = {}
    with open(results, "r") as infile:
        for line in infile:
            items = line.strip().split("\t")
            key = items[1]
            parsed_results[key] = items
    return parsed_results

def parse_gff(gff):
    parsed_gff = {}i
    line_num = 1
    with open(gff, "r") as infile:
        for line in infile:
            items = line.strip().split("\t")
            key = line_num
            parsed_gff[key] = items
            line_num += 1
    reutrn parsed_gff

def within_range(a, b):
    # a and b are tuples, the first item being a minimum value and the second
    # item being a maximum value
    if a[0] >= b[0] or a[1] <= b[1]:
        return True
    else:
        return False

def find_matches(blast_results, gff_results):
    matches = []


def make_new_gff():
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast",
            type = str,
            help = "Input file of filtered BLAST results")
    parser.add_argument("-g", "--gff",
            type = str,
            help = "Input GFF file")
    parser.add_argument("-o", "--output",
            type = str,
            help = "Output file name")
    args = parser.parse_args()
    results = parse_filtered_blast(args.blast)


if __name__=="__main__":
    main()
