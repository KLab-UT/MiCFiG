#
#
# process_blast.py 
# script to process xml blast results and make a bed file for the top hits, also creates a csv for query_gene, pident
# 
# 
import argparse, os
import xml.etree.ElementTree as ET
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
import re
import csv

class BlastHit:
    """
    Manage blast hits from ET element
    """
    def __init__(self, hit_element, query_gene, query_id, query_length):
        self.query_gene = query_gene
        self.query_id = query_id
        self.query_len = int(query_length)
        self.accession = hit_element.find("Hit_accession").text
        self.definition = hit_element.find("Hit_def").text
        self.length = int(hit_element.find("Hit_len").text)

        print(f"Hit number: {(hit_element.find('Hit_num').text)}")
        self.hit_num = int((hit_element.find('Hit_num').text))
        

        # Iterate over Hsp elements and create Hsp objects if #1 hit
        self.top_hsp = None
        if self.hit_num == 1: 
            for hsp_element in hit_element.findall("Hit_hsps/Hsp"):
                hsp = self.create_hsp_object(hsp_element)
                if hsp.hsp_num == 1:
                    self.top_hsp=hsp
                    break
        self.overall_identity = self.calculate_percent_identity() # TOP HSP identity
        
    
    def create_hsp_object(self, hsp_element):
        """
        Create Hsp object for blast hits from Hsp XML element
        """
        return Hsp(
            hsp_num=int(hsp_element.find("Hsp_num").text),
            bit_score=float(hsp_element.find("Hsp_bit-score").text),
            score=int(hsp_element.find("Hsp_score").text),
            evalue=float(hsp_element.find("Hsp_evalue").text),
            query_start=int(hsp_element.find("Hsp_query-from").text),
            query_end=int(hsp_element.find("Hsp_query-to").text),
            hit_start=int(hsp_element.find("Hsp_hit-from").text),
            hit_end=int(hsp_element.find("Hsp_hit-to").text),
            query_frame=int(hsp_element.find("Hsp_query-frame").text),
            hit_frame=int(hsp_element.find("Hsp_hit-frame").text),
            identity=int(hsp_element.find("Hsp_identity").text),
            positive=int(hsp_element.find("Hsp_positive").text),
            gaps=int(hsp_element.find("Hsp_gaps").text),
            align_len=int(hsp_element.find("Hsp_align-len").text),
            qseq=hsp_element.find("Hsp_qseq").text,
            hseq=hsp_element.find("Hsp_hseq").text,
            midline=hsp_element.find("Hsp_midline").text
        )
    def calculate_percent_identity(self):
        """
        Calculate the percentage identity between query and hit sequences from top HSP
        """
        qseq = self.top_hsp.qseq
        hseq = self.top_hsp.hseq
    
        if len(qseq) != len(hseq):
            raise ValueError("Query and hit sequences must be of the same length")
        
        identity_count = sum(1 for q, h in zip(qseq, hseq) if q == h)
        percent_identity = (identity_count / len(qseq)) * 100
        return percent_identity

    def make_bed(self, output_path):
        """
        Make bed file for blast hit, with an entry for each hsp
        """
        with open(output_path, 'w') as bed_out:
            bed_line = "\t".join([self.definition, str(self.top_hsp.hit_start), str(self.top_hsp.hit_end), self.query_id])
            bed_out.write(bed_line)
    
    def get_query_gene(self, gene_ids_csv):
        """
        Uses the query id to find gene name in csv with genes and ids. Gene needs to be first field in csv w/ ids
        """
        with open(gene_ids_csv, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if self.query_id in row:
                    return row[0]  # Return the first field of the row where the ID is found
        return None  # Return None if the query ID is not found in the CSV


class Hsp:
    """
    Manage HSP information
    An HSP is a high scoring pair representing a specific aligngment within a sequence hit by blast
    """
    def __init__(self, hsp_num, bit_score, score, evalue, query_start, query_end,
                 hit_start, hit_end, query_frame, hit_frame, identity,
                 positive, gaps, align_len, qseq, hseq, midline):
        self.hsp_num = hsp_num
        self.bit_score = bit_score
        self.score = score
        self.evalue = evalue
        self.query_start = query_start
        self.query_end = query_end
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.query_frame = query_frame
        self.hit_frame = hit_frame
        self.identity = identity
        self.positive = positive
        self.gaps = gaps
        self.align_len = align_len
        self.qseq = qseq
        self.hseq = hseq
        self.midline = midline
    

def get_query_gene(query_id, gene_ids_csv):
    """
    Uses the query id to find gene name in csv with genes and ids. Gene needs to be first field in csv w/ ids
    """
    with open(gene_ids_csv, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if query_id in row:
                return row[0]  # Return the first field of the row where the ID is found
    return None  # Return None if the query ID is not found in the CSV

def log_identity(blast_hit, identities_csv):
    """
    Make a csv containing the overall identity for the top Blast hit of a gene
    """
    entry = (blast_hit.query_gene, blast_hit.overall_identity)
    # Append the entry to the CSV file
    with open(identities_csv, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(entry)


def parse_blast(tree_root, gene_ids_csv, bed_output_path, identities_output_csv):
    """
    """
    query_id = tree_root.find(".//Iteration_query-def").text
    query_gene_name = get_query_gene(query_id, gene_ids_csv)
    query_length = tree_root.find(".//Iteration_query-len").text

    #Process blast hits
    top_hit = None
    for hit_element in tree_root.findall(".//Hit"):
        hit = BlastHit(hit_element, query_gene_name, query_id, query_length)
        if hit.hit_num == 1:
            top_hit = hit
            top_hit.make_bed(bed_output_path)
            log_identity(top_hit, identities_output_csv)
            break
    # If top_hit is None, save an empty file at bed_output_path
    if top_hit is None:
        with open(bed_output_path, 'w') as empty_file:
            pass

    

def main(xml_file, gene_ids_csv, bed_output_path, identities_output_csv):
    # Read XML blast results 
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Process blast hits
    parse_blast(root, gene_ids_csv, bed_output_path, identities_output_csv)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process XML blast results and create a BED file for the top hits, also creates a CSV for query_gene, pident")
    parser.add_argument("xml_file", help="Path to the XML file containing blast results")
    parser.add_argument("gene_ids_csv", help="Path to the CSV file containing gene IDs")
    parser.add_argument("--bed_output_path", default=".", help="path to save the BED files (default: current directory)")
    parser.add_argument("--identities_output_csv", default="identities.csv", help="Path to the output CSV file for gene identities (default: identities.csv)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    main(args.xml_file, args.gene_ids_csv, args.bed_output_path, args.identities_output_csv)