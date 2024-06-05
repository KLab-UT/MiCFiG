# MiCFiG

MiCFiG (MitoCarta Filter GFF) is a bioinformatics tool that uses the MitoCarta genomic database to to filter GFF files of nonmodel organisms.

# Contents

- [Context](#context)
- [Methods](#methods)
  - [Scripts](#scirpts)
  - [Pipeline](#pipeline)
- [Walkthrough](#walkthrough)

## Context

A common practice in bioinformatics is to annotate genomic files. Annotation files (e.g. BED & GFF) contain the locations of known genes in genomic files, acting as a road map for genes of interest. By their definiton, the genomes model organisms are incredibly well-annotated - lots of people have taken to the call and had their work double- & triple-checked. However, the pool of non-model organisms without fully annotated genomes is vast. 

MiCFiG aims to make annotation of species with incomprehensive genetic labeling a much easier process. It does so by using BLAST to identify & filter genetic features in non-model organisms based on homologies with model species. The MitoCarta database, a catalog of mitochondrial & nuclear protein-encoding genes in *Homo sapiens* & *Mus musculus*, provides a diving board for such comparisons.

## Methods

The MiCFiG pipeline is intended for non-model organisms with genomes that are not thoroughly annotated with informative gene names. By using the MitoCarta database, the seqeunces of known human genes are used as queries for non-model organisms to identify mitochondrial-targeting genes. These are then extracted in a filtered annotation file & named, providing annotated genes for downstream analysis.

### Scripts

* split_fasta.py
  * Takes a comprehensive FASTA file (in this case, taken from the MitoCarta database) and splits it into individual FASTA files - one for each sequence

* blast_script.sh
  * Loops through FASTA files and performs a BLAST search

* process_blast_results.sh
  * For each BLAST hit, grabs the minimum subject start and maximum subject end from the top rated HSP score
  * process_blast.py
    * Processes BLAST results and makes a BED file for the top hits based on HSP score

* filter_gff.py
  * Identifies overlaps between BED & GFF files
  * Determines which gene from model organism belongs with the genes in non-model organism

* add_gene_name.sh ( I think this goes here in the pipeline )
  * Adds MitoCarta gene name to GFF file
    * Accounts for potential differences in genetic nomenclature between species
  * extract_gene_name.py
    * Extracts gene name & sequence ID from MitoCarta-produced FASTA file & puts it in a lovely text file

* process_gff.py
  * Creates complete & annotated GFF file for non-model organism
  * Keeps track of genes not found to have a sister homology

### Pipeline

![alt text](Pipeline.png "Generalized MiCFiG Pipeline")

## Walkthrough

Work in progress - TGK
