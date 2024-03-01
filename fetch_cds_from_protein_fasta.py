import argparse
import os 
from Bio import SeqIO, Entrez
import time 
import csv

def fetch_protein_cds(protein_accession, output_directory, max_retries=3, retry_delay=1):
    """
    Fetch the coding sequences for the provided protein accession using NCBI Entrez
    Downloads fasta and returns the accession id to the CDS corresponding to the protein accession
    Sequences that cant be fetched will be logged in problem seqs

    """

    retry_count = 0
    while retry_count < max_retries:
        try:
            handle = Entrez.efetch(db="protein", id=protein_accession, rettype="gb", retmode="text")
            protein_record = SeqIO.read(handle, "genbank")
            handle.close()
            break  # If successful, exit the loop
        except Exception as e:
            print(f"Error fetching protein record (Attempt {retry_count + 1}): {e}")
            retry_count += 1
            time.sleep(retry_delay)

    if retry_count == max_retries:
        print(f"Max retries reached for {protein_accession}, appending to problem_seqs.")
        with open("problem_seqs", "a") as problem_seqs:
            problem_seqs.write(f"{protein_accession}\n")
        return None

    # Get accession ID encoding sequence associated with protein accession
    for feature in protein_record.features:
        #print(f"protein feature: {feature}\nfeature type: {feature.type}\nfeature qualifiers: {feature.qualifiers}")
        if feature.type == "CDS" and "coded_by" in feature.qualifiers:
            coding_accession = feature.qualifiers['coded_by'][0].split(":")[0]


            #fetch the encoding sequence
            retry_count = 0
            while retry_count < max_retries:
                try:
                    coding_handle = Entrez.efetch(db="nuccore", id=coding_accession, rettype="fasta", retmode="text")
                    coding_record = SeqIO.read(coding_handle, "fasta")
                    coding_handle.close()

                    # Write sequence
                    output_file_path = os.path.join(output_directory, f"{coding_record.id}.fasta")
                    print(f"Writing cds: {output_file_path}")
                    SeqIO.write(coding_record, output_file_path, 'fasta')
                    return coding_accession
                    
                    break

                # Retry if error
                except Exception as e:
                    print(f"Error fetching cds record (Attempt {retry_count + 1}): {e}")
                    retry_count += 1
                    time.sleep(retry_delay)

            # Log protein accession to problem seqs if not fetched after max retries
            if retry_count == max_retries:
                print(f"Max retries reached for {protein_accession}, appending to problem_seqs.")
                with open("problem_seqs", "a") as problem_seqs:
                    problem_seqs.write(f"{protein_accession}\n")
                return None

        
def get_accesions(fasta_file, save_directory=None):
    """
    Return a list of accession ids for all sequences in provided fasta_file and saves sequences as seperate files if save_directory not None
    """
    records = SeqIO.parse(fasta_file, "fasta")
    accessions = set()

    for record in records:
        accessions.add(record.id)

        if save_directory:
            output_filepath = os.path.join(save_directory, f"{record.id}.fasta")
            print(f"Saving protein fasta for: {record.id}")
            SeqIO.write(record, output_filepath, 'fasta')
    return accessions


def save_records(records, output_directory):
    """ 
    Save a list of SeqIO records as individual fastas in output_directory named after the record id
    """
    for record in records:
        output_file_path = os.path.join(output_directory, f"{record.id}.fasta")
        with open(output_file_path, 'w') as output_handle:
            SeqIO.write(record, output_handle, 'fasta')

def log_accessions(accessions_dict):
    # Log protein and corresponding cds accesions
    with open("protein_cds_accessions.csv", 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in accessions_dict.items():
            writer.writerow([key, value])
        



def main():
    parser = argparse.ArgumentParser(description="Script used to download CDS fasta sequences as seperate files for every protein sequence in input fasta")
    parser.add_argument("fasta_file", type=str, help="Input fasta containg protein sequences to download CDS for")
    parser.add_argument("fastas_directory", type=str, help="Directory to save fasta files for fetched CDS")
    parser.add_argument("email", type=str, help="Email for NCBI Entrez.")

    args = parser.parse_args()
    Entrez.email = args.email
    fasta_file = args.fasta_file
    fastas_dir = args.fastas_directory

    # Make directories to save seperate fastas
    cds_fastas_dir = os.path.join(fastas_dir, "cds_fastas")
    protein_fastas_dir = os.path.join(fastas_dir, "protein_fastas")
    os.makedirs(cds_fastas_dir, exist_ok=True)
    os.makedirs(protein_fastas_dir, exist_ok=True)


    # Get CDS accessions
    protein_accessions = get_accesions(fasta_file, save_directory=protein_fastas_dir)

    # initilize dictionary to store corresponding accession ids
    protein_cds_accessions = {p_accession: None for p_accession in protein_accessions}



    # Get CDS from protein accessions
    for accession in protein_accessions:
        cds_accession = fetch_protein_cds(accession, cds_fastas_dir)

        # Update dict
        protein_cds_accessions[accession] = cds_accession

    # Log corresponding accessions
    log_accessions(protein_cds_accessions)

    




if __name__ == "__main__":
    main()
