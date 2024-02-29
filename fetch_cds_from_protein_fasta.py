import argparse
import os 
from Bio import SeqIO, Entrez
import time 

def fetch_protein_cds(protein_accession, output_directory, max_retries=3, retry_delay=1):
    """
    Fetch the coding sequences for the provided protein accession using NCBI Entrez
    """

    retry_count = 0
    while retry_count < max_retries:
        try:
            handle = Entrez.efetch(db="protein", id=protein_accession, rettype="gb", retmode="text")
            protein_record = SeqIO.read(handle, "genbank")
            handle.close()
            print(f"Successfully fetched protein record for {protein_accession}")
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
    print(f"Features in protein record: {len(protein_record.features)}")
    for feature in protein_record.features:
        #print(f"protein feature: {feature}\nfeature type: {feature.type}\nfeature qualifiers: {feature.qualifiers}")
        if feature.type == "CDS" and "coded_by" in feature.qualifiers:
            coding_accession = feature.qualifiers['coded_by'][0].split(":")[0]
            #print(f"Coding accesion: {coding_accession}")

            #fetch the encoding sequence
            retry_count = 0
            while retry_count < max_retries:
                try:
                    coding_handle = Entrez.efetch(db="nuccore", id=coding_accession, rettype="fasta", retmode="text")
                    coding_record = SeqIO.read(coding_handle, "fasta")
                    coding_handle.close()

                    # nucleotides = {"A", "T", "G", "C"}
                    #sequence = str(coding_record.seq)

                    # # Replace characters not in nucleotides or '-' with '-'
                    # sanitized_sequence = ''.join(c if c in nucleotides or c == '-' else '-' for c in sequence)

                    # Write sequence
                    output_file_path = os.path.join(output_directory, f"{coding_record.id}.fasta")
                    SeqIO.write(coding_record, output_file_path, 'fasta')

                except Exception as e:
                    print(f"Error fetching protein record (Attempt {retry_count + 1}): {e}")
                    retry_count += 1
                    time.sleep(retry_delay)
            if retry_count == max_retries:
                print(f"Max retries reached for {protein_accession}, appending to problem_seqs.")
                with open("problem_seqs", "a") as problem_seqs:
                    problem_seqs.write(f"{protein_accession}\n")
                return None

        
def get_accesions(fasta_file):
    """
    Return a list of accession ids for all sequences in provided fasta_file
    """
    records = SeqIO.parse(fasta_file, "fasta")
    accessions = set()

    for record in records:
        accessions.add(record.id)
    return accessions


def save_records(records, output_directory):
    """ 
    Save a list of SeqIO records as individual fastas in output_directory named after the record id
    """
    for record in records:
        output_file_path = os.path.join(output_directory, f"{record.id}.fasta")
        with open(output_file_path, 'w') as output_handle:
            SeqIO.write(record, output_handle, 'fasta')
        


def main():
    parser = argparse.ArgumentParser(description="Script used to download CDS fasta sequences as seperate files for every protein sequence in input fasta")
    parser.add_argument("fasta_file", type=str, help="Input fasta containg protein sequences to download CDS for")
    parser.add_argument("fastas_directory", type=str, help="Directory to save fasta files for fetched CDS")
    parser.add_argument("email", type=str, help="Email for NCBI Entrez.")

    args = parser.parse_args()
    Entrez.email = args.email
    fasta_file = args.fasta_file
    fastas_dir = args.fastas_directory


    # Get CDS accessions
    protein_accessions = get_accesions(fasta_file)
    cds_records = []


    # Get CDS from CDS accessions
    for accession in protein_accessions:
        fetch_protein_cds(accession, fastas_dir)




if __name__ == "__main__":
    main()
