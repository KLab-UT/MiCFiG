import argparse
import os 
from Bio import SeqIO, Entrez

def fetch_protein_cds(protein_accession):
    """
    Fetch the coding sequences for the provided protein accession using NCBI Entrez
    """

    # Fetch the protein record using the protein accession ID
    try:
        handle = Entrez.efetch(db="protein", id=protein_accession, rettype="gb", retmode="text")
        protein_record = SeqIO.read(handle, "genbank")
        #print(f"Protein Record: {protein_record}")
    except Exception as e:
        print(f"Error fetching protein record: {e}")
        return None

    # Get accession ID encoding sequence associated with protein accession
    print(f"Features in protein record: {len(protein_record.features)}")
    for feature in protein_record.features:
        #print(f"protein feature: {feature}\nfeature type: {feature.type}\nfeature qualifiers: {feature.qualifiers}")
        if feature.type == "CDS" and "coded_by" in feature.qualifiers:
            coding_accession = feature.qualifiers['coded_by'][0].split(":")[0]
            #print(f"Coding accesion: {coding_accession}")

            #fetch the encoding sequence
            try:
                coding_handle = Entrez.efetch(db="nuccore", id=coding_accession, rettype="fasta", retmode="text")
                coding_record = SeqIO.read(coding_handle, "fasta")
                coding_handle.close()

                # nucleotides = {"A", "T", "G", "C"}
                #sequence = str(coding_record.seq)

                # # Replace characters not in nucleotides or '-' with '-'
                # sanitized_sequence = ''.join(c if c in nucleotides or c == '-' else '-' for c in sequence)

                return coding_record

            except Exception as e:
                print(f"An error occurred: {e}")
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
        cds_record = fetch_protein_cds(accession)
        #print(f"Record ID: {cds_record.id}\nRecord Description: {cds_record.description}\nRecord seq: {cds_record.seq}")
        cds_records.append(cds_record)

    # Save CDS as individual fastas
    save_records(cds_records, fastas_dir)

    




if __name__ == "__main__":
    main()
