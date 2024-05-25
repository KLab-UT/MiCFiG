import argparse, os, time
from Bio import Entrez, SeqIO, Seq


class Gene:
    def __init__(self, protein_record):
        """
        Takes a SeqIO protein record and fetches related info 
        """

        self.p_record = protein_record
        self.gene_name = protein_record.description.split()[2]
        self.protein_id = protein_record.id
        
        # Fetch data for cds
        self.cds_record = fetch_protein_cds_record(self.protein_id)
        if self.cds_record is not None:  # Check if cds_record is not None
            self.cds_id = self.cds_record.id
        else:
            self.cds_id = None  # Set cds_id to None if cds_record is None
        # Check data
        self.print_attributes()

    def print_attributes(self):
        """
        Print all attributes of a Gene object
        """
        print("Gene attributes:")
        print(f"Gene Name: {self.gene_name}")
        print(f"Protein ID: {self.protein_id}")
        print(f"CDS ID: {self.cds_id}")
    
def fetch_protein_cds_record(protein_accession, max_retries=3, retry_delay=1):
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
        print(f"Max retries reached for {protein_accession}, appending to problem_seqs.txt")
        with open("problem_seqs.txt", "a") as problem_seqs:
            problem_seqs.write(f"{protein_accession}\n")
        return None

    # Get accession ID encoding sequence associated with protein accession
    for feature in protein_record.features:
        #print(f"protein feature: {feature}\nfeature type: {feature.type}\nfeature qualifiers: {feature.qualifiers}")
        if feature.type == "CDS" and "coded_by" in feature.qualifiers:
            coding_accession = feature.qualifiers['coded_by'][0].split(":")[0]

            # Fetch the cds id and cds sequence
            retry_count = 0
            while retry_count < max_retries:
                try:
                    coding_handle = Entrez.efetch(db="nuccore", id=coding_accession, rettype="fasta", retmode="text")
                    coding_record = SeqIO.read(coding_handle, "fasta")
                    coding_handle.close()
                    return coding_record # return SeqIO record for coding sequence
                    break

                # Retry if error
                except Exception as e:
                    print(f"Error fetching cds record (Attempt {retry_count + 1}): {e}")
                    retry_count += 1
                    time.sleep(retry_delay)

            # Log protein accession to problem seqs if not fetched after max retries
            if retry_count == max_retries:
                print(f"Max retries reached for {protein_accession}, appending to problem_seqs.")
                with open("problem_seqs.txt", "a") as problem_seqs:
                    problem_seqs.write(f"{protein_accession}\n")
                return None

def save_records(records, output_directory):
    """ 
    Save a list of SeqIO records as individual fastas in output_directory named after the record id
    """
    for record in records:
        if record is not None:
            output_file_path = os.path.join(output_directory, f"{record.id}.fasta")
            with open(output_file_path, 'w') as output_handle:
                SeqIO.write(record, output_handle, 'fasta')

def log_gene(gene, log_file):
    line = f"{gene.gene_name},{gene.protein_id},{gene.cds_id}"
    with open(log_file, "a") as file:
        file.write(line + "\n")

def main():
    parser = argparse.ArgumentParser(description="Script used to download CDS fasta sequences as seperate files for every protein sequence in input fasta")
    parser.add_argument("fasta_file", type=str, help="Input fasta containg all protein sequences to download CDS for")
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

    # Iterate over each protein sequence in the input fasta and making a Gene object for each
    genes = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene = Gene(record)
        genes.append(gene)
    
    # Save protein and cds fastas for genes
    protein_records = []
    cds_records = []
    for gene in genes: 
        protein_records.append(gene.p_record)
        cds_records.append(gene.cds_record)
    save_records(protein_records, protein_fastas_dir)
    save_records(cds_records, cds_fastas_dir)

    # Log genes
    log_file = "genes_ids.csv"
    for gene in genes:
        log_gene(gene, log_file)

if __name__ == "__main__":
    main()



