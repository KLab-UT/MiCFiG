from Bio import Entrez
import os
import subprocess
import argparse

def download_assemblies_for_order(tax_id, email, save_path):
    """
    Download assemblies and GFF files for a specified taxonomic group using Entrez from NCBI.

    Args:
        tax_id (int): Taxonomic identifier for the group.
        email (str): Your email address for NCBI API usage.
        save_path (str): Path to save the order directory.
    """
    Entrez.email = email

    # Set the batch size and initial start index
    retmax = 500
    retstart = 0

    # Create a directory for the group
    order_directory = os.path.join(save_path, f"Assemblies_{tax_id}")
    if not os.path.exists(order_directory):
        os.makedirs(order_directory)

    while True:
        try:
            # Search for the assemblies associated with the specified tax_id
            id_list = fetch_assembly_ids(tax_id, retmax, retstart)
            if not id_list:
                break

            # Fetch the FTP paths for the assemblies and download all available files
            download_assemblies(id_list, order_directory)
            retstart += retmax
        except Exception as e:
            print(f"Error: {e}")

def fetch_assembly_ids(tax_id, retmax, retstart):
    """
    Fetch assembly IDs for a specified taxonomic group.

    Args:
        tax_id (int): Taxonomic identifier for the group.
        retmax (int): Maximum number of records to retrieve.
        retstart (int): Index of the first record to retrieve.

    Returns:
        list: List of assembly IDs.
    """
    search_term = f"txid{tax_id}[Organism] AND latest[filter] AND (latest_refseq[filter] OR latest_genbank[filter])"
    handle = Entrez.esearch(db="assembly", term=search_term, retmax=retmax, retstart=retstart)
    record = Entrez.read(handle)
    return record["IdList"]

def download_assemblies(id_list, order_directory):
    """
    Download assemblies and GFF files for a list of assembly IDs.

    Args:
        id_list (list): List of assembly IDs.
        order_directory (str): Path to the order directory.
    """
    for assembly_id in id_list:
        try:
            # Fetch information about the assembly
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            record = Entrez.read(handle)
            ftp_path_genbank = record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]
            ftp_path_refseq = record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
            assembly_name = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyName"]

            # Create a directory for the assembly
            directory_name = os.path.join(order_directory, assembly_name.replace(" ", "_"))
            if not os.path.exists(directory_name):
                os.makedirs(directory_name)

            # Download GenBank files using wget
            os.chdir(directory_name)
            download_command_genbank = f"wget -r -nd -A '*' {ftp_path_genbank}"
            subprocess.call(download_command_genbank, shell=True)

            # Download GFF files using wget
            download_command_gff = f"wget -r -nd -A '*.gff*' {ftp_path_refseq}"
            subprocess.call(download_command_gff, shell=True)

            os.chdir('..')
        except Exception as e:
            print(f"Error downloading assembly {assembly_id}: {e}")

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description='Download assemblies and GFF files for a specified taxonomic group using Entrez from NCBI.')
    parser.add_argument('--tax_id', type=int, required=True, help='Taxonomic identifier for the group')
    parser.add_argument('--email', type=str, required=True, help='Your email address')
    parser.add_argument('--save_path', type=str, required=True, help='Path to save the order directory')
    args = parser.parse_args()

    # Main execution
    download_assemblies_for_order(args.tax_id, args.email, args.save_path)

