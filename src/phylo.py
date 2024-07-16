#############
## Imports ##
#############

import argparse
import datetime
import logging
import subprocess
from pathlib import Path

##############
## Function ##
##############

def get_representative_sequences(clusters_file):
    """Reads a clusters.tsv file and retrieves the identifiers in the first column. 

    Args:
        clusters_file (Path): the clusters.tsv file

    Returns:
        representative_seqs (set): the set of representative sequences
    """
    
    representative_seqs = set()
    
    with open(clusters_file, "r") as f:
        for line in f:
            representative_seqs.add(line.split()[0])
            
    return representative_seqs

def filter_fasta(fasta_file, representative_seqs):
    """Filters fasta file based on a set of ids

    Args:
        fasta_file (Path): the fasta file
        representative_seqs (set): the set of representative sequences

    Returns:
        filtered_seqs (dict): the filtered sequences
    """
    
    filtered_seqs = {}
    
    with open(fasta_file, "r") as f:
        seq_id = ""
        for line in f:
            if line.startswith(">"):
                seq_id = line.split()[0][1:]
                if seq_id in representative_seqs:
                    filtered_seqs[seq_id] = line
            
            else:
                if seq_id in representative_seqs:
                    filtered_seqs[seq_id] += line
                    
    return filtered_seqs

def write_sequences(filtered_seqs, outdir):
    """Writes sequences

    Args:
        filtered_seqs (dict): the filtered sequences
        outdir (Path): output directory
        
    Returns:
        new_fasta (Path): the ouput fasta file
    """
    
    new_fasta = outdir / "representatives_seqs.fasta"
    
    with new_fasta.open("w") as f:
        for seq in filtered_seqs:
            f.write(filtered_seqs[seq])

    return new_fasta

def run_mafft(fasta_file, outdir):
    """Runs Mafft to generate a multiple sequences alignemnt
    Args:
        fasta_file (Path): the fasta file
        outdir (Path): the output directory

    Returns:
        msa_output (Path): the output msa file
    """
    
    msa_output = outdir / (fasta_file.stem + "_msa.fasta")
    
    MAFFT = "mafft"
    
    command = f"{MAFFT} --auto {fasta_file}"
    
    with msa_output.open("w") as f:
        subprocess.run(command.split(), stdout=f, text=True)
        
    return msa_output

##########
## MAIN ##
##########

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s - %(asctime)s - %(message)s ")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str, default="./", metavar="",
                        help="output directory [default: ./]")
    parser.add_argument("-f", "--fasta", type=str, metavar="", required=True,
                        help="fasta file containing all the candidate sequences")
    parser.add_argument("--clusters", type=str, metavar="",
                        help="clusters tsv file")
    
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    
    start = datetime.datetime.now()
    
    if args.clusters is not None:
        logging.info(f"Reads {args.clusters}")
        representative_seqs = get_representative_sequences(Path(args.clusters))
        logging.info("Retrieves only the representative sequence of each cluster")
        filtered_seqs = filter_fasta(Path(args.fasta), representative_seqs)
        logging.info("Writes the fasta file")
        filtered_fasta = write_sequences(filtered_seqs, outdir)
        
        logging.info("Build MSA with Mafft")
        msa_file = run_mafft(filtered_fasta, outdir)
    else:
        logging.info("Build MSA with Mafft")
        msa_file = run_mafft(Path(args.fasta), outdir)
        
    
    
    end = datetime.datetime.now()