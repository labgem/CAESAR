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
        pass
    
    end = datetime.datetime.now()