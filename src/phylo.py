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