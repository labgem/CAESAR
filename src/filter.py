#############
## Imports ##
#############

import argparse
import logging
import re
from pathlib import Path

###############
## Functions ##
###############

def filter_sequence_properties(blast_file, pid, cov, min_len, max_len):
    """Filters sequences based on length, %id and %cov

    Args:
        blast_file (Path): The path of blastp output
        pid (float): %id treshold
        cov (float): %cov treshold
        min_len (int): minmum length
        max_len (int): maximum length

    Returns:
        map_seq (dict): as key the db where come from the sequences and a the
        set of sequence id as value
        
        blastp_lines (str): lines of sequences which pass the filter
    """
    
    blastp_lines = ""
    map_seq = {}
    
    with open(blast_file, "r") as f:
        for line in f:
            ls = line.split()
            p = float(ls[5])
            c = float(ls[6])
            size = int(ls[3])
            
            # filter by %identity
            if p == 100 or p < pid:
                continue
            
            # filter by %coverage
            if c < cov:
                continue
            
            # filter by sequence length
            if size > max_len or size < min_len:
                continue
            
            # sequences from uniprot
            if 'sp|' in ls[2] or 'tr|' in ls[2]:
                seq_id = re.search("\\|(\\w+)\\|", ls[2]).group(1)
                try:
                    map_seq["uniprot"].add(seq_id)
                except KeyError:
                    map_seq["uniprot"] = set()
                    map_seq["uniprot"].add(seq_id)
            
            else:
                seq_id = ls[2]
                source = re.search("matches_(.+?)\.tsv", blast_file.name).group(1)
                
                try:
                    map_seq[source].add(seq_id)
                except KeyError:
                    map_seq[source] = set()
                    map_seq[source].add(seq_id)
            
            blastp_lines += line
    
    return map_seq, blastp_lines

##########
## MAIN ##
##########

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s - %(asctime)s - %(message)s ")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str, default="./", metavar="",
                        help="output directory [default: ./]")
    parser.add_argument("-d", "--data", type=str, metavar="", required=True,
                        help="directory containing blast output")
    parser.add_argument("-q", "--query", type=str, metavar="", required=True,
                        help="set of reference sequences")
    parser.add_argument('--id', type=float, metavar="", default=30.0,
                        help=r"% identity threshold"" [default: 30.0]")
    parser.add_argument("--cov", type=float, metavar="", default=80.0,
                        help=r"% coverage threshold"" [default: 80.0]")
    parser.add_argument("--min", type=int, metavar="", default=200,
                        help="minimum candidate sequence length [default: 200]")
    parser.add_argument("--max", type=int, metavar="", default=1000,
                        help="maximum candidate sequence length [default: 1000]")
    parser.add_argument("--tax", type=str, metavar="", default="ABE",
                        help="selected superkingdoms to filter candidate sequences"
                        ", A=Archae, B=Bacteria, E=Eukaryota [default: ABE]")
    
    args = parser.parse_args()

    blastp_list_file = Path(args.data).glob("matches*.tsv")
    map_seq = {}
    
    
    for blastp_file in blastp_list_file:
        logging.info(blastp_file.name)
        seq, blastp_lines = filter_sequence_properties(blast_file=blastp_file,
                                                           pid=args.id,
                                                           cov=args.cov,
                                                           max_len=args.max,
                                                           min_len=args.min)
        
        map_seq.update(seq)