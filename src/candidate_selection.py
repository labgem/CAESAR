#############
## Imports ##
#############

import argparse
import logging
import re
import yaml
from pathlib import Path

##############
## Function ##
##############

def read_yaml_config(config_path):
    
    with open(config_path, "r") as f:
        yml = yaml.safe_load(f)
        
    db_path = {}
    
    for key in yml:
        if "db" in key:
            db_name = key[:-3].lower()
            db_path[db_name] = {}
            
            for db_format in yml[key]:
                db_path[db_name].update(db_format)
    
    return db_path, yml

def read_clusters(cluster_file):
    
    clusters = {}
    
    with open(cluster_file, "r") as f:
        for line in f:
            split_line = line.split()
            
            # If the cluster name it's in the uniprot header format
            if "sp|" in split_line[0] or "tr|" in split_line[0]:
                ref = re.search("\\|(\\w+)\\|", split_line[0]).group(1)
                if clusters.get(ref) is None:
                    clusters[ref] = []
            
            else:
                ref = split_line[0]
                if clusters.get(ref) is None:
                    clusters[ref] = []
            
            # If the member name it's in the uniprot header format      
            if "sp|" in split_line[1] or "tr|" in split_line[1]:
                member = re.search("\\|(\\w+)\\|", split_line[1]).group(1)
                clusters[ref].append(member)
            else:
                clusters[ref].append(split_line[1])
    
    return clusters

##########
## MAIN ##
##########

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s - %(asctime)s - %(message)s ")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str, default="./", metavar="",
                        help="output directory [default: ./]")
    parser.add_argument("-c", "--config", type=str, metavar="", required=True,
                        help="the configuration file in yaml format")
    parser.add_argument("-f", "--fasta", type=str, metavar="", required=True,
                        help="fasta file containing all the candidate sequences")
    parser.add_argument("--clusters", type=str, metavar="", required=True,
                        help="clusters tsv file")
    parser.add_argument("--sources", type=str, metavar="", required=True,
                        help="sources file indicating the source database of each"
                        "sequences")
    
    args = parser.parse_args()
    
    config_path = Path(args.config)
    db_path, yml = read_yaml_config(config_path)
    
    cluster_file = Path(args.clusters)
    clusters = read_clusters(cluster_file)