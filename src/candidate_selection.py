#############
## Imports ##
#############

import argparse
import logging
import re
import yaml
from pathlib import Path

###########
## Class ##
###########

class Candidate():
    
    def __init__(self, name):
        self.name = name
        self.protein_fasta = ""
        self.nucleic_fasta = ""
        self.os = None
        self.ox = None
        self.embl_id = None
        self.gc = 0.0
        self.gc_diff = 100.0
        self.source = None
        self.query = None
        
    def __str__(self):
        return f"{self.name} {self.os} {self.ox} {self.source}"

    def set_organism_name(self, os):
        self.os = os
        
    def set_organism_identifier(self, ox):
        self.ox = ox
        
    def set_embl_id(self, embl_id):
        self.embl_id = embl_id
        
    def set_query_info(self, query_info):
        self.query = query_info
        
    def set_source(self, source):
        self.source = source
        
    def set_gc(self, gc, gc_goal):
        self.gc = gc
        self.gc_diff = abs(gc_goal - gc)
        
    def update_protein_fasta(self, s):
        self.protein_fasta += s
        
    def update_nucleic_sequence(self, s):
        self.nucleic_fasta += s

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

def read_sources_file(sources_file):
    
    sources = {}
    
    with open(sources_file, "r") as f:
        for line in f:
            split_line = line.split()
            if "sp|" in split_line[0] or "tr|" in split_line[0]:
                seq_id = re.search("\\|(\\w+)\\|", split_line[0]).group(1)
                sources[seq_id] = split_line[1]
            else:
                sources[split_line[0]] = split_line[1]
                
    return sources

def read_strain_library(strain_library_path):
    
    strain_library = {"name":[], "tax_id":[], "ressource":[]}
    mda = False
    with open(strain_library_path, "r") as fin:
        for i, line in enumerate(fin):
            if i == 0:
                if "MDA" in line.strip().split("\t"):
                    mda = line.strip().split("\t").index("MDA")
                    strain_library["mda"] = []
                continue
            else:
                split_line = line.split("\t")
                strain_library["name"].append(split_line[0].strip())
                strain_library["tax_id"].append(split_line[1].strip())
                strain_library["ressource"].append(f"{split_line[2].strip()},"
                                                   f"{split_line[3].strip()}")
                if mda != False:
                    strain_library["mda"].append(split_line[mda].strip())
                    
    return strain_library

def read_fasta_candidiates(fasta_file):

    all_seq = {}
    
    with open(fasta_file, "r") as f:
        seq_id = ""
        for line in f:
            if line.startswith(">"):
                if "sp|" in line or "tr|" in line:
                    seq_id = re.search("\\|(\\w+)\\|", line).group(1)
                # If other header
                else:
                    seq_id = line[1:].split()[0]
                    
                all_seq[seq_id] = Candidate(seq_id)
                all_seq[seq_id].update_protein_fasta(line)
                
                try:
                    os = re.search("OS=(.+?)(.OX|$)", line).group(1)
                    all_seq[seq_id].set_organism_name(os)
                except AttributeError:
                    continue
                
                try:
                    ox = re.search("OX=(\\d+)", line).group(1)
                    all_seq[seq_id].set_organism_identifier(ox)
                except AttributeError:
                    continue
                
            else:
                all_seq[seq_id].update_protein_fasta(line)
                
    return all_seq
        
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
                        help="sources file indicating the sources database of each"
                        "sequences")
    
    args = parser.parse_args()
    
    config_path = Path(args.config)
    db_path, yml = read_yaml_config(config_path)
    
    cluster_file = Path(args.clusters)
    clusters = read_clusters(cluster_file)
    
    sources_file = Path(args.sources)
    sources = read_sources_file(sources_file)
    
    strain_library_file = yml["strain_library"]
    strain_library = read_strain_library(strain_library_file)
    
    fasta_file = Path(args.fasta)
    all_seq = read_fasta_candidiates(fasta_file)