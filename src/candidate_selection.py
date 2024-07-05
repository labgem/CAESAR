#############
## Imports ##
#############

import argparse
import datetime
import logging
import re
import requests
import time
import traceback
import urllib
import yaml
from Bio import Entrez
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from filter import run_seqkit

###########
## Class ##
###########

class Candidate():
    """A class to represent a candidate sequence
    
    Attributes
    ----------
    name (str): sequence id
    protein_fasta (str): protein sequence in fasta format
    nucliec_fasta (str): nucleic sequence in fasta with the header of protein_fasta
    os (str): organism name
    ox (str): tax id
    cds_id (str): EMBL-GenBank-DDBJ_CDS id
    gc (float): %GC
    gc_diff (float): difference between the target %GC and the %GC of the sequence
    source (str): string to indicate the origin database
    query (tuple): information from blast output
    """
    
    def __init__(self, name):
        """Initialize the object

        Args:
            name (str): sequence id
        """
        self.name = name
        self.protein_fasta = ""
        self.nucleic_fasta = ""
        self.os = None
        self.ox = None
        self.cds_id = None
        self.gc = 0.0
        self.gc_diff = 100.0
        self.source = None
        self.query = None
        
    def __str__(self):
        """return a human-readable string to represent the object

        Returns:
            str: formatted string
        """
        return f"{self.name} {self.os} {self.ox} {self.source}"

    def set_organism_name(self, os):
        """Sets the organism name

        Args:
            os (str): organism name
        """
        self.os = os
        
    def set_organism_identifier(self, ox):
        """Sets the tax id

        Args:
            ox (str): tax id
        """
        self.ox = ox
        
    def set_cds_id(self, cds_id):
        """Sets the EMBL-GenBank-DDBJ_CDS id

        Args:
            cds_id (str): EMBL-GenBank-DDBJ_CDS id
        """
        self.cds_id = cds_id
        
    def set_query_info(self, query_info):
        """Sets the information on the query sequence that best matches the
        candidate

        Args:
            query_info (tuple): tuple containing info from the blast output
        """
        self.query = query_info
        
    def set_source(self, source):
        """Sets the source (origin) database

        Args:
            source (str): source (origin) database
        """
        self.source = source
        
    def set_gc(self, gc, gc_target):
        """Sets the %GC of the nucleic sequence

        Args:
            gc (float): %GC of the candidate
            gc_target (float): target %GC
        """
        self.gc = gc
        self.gc_diff = abs(gc_target - gc)
        
    def update_protein_fasta(self, s):
        """Sets or update fasta of protein sequence

        Args:
            s (str): string
        """
        self.protein_fasta += s
        
    def update_nucleic_fasta(self, s):
        """Sets or update fasta of nucleic sequence

        Args:
            s (st): string
        """
        self.nucleic_fasta += s

##############
## Function ##
##############

def read_yaml_config(config_path):
    """Reads the config file in yaml format

    Args:
        config_path (Path): path of the config file

    Returns:
        db_path (dict): path of the databases
        yml (dict): the config file loaded in a dictionnary
    """
    
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

def read_update_file(update_file):
    """Reads the file containing the list of identifiers to exclude,
    together with all their clusters.

    Args:
        update_file (Path): path of the update file

    Returns:
        exclude_cand (set): set containing the ids present in the file
    """
    
    exclude_cand = set()
    with open(update_file, "r") as f:
        for line in f:
            seq_id = ""
            # uniprot identifier
            if "tr|" in line.strip() or "sp|" in line.strip():
                seq_id = re.search("\\|(\\w+)\\|", line).group(1)
            else:
                seq_id = line.strip()
                
            exclude_cand.add(seq_id)
            
    return exclude_cand

def read_clusters(cluster_file, exclude_cand, sources):
    """Reads the cluster tsv file

    Args:
        cluster_file (Path): path of the cluster tsv file
        exclude_cand (set): set containing the ids to exlude

    Returns:
        clusters (dict): the cluster name as key and the list of members as
        value 
    """
    
    clusters = {}
    clusters_sources = {}
    exclude_cluster = set()
    
    with open(cluster_file, "r") as f:
        for line in f:
            split_line = line.split()
            
            # If the cluster name it's in the uniprot identifier format
            if "sp|" in split_line[0] or "tr|" in split_line[0]:
                ref = re.search("\\|(\\w+)\\|", split_line[0]).group(1)
            
            else:
                ref = split_line[0]
            
            # If the member name it's in the uniprot identifier format      
            if "sp|" in split_line[1] or "tr|" in split_line[1]:
                member = re.search("\\|(\\w+)\\|", split_line[1]).group(1)
            else:
                member = split_line[1]
            
            if ref not in clusters:
                clusters[ref] = []
                clusters_sources[ref] = set()
            
            clusters[ref].append(member)
            clusters_sources[ref].add(sources[member])
            
            if ref in exclude_cand:
                exclude_cluster.add(ref)
            elif member in exclude_cand:
                exclude_cluster.add(ref)
    
    # Removes the clusters with their name in exclude_cluster
    for clust in exclude_cluster:
        del clusters[clust]
    
    return clusters, clusters_sources, exclude_cluster

def max_candidate_per_cluster(clusters, value, exclude_cluster):
    """Defines the maximum number of candidate to selection per cluster

    Args:
        clusters (dict): as key the cluster name and as value the set of members
        value (int | float):  the maximum number of candidate or a percentage of
        the cluster size

    Raises:
        TypeError: raised if value isn't an int or a float

    Returns:
        max_cand (dict): as key the cluster name and as value the maximum number
        text (str): info on the clustering
    """
    
    max_cand = {}
    sum_len = 0
    sum_without_singleton = 0
    nb_singleton = 0
    max_len = 0
    max_name = ""
    L = len(clusters)
    
    if isinstance(value, int):
        max_cand = {k:value for k in clusters}
        for clust in clusters:
            max_cand[clust] = value
            l = len(clusters[clust])
            sum_len += l
            if l == 1:
                nb_singleton += 1
            else:
                sum_without_singleton += l
                if l > max_len:
                    max_len = l
                    max_name = clust
    
    elif isinstance(value, float):
        for clust in clusters:
            n = int(round(len(clusters[clust]) / (100/value), 0))
            if n == 0:
                n = 1
            
            max_cand[clust] = n
        
    else:
        raise TypeError(f"'value' should be an int or a float not a: {type(value)}")
    
    mean_len = round(sum_len / L, 1)
    try:
        mean_without_singleton = round((sum_without_singleton / (L - nb_singleton)),1)
    except ZeroDivisionError:
        mean_without_singleton = "NA"
    p_singleton = round((nb_singleton / L) * 100, 1)
    
    text = f"Number of clusters: {L}\n"
    text += f"Number of singleton (cluster with 1 sequence): {nb_singleton}\n"
    text += f"Proportion of singleton: {p_singleton}\n"
    text += f"Largest cluster: {max_name} size: {max_len}\n"
    text += f"Mean cluster size with singleton: {mean_len}\n"
    text += f"Mean cluster size without singleton: {mean_without_singleton}"
    
    logging.info(f"Clusters Statistics:\n{text}")
    if len(exclude_cluster) != 0:
        exclude_str = '\n'.join(exclude_cluster)
        logging.info(f"List of excluded clusters due to the --update option:"
                     f"\n{exclude_str}")
        text += "List of excluded clusters due to the --update option:"
        text += f"\n{exclude_str}"

    text = "## Clustering ##\n" + text + "\n\n"
    
    return max_cand, text

def read_sources_file(sources_file):
    """Reads sources file

    Args:
        sources_file (Path): path of the sources file

    Returns:
        sources (dict): sequence id as key and the source database as value
    """
    
    sources = {}
    
    with open(sources_file, "r") as f:
        for line in f:
            split_line = line.split()
            # uniprot identifier
            if "sp|" in split_line[0] or "tr|" in split_line[0]:
                seq_id = re.search("\\|(\\w+)\\|", split_line[0]).group(1)
                sources[seq_id] = split_line[1]
            else:
                sources[split_line[0]] = split_line[1]
                
    return sources

def read_strain_library(strain_library_path):
    """Reads strain library file

    Args:
        strain_library_path (Path): path of the strain library file

    Returns:
        strain_library (dict): each key have as value a list corresponding to 
        one or two columns of the strain library file
    """
    
    strain_library = {"name":[], "tax_id":[], "resource":[]}
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
                strain_library["resource"].append(f"{split_line[2].strip()},"
                                                   f"{split_line[3].strip()}")
                if mda != False:
                    strain_library["mda"].append(split_line[mda].strip())
                    
    return strain_library

def read_fasta_candidates(fasta_file):
    """Reads a multi fasta file
    
    For each sequence, a Candidate object is created

    Args:
        fasta_file (Path): path of the fasta file

    Returns:
        all_seq (dict): as key the sequence id and as value a Candidate object
    """

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
                
                # Instantiate a Candidate object
                all_seq[seq_id] = Candidate(seq_id)
                # Adds the line to the protein_fasta attribute
                all_seq[seq_id].update_protein_fasta(line)
                
                try:
                    os = re.search("OS=(.+?)(.OX|$)", line).group(1)
                    # Sets the organism name
                    all_seq[seq_id].set_organism_name(os)
                except AttributeError:
                    continue
                
                try:
                    ox = re.search("OX=(\\d+)", line).group(1)
                    # Sets the tax id
                    all_seq[seq_id].set_organism_identifier(ox)
                except AttributeError:
                    continue
                
            else:
                all_seq[seq_id].update_protein_fasta(line)
                
    return all_seq

def get_query_information(all_seq, data_file):
    """Reads blast outputs to get information between candiates and queries
    (reference sequences)

    Args:
        all_seq (dict): as key the sequence id and as value a Candidate object
        data_file (Path): path of the file containing the data

    Returns:
        all_seq (dict): the dict with the updated object
        ref_candidate_count (dict): for each reference, indicate the number of 
        candidate who matched with it
    """
    
    ref_candidate_count = {}

    with open(data_file) as f:
        for line in f:
            split_line = line.split()
            
            seq_id = ""
            query = split_line[0]
            pident = split_line[5]
            qcobhsp = split_line[6]
            positives = round((int(split_line[7]) / int(split_line[4])) *100, 1)
            mismatch = round((int(split_line[8]) / int(split_line[4])) *100, 1)
            gaps = round((int(split_line[9]) / int(split_line[4])) *100, 1)
            e_value = split_line[10]
            query_info = (query, pident, qcobhsp, positives, mismatch,
                          gaps, e_value)
            
            
            if query not in ref_candidate_count:
                ref_candidate_count[query] = {"total":0, "selected":0}
                
            if "sp|" in split_line[2] or "tr|" in split_line[2]:
                seq_id = re.search("\\|(\\w+)\\|", split_line[2]).group(1)
            else:
                seq_id = split_line[2]
            
            if seq_id == "A0A7L4X2M8":
                print(line)
            
            if all_seq[seq_id].query is None:
                
                all_seq[seq_id].set_query_info(query_info)
                ref_candidate_count[query]["total"] += 1
            
            # If a candidate has matched with more than one reference, they are
            # only counted for the reference for whixh they obtained the best scpre
            else:
                queryB = all_seq[seq_id].query[0]
                pidentB = all_seq[seq_id].query[1]
                qcobhspB = all_seq[seq_id].query[2]
                e_valueB = all_seq[seq_id].query[-1]
                
                if float(e_value) < float(e_valueB):
                    all_seq[seq_id].set_query_info(query_info)
                    ref_candidate_count[queryB]["total"] -= 1
                    ref_candidate_count[query]["total"] += 1
                    
                elif float(e_value) == float(e_valueB):
                    if pident > pidentB:
                        all_seq[seq_id].set_query_info(query_info)
                        ref_candidate_count[queryB]["total"] -= 1
                        ref_candidate_count[query]["total"] += 1
                        
                    elif pident == pidentB:
                        if qcobhsp > qcobhspB:
                            all_seq[seq_id].set_query_info(query_info)
                            ref_candidate_count[queryB]["total"] -= 1
                            ref_candidate_count[query]["total"] += 1
                            
    return all_seq, ref_candidate_count

def preselect_candidates(all_seq, clusters, sources, strain_library):
    """Selects candidates according to the content of the strain library

    Args:
        all_seq (dict): as key the sequence id and as value a Candidate object
        clusters (dict): the cluster name as key and the list of members as value 
        sources (dict): sequence id as key and the source database as value
        strain_library (dict): each key have as value a list corresponding to 
        one or two columns of the strain library file

    Returns:
        presel (dict): as key the cluster name and as value a dict with the
        list of selected candidates for each category of selection
        finded (dict): as key the source db and as value the set of preselected 
        candidate ids
        not_finded (set): the candidate with no correspondance in the strain library
    """
    
    presel = {}
    not_finded = set()
    finded = {elem:set() for elem in set(sources.values())}
    
    alt_db = set(elem for elem in finded if elem not in ["uniprot", "nr"])
    
    for clust in clusters:
        # 1st : tax id
        # 2nd : DSM or ATCC id
        # 3rd : organism name
        # 4th : alt db
        
        # the different type of selection
        clust_cand = {"tax_id":[], "home_strain":[], "species":[], "order":[]}
        for db in alt_db:
            clust_cand[db] = []
            
        for cand in clusters[clust]:

            os = all_seq[cand].os
            ox =  all_seq[cand].ox

            atcc_id = None
            dsm_id = None
            
            try:
                atcc_id = "ATCC," + re.search("ATCC(.+?)([^\\s]+)", os).group(2)
            except AttributeError:
                pass
            except TypeError:
                pass
            
            try:
                dsm_id = "DSMZ," + re.search("DSM(.+?)([^\\s]+)", os).group(2)
            except AttributeError:
                pass
            except TypeError:
                pass
            
            # Don't selects organism with at least one of this key word in the name
            p = "unclassified|bacterium|archaeon|fragment|uncultured|_"
            
            find = False
            for key in ["tax_id", "resource", "name"]:
                
                if key == "tax_id":
                    for i in range(len(strain_library[key])):
                        if ox == strain_library[key][i]:
                            clust_cand["tax_id"].append((cand,
                                                         strain_library["name"][i],
                                                         strain_library["tax_id"][i],
                                                         strain_library["resource"][i]))

                            find = True
                            break
                
                elif key == "resource":
                    for i in range(len(strain_library[key])):
                        # DSM
                        if dsm_id == strain_library[key][i]:
                            clust_cand["home_strain"].append((cand,
                                                              strain_library["name"][i],
                                                              strain_library["tax_id"][i],
                                                              strain_library["resource"][i]))
                            find = True
                            break
                        
                        # ATCC
                        elif atcc_id == strain_library[key][i]:
                            clust_cand["home_strain"].append((cand,
                                                              strain_library["name"][i],
                                                              strain_library["tax_id"][i],
                                                              strain_library["resource"][i]))
                            find = True
                            break
            
                elif key == "name":
                    # search for banned keyword
                    if os is not None and re.search(p, os):
                        break
                    
                    elif os is not None:
                        # traditional organism is in the format : genus species
                        os_split = os.split()
                        genus = os_split[0]
                        specie = os_split[1]

                        name = f"{genus}\\s{specie}"
                        for i in range(len(strain_library[key])):

                            # Checks exact match
                            if os == strain_library[key][i]:
                                clust_cand["species"].append((cand,
                                                            strain_library["name"][i],
                                                            strain_library["tax_id"][i],
                                                            strain_library["resource"][i]))
                                find = True
                                break
                            
                            # Checks if the name pattern match
                            elif re.search(name, strain_library[key][i]):
                                clust_cand["species"].append((cand,
                                                                strain_library["name"][i],
                                                                strain_library["tax_id"][i],
                                                                strain_library["resource"][i]))
                                find = True
                                break
                        
                        # No correspondance in the strain library         
                        if find is False:
                            # If a DSM or ATCC id is present, the candidate can
                            # be purchased from the collections
                            if dsm_id is not None:
                                clust_cand["order"].append((cand,os,ox,dsm_id))
                                find = True
                            
                            elif atcc_id is not None:
                                clust_cand["order"].append((cand,os,ox,atcc_id))
                                        
                                find = True
                    
                    # If 'OS=' isn't in the header, this means that the sequence
                    # comes from a db other than uniprot or nr 
                    # (the header format is modified by CAESAR)
                    elif len(alt_db) != 0:
                        clust_cand[sources[cand]].append((cand,None,None,None))
                        find = True
                
                # Adds the candidates in the finded dict  
                if find is True:
                    finded[sources[cand]].add(cand)
                    break
                
            # Adds the candidates in the not_finded set
            if find is False:
                not_finded.add(cand)
            
        presel[clust] = clust_cand
    
    return presel, finded, not_finded

def uniprot_id_mapping(query):
    """Sends a request at the ID Mapping Tool of uniprot

    Map Uniprot id to EMBL-GenBank-DDBJ_CDS

    Args:
        query (list): list of ids

    Raises:
        ValueError: raised if more than 100 000 ids is provided
        RuntimeError: raised if jobStatus is present with a status different from
        'RUNNING' or 'NEW'
        TypeError: raised if the job is completed but no result is returned

    Returns:
        result (requests.Response): the response
    """
    
    if len(query) > 100_000:
        raise ValueError(f"number of ids must be less than 100 000")
    
    POOLING_INTERVAL = 5
    API_URL = "https://rest.uniprot.org/idmapping/"
    
    query_str = ",".join(query)
    
    rpost = requests.post(f"{API_URL}run",
                          data={"from":"UniProtKB_AC-ID",
                                "to":"EMBL-GenBank-DDBJ_CDS",
                                "ids":query_str})
    
    job_id = False
    
    # Get job Id
    while isinstance(job_id, bool):
        try:
            job_id = rpost.json()["jobId"]
        except KeyError:
            time.sleep(POOLING_INTERVAL)
        
    result = False
    
    # Infinite loop
    while True:
        job = requests.get(f"{API_URL}status/{job_id}").json()
        
        # if jobStatus is present, this means that the job has not been completed
        if "jobStatus" in job:
            if job["jobStatus"] == "RUNNING" or job["jobStatus"] == "NEW":
                logging.info(f"Checks job status...")
                time.sleep(POOLING_INTERVAL)
                
            else:
                logging.error(f"ID MAPPING TOOL ERROR:\n\t{job['jobStatus']}")
                raise RuntimeError
            
        else:
            url = f"{API_URL}stream/{job_id}?format=tsv"
            result = requests.get(url)
            break
    
    if isinstance(result, bool):
        raise TypeError(f"Job completed but no results found")
    
    return result
                    
                    
def read_id_mapping_tool_result(result):
    """Reads the response of ID Mapping Tool

    Args:
        result (requests.Response): the response

    Returns:
        cds_map (dict): as key the EMBL-GenBank-DDBJ_CDS id and as value the
        uniport id
    """
    
    cds_map = {}
    
    mapping = result.text.split("\n")[1:-1]
    for line in mapping:
        split_line = line.split()
        cds_map[split_line[1]] = split_line[0]
        
    return cds_map

def efecth_fasta_cds_na(query, mail):
    """Fetchs the nucleic fasta of the cds from the ncbi database

    Args:
        query (list): list of ids
        mail (str): mail

    Raises:
        ValueError: raised if more than 200 ids is provided

    Returns:
        handle (urllib.request.Request): handle to the results
    """
    
    Entrez.email = mail
    
    if len(query) > 200:
        raise ValueError("number of ids must be less or equal than 200")
    
    query_str = ",".join(query)
    
    handle = Entrez.efetch(db="protein", rettype="fasta_cds_na", retmode="text",
                           id=query_str)
    
    return handle

def read_cds_na(handle, all_seq, uniprot_cds_map, target_gc):
    """Reads the handle containing the fasta

    Args:
        handle (urllib.request.Request): _description_
        all_seq (dict): as key the sequence id and as value a Candidate object
        uniprot_cds_map (dict): as key the EMBL-GenBank-DDBJ_CDS id and as value
        the uniport id
        target_gc (float): the target of %GC

    Returns:
        all_seq (dict): the dict with the updated Candidate object
    """
    
    for fasta in handle.read().split(">")[1:]:
        try:
            protein_id = re.search("protein_id=([^\\s]+)\\]", fasta).group(1)
        except:
            continue
        
        # If source db = uniprot
        if protein_id in uniprot_cds_map:
            uniprot_id = uniprot_cds_map[protein_id]
            protein_header = all_seq[uniprot_id].protein_fasta.split("\n")[0]
            # replace the header by the same as in .protein_fasta
            nucleic_fasta = re.sub("^(.+?)\n", protein_header+"\n", fasta)
            all_seq[uniprot_id].update_nucleic_fasta(nucleic_fasta)
            all_seq[uniprot_id].set_cds_id(protein_id)
            # Calculate the %GC
            try:
                gc = compute_gc(nucleic_fasta)
                all_seq[uniprot_id].set_gc(gc, 50.0)
            except ZeroDivisionError:
                continue
        else:
            try:
                protein_header = all_seq[protein_id].protein_fasta.split("\n")[0]
            except Exception as err:
                raise err
            nucleic_fasta = re.sub("^(.+?)\n", protein_header+"\n", fasta)
            all_seq[protein_id].update_nucleic_fasta(nucleic_fasta)
            all_seq[protein_id].set_cds_id(protein_id)
            try:
                gc = compute_gc(nucleic_fasta)
                all_seq[protein_id].set_gc(gc, target_gc)
            except ZeroDivisionError:
                continue
        
    return all_seq

def worker(data):
    """Function run by a worker of ThreadPoolExecutor

    Args:
        data (list): list of all necessary variables

    Returns:
        all_seq (dict): the dict with the updated Candidate object
    """
    query = data[0]
    mail = data[1]
    uniprot_cds_map = data[2]
    all_seq = data[3]
    target_gc = data[4]
    
    handle = efecth_fasta_cds_na(query, mail)
    try:
        all_seq = read_cds_na(handle, all_seq, uniprot_cds_map, target_gc)
    except Exception as err:
        raise err
    
    return all_seq


def compute_gc(fasta):
    """Calculates the %GC

    Args:
        fasta (str): sequence in fasta format

    Returns:
        gc (float): the %GC
    """
    
    seq = "".join(fasta.split("\n")[1:])
    
    count = Counter(seq)
    
    gc = ((count["G"] + count["C"]) /
          (count["A"] + count["T"] + count["G"] + count["C"])) * 100
        
    return round(gc,1)

def select_candidate(presel, all_cand, yml, sources_db, max_cand):
    """Selection the candidates

    Args:
        presel (dict): the preselected candidates for each cluster for each 
        category of selection
        all_cand (dict): as key the sequence id and as value a Candidate object
        yml (dict): the config file loaded in a dictionnary
        sources_db (dict): path of the databases
        max_cand (dict): for each cluster, indicates the maximum number of 
        candidates to selection

    Returns:
        selected_cand_per_clust (dict): as key the cluster name and as value a
        list with a 5-tuple per candidate
    """
    
    selected_cand_per_clust = {}
    
    # Define a priority order for the categories
    category_rank = []
    
    if "candidate_selection" in yml:
        for category in yml["candidate_selection"]:
            
            if category.lower() == "strain_library":
                category_rank.extend(["tax_id", "home_strain", "species"])
            
            elif category.lower() == "order":
                category_rank.append("order")
                
            else:
                if category.lower() in sources_db:
                    category_rank.append(category.lower())
    
    else:
        category_rank.extend(["tax_id", "home_strain", "species", "order"])
        category_rank.extend([cat.lower() for cat in category_rank if cat
                              not in ["uniprot", "nr"]])
    
    # For each cluster selects candidates
    for clust_name, clust_cand in presel.items():
        nb = 0
        selected_cand = []
        
        for category in category_rank:
            if category not in clust_cand:
                continue
            
            # selects the candidates based on the difference between their %GC
            # and a taget %GC
            gc_sorted = [(cand[0], all_cand[cand[0]].gc_diff, cand[1], cand[2],
                          cand[3]) for cand in clust_cand[category]]
            gc_sorted = sorted(gc_sorted, key=lambda x: x[1])
            for cand in gc_sorted:
                selected_cand.append((cand[0], category, cand[2],
                                      cand[3], cand[4]))
                nb += 1
                if nb == max_cand[clust_name]:
                    break
            
            if nb == max_cand[clust_name]:
                break
        
        selected_cand_per_clust[clust_name] = selected_cand
        
    return selected_cand_per_clust

def write_results(selected_cand_per_clust, clusters, clusters_sources, all_cand, outdir):
    """Writes the outputs

    Args:
        selected_cand_per_clust (dict): The selected candiate(s) per cluster
        all_cand (_type_): _description_
        outdir (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    text = "## Candidate Selection ##\nCategory\tNumber of candidates\n"
    cluster_without_cand = 0
    results_categories = {"strain_library":{"table":"", "faa":"", "fna":""},
                          "order":{"table":"", "faa":"", "fna":""}}
    
    nb_cand_per_cat = {c:0 for c in results_categories}
    
    # Loop to complete results_categories
    for cluster in selected_cand_per_clust:
        
        if len(selected_cand_per_clust[cluster]) == 0:
            cluster_without_cand += 1
            continue
        
        sources_in_clust = " ".join(clusters_sources[cluster])
        
        for cand in selected_cand_per_clust[cluster]:
            n_cluster = len(clusters[cluster])
            name = cand[0]  # sequence id
            category = cand[1]  # category/type of selection, e.g : tax_id, home_strain...
            os = all_cand[name].os  # organisme name
            ox = all_cand[name].ox  # tax id
            gc = all_cand[name].gc  # %GC
            cds_id = all_cand[name].cds_id  # EMBL-GenBank-DDBJ_CDS
            
            if all_cand[name].query is not None:
                query_name = all_cand[name].query[0]  # reference sequence name
                pident = all_cand[name].query[1]  # %id
                qcovhsp = all_cand[name].query[2]  # %cov
                positives = all_cand[name].query[3]  # %positives
                mismatch = all_cand[name].query[4]  # %mismatch
                gaps = all_cand[name].query[5]  # %gaps
                e_value = all_cand[name].query[6]  # e-value
            
            # provide a data file with the -d/--data flag is optional
            else:
                query_name, pident, qcovhsp, positives = None, None, None, None
                mismatch, gaps, e_value = None, None, None
            
            sl_org = cand[2]  # organism name of the line that match in strain library
            sl_tax_id = cand[3]  # tax_id of the line that match in strain library
            
            if cand[4] is not None:
                sl_resource = cand[4].split(",")[0]  # ATCC or DSM
                sl_resource_id = cand[4].split(",")[1]  # id in the collection
            else:
                sl_resource = None
                sl_resource_id = None
            
            line = f"{name}\t{cluster}\t{n_cluster}\t{sources_in_clust}\t{os}\t{ox}\t{cds_id}\t{gc}\t"
            line += f"{query_name}\t{pident}\t{qcovhsp}\t{positives}\t"
            line += f"{mismatch}\t{gaps}\t{e_value}\t{category}\t"
            line += f"{sl_org}\t{sl_tax_id}\t{sl_resource}\t{sl_resource_id}\n"
            
            if category in ["tax_id", "home_strain", "species"]:
                nb_cand_per_cat["strain_library"] += 1
                results_categories["strain_library"]["table"] += line
                results_categories["strain_library"]["faa"] += all_cand[name].protein_fasta
                results_categories["strain_library"]["fna"] += all_cand[name].nucleic_fasta
                
            elif category == "order":
                nb_cand_per_cat["order"] += 1
                results_categories["order"]["table"] += line
                results_categories["order"]["faa"] += all_cand[name].protein_fasta
                results_categories["order"]["fna"] += all_cand[name].nucleic_fasta
                
            else:
                if category not in results_categories:
                    results_categories[category] = {"table":"", "faa":"",
                                                    "fna":""}
                    nb_cand_per_cat[category] = 0
                
                nb_cand_per_cat[category] += 1
                results_categories[category]["table"] += line
                results_categories[category]["faa"] += all_cand[name].protein_fasta
                results_categories[category]["fna"] += all_cand[name].nucleic_fasta
    
    # Loop to writes results per category       
    for category in results_categories:
        text += f"{category}\t{nb_cand_per_cat[category]}\n"
        if len(results_categories[category]["table"]) == 0:
            continue
        
        # subdirectories
        category_dir = Path.joinpath(outdir, category)
        if not category_dir.exists():
            category_dir.mkdir()
        
        table_tsv = Path.joinpath(category_dir, "all_candidates.tsv")
        
        header = "Candidate\tCluster\tCluster_size\tSources\tOrganism\tTax_id\tEMBL-GenBank-DDBJ_CDS\t"
        header += "GC\tQuery\tid\tcov\tpositives\tmismatch\tgaps\te-value\t"
        header += "Selection_type\tStrain_library_organism\tStrain_library_tax_id\t"
        header += "Collection\tCollection_id\n"
        
        with open(table_tsv, "w") as f:
            f.write(header)
            f.write(results_categories[category]["table"])
            
        if len(results_categories[category]["faa"]) == 0:
            pass
        else:
            faa_file = Path.joinpath(category_dir, "all_candidates.faa")
            faa_file.write_text(results_categories[category]["faa"])
            
        if len(results_categories[category]["fna"]) == 0:
            pass
        else:
            faa_file = Path.joinpath(category_dir, "all_candidates.fna")
            faa_file.write_text(results_categories[category]["fna"])
            
    text += f"Number of cluster without candidates: {cluster_without_cand}"
    
    summary_file = outdir / "summary.out"
    with summary_file.open("a") as f:
        f.write(text)
    
    return 0

def add_mda(outdir, strain_library):
    
    sl_dir = [d for d in outdir.iterdir() if d.match("strain_library")][0]
    
    sl_tsv = Path.joinpath(sl_dir, "all_candidates.tsv")
    new_line = ""
    with sl_tsv.open("r") as f:
        for i, line in enumerate(f):
            if i != 0:
                taxid = line.split("\t")[4]
                collection = line.split("\t")[-2]
                collection_id = line.split("\t")[-1].strip()

                try:
                    tax_index = strain_library["tax_id"].index(taxid)
                    new_line += line.strip() + "\t" + f"{strain_library['mda'][tax_index]}\n"
                except ValueError:
                    
                    try:
                        collection_str = f"{collection},{collection_id}"
                        ressource_index = strain_library["resource"].index(collection_str)
                        new_line += line.strip() + "\t" + f"{strain_library['mda'][ressource_index]}\n"
                    except ValueError:
                       new_line += line.strip() + "\t" + "None\n"
                       logging.error("did not find this tax_id or collection_id"
                                     " in the strain_library to add the mda"
                                     " column to the output for strain_library"
                                     " candidates")
                        
            else:
                new_line += line.strip() + "\t" + "MDA\n"
                
    sl_tsv.write_text(new_line)
    
    return 0            
        
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
    parser.add_argument("-d", "--data", type=str, metavar="",
                        help="filtered blast output")
    parser.add_argument("-u", "--update", type=str, metavar="",
                        help="file containing a list of identifiers, the groups"
                        " in which they are found will be excluded")
    parser.add_argument("--gc", type=float, metavar="", default=50.0,
                        help="Target GC percentage to decide between candidates")
    ncand_group = parser.add_mutually_exclusive_group()
    ncand_group.add_argument("-n", "--nb-cand", type=int, metavar="", default=1,
                        help="maximum number of candidate per cluster [default: 1]")
    ncand_group.add_argument("--cov-per-cluster", type=float, metavar="",
                             help="Uses a percentage of each cluster as "
                             "maximum number of candidates rather than a given"
                             " number")
    
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    
    start = datetime.datetime.now()
    
    config_path = Path(args.config)
    db_path, yml = read_yaml_config(config_path)
    
    if args.update is not None:
        update_file = Path(args.update)
        exclude_cand = read_update_file(update_file)
    else:
        exclude_cand = set()
        
    sources_file = Path(args.sources)
    sources = read_sources_file(sources_file)
    
    cluster_file = Path(args.clusters)
    clusters, clusters_sources, exclude_cluster = read_clusters(cluster_file,
                                                                exclude_cand,
                                                                sources)
    
    if args.cov_per_cluster is not None:
        max_cand, text = max_candidate_per_cluster(clusters, args.cov_per_cluster,
                                             exclude_cluster)
    else:
        max_cand, text = max_candidate_per_cluster(clusters, args.nb_cand,
                                             exclude_cluster)
        
    summary_file = outdir / "summary.out"
    if summary_file.exists:
        with summary_file.open("a") as f:
            f.write(text)
    else:
        summary_file.write_text(text)
    
    strain_library_file = yml["strain_library"]
    strain_library = read_strain_library(strain_library_file)
    
    fasta_file = Path(args.fasta)
    all_seq = read_fasta_candidates(fasta_file)
    
    if args.data is not None:
        data_file = Path(args.data)
    
        all_seq, ref_candidate_count = get_query_information(all_seq, data_file)

    logging.info("Preselect candidate based on strain library")
    start_presel = datetime.datetime.now()
    presel, finded, not_finded = preselect_candidates(all_seq, clusters,
                                                      sources, strain_library)
    end = datetime.datetime.now()
    logging.info(f"elapsed time: {end - start_presel}")

    for cand in not_finded:
        del all_seq[cand]
    
    del not_finded
    
    if "nr" not in finded:
        finded["nr"] = set()
        
    uniprot_cds_map = {}
    
    # Uniprot ID Mapping Tool
    if "uniprot" in finded:
        logging.info("Uniprot ID Mapping Tool")
        start_mapping = datetime.datetime.now()
        if len(finded["uniprot"]) > 100_000:
            finded["uniprot"] = list(finded["uniprot"])
            for i in range(len(finded["uniprot"]), 100_000):
                query = finded["uniprot"][i:i+100_000]
                result = uniprot_id_mapping(query)
                cds_map = read_id_mapping_tool_result(result)
                uniprot_cds_map.update(cds_map)
                finded["nr"].update(set(cds_map.keys()))
                
        else:
            result = uniprot_id_mapping(finded["uniprot"])
            cds_map = read_id_mapping_tool_result(result)
            uniprot_cds_map.update(cds_map)
            finded["nr"].update(set(cds_map.keys()))
    
    end = datetime.datetime.now()
    logging.info(f"elapsed time: {end - start_mapping}")
    
    # NCBI Entrez to retrieves nucleic sequences
    if "nr" in finded:
        
        finded["nr"] = list(finded["nr"])
        n_seq_id = len(finded["nr"])
        
        data = []
        for i in range(0, n_seq_id, 200):
            query = finded["nr"][i:i+200]
            query_cand = {k:all_seq[k] for k in all_seq if k in query}
            query_map_cds = {}
            
            for k, v in uniprot_cds_map.items():
                if k in query:
                    query_cand[v] = all_seq[v]
                    query_map_cds[k] = v
            data.append((query, yml["mail"], query_map_cds, query_cand, args.gc))

        logging.info(f"NCBI efetch fasta_cds_na")
        start_efecth = datetime.datetime.now()
        with ThreadPoolExecutor(max_workers=3) as executor:
            i = 0
            logging.info(f"{i}/{n_seq_id}")
            
            
            # submit tasks
            # a task = fetch 200 fasta_cds_na and parse them to modify the header
            # and add the fasta in the corresponding candidate object
            futures = [executor.submit(worker, elem) for elem in data]
            
            # iterates over results as soon they are completed
            for future in as_completed(futures):
                try:
                    updated_candidate = future.result()
                    all_seq.update(updated_candidate)
                    i += 200
                    if i <= n_seq_id:
                        logging.info(f"{i}/{n_seq_id}")
                    else:
                        logging.info(f"{n_seq_id}/{n_seq_id}")
                except urllib.error.HTTPError as err:
                    time.sleep(0.5)
                    i += 200
                    if i <= n_seq_id:
                        logging.error(f"{i}/{n_seq_id} {err}")
                    else:
                        logging.error(f"{n_seq_id}/{n_seq_id} {err}")
                except Exception as err:
                    i += 200
                    if i <= n_seq_id:
                        logging.error(f"{i}/{n_seq_id} An error has occured:\n"
                                      f"{traceback.format_exc()}")
                    else:
                        logging.error(f"{n_seq_id}/{n_seq_id} An error has "
                                      f"occured:\n{traceback.format_exc()}")
        
        end = datetime.datetime.now()
        logging.info(f"elapsed time: {end - start_efecth}")
        del data
    
    for key in finded:
        if key in ["uniprot", "nr"]:
            continue
        
        elif "fna" not in db_path[key]:
            continue
        
        key_file_id = outdir / f"{key}_ids.txt"
        if key != "tara":
            key_file_id.write_text("\n".join(finded[key]))
        else:
            key_file_id.write_text("\n".join([s[:-2] for s in finded[key]]))
        
        logging.info(f"run seqkit on {key} to find nucleic sequences")
        start_seqkit = datetime.datetime.now()
        result = run_seqkit(key_file_id, db_path[key]["fna"])
        end = datetime.datetime.now()
        logging.info(f"elapsed time: {end - start_seqkit}")
        
        if result.returncode != 0:
            logging.info("An error has occured during seqkit process to obtain",
                         f" nucleic sequences for {key}:",
                         f" {result.returncode}\n{result.stderr.decode('utf-8')}")
        else:
            for fasta in result.stdout.decode('utf-8').split(">")[1:]:
                seq_id = fasta.split()[0]
                try:
                    all_seq[seq_id].update_nucleic_fasta(fasta)
                    gc = compute_gc(fasta)
                    all_seq[seq_id].set_gc(gc, 50.0)
                except KeyError:
                    seq_id = seq_id + "_1"
                    all_seq[seq_id].update_nucleic_fasta(fasta)
                    gc = compute_gc(fasta)
                    all_seq[seq_id].set_gc(gc, 50.0)
            
            key_file_id.unlink()
    
    logging.info("Selects candidates")
    start_selection = datetime.datetime.now()
    selected_cand_per_clust = select_candidate(presel,
                                               all_seq,
                                               yml,
                                               db_path,
                                               max_cand)
    
    write_results(selected_cand_per_clust, clusters, clusters_sources, all_seq,
                  outdir)
    
    end = datetime.datetime.now()
    logging.info(f"elapsed time: {end - start_selection}")    
    if "mda" in strain_library:
        add_mda(outdir, strain_library)
        
    end = datetime.datetime.now()
    logging.info(f"candidate_selection.py elapsed time: {end - start}")