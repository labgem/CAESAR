#############
## Imports ##
#############

import argparse
import logging
import re
import requests
import time
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
    
    def __init__(self, name):
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
        return f"{self.name} {self.os} {self.ox} {self.source}"

    def set_organism_name(self, os):
        self.os = os
        
    def set_organism_identifier(self, ox):
        self.ox = ox
        
    def set_cds_id(self, cds_id):
        self.cds_id = cds_id
        
    def set_query_info(self, query_info):
        self.query = query_info
        
    def set_source(self, source):
        self.source = source
        
    def set_gc(self, gc, gc_goal):
        self.gc = gc
        self.gc_diff = abs(gc_goal - gc)
        
    def update_protein_fasta(self, s):
        self.protein_fasta += s
        
    def update_nucleic_fasta(self, s):
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

def get_query_information(all_seq, data_file):
    
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
            
            if all_seq[seq_id].query is None:
                
                all_seq[seq_id].set_query_info(query_info)
                ref_candidate_count[query]["total"] += 1
                
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
    
    presel = {}
    not_finded = set()
    finded = {elem:set() for elem in set(sources.values())}
    
    alt_db = set(elem for elem in finded if elem not in ["uniprot", "nr"])
    
    for clust in clusters:
        # 1st : tax id
        # 2nd : organism name (excact match, include strain etc)
        # 3rd : DSM or ATCC id in strain_library
        # 4th : organism name (same specie)
        # 5th : alt db
        
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
                        if dsm_id == strain_library[key][i]:
                            clust_cand["home_strain"].append((cand,
                                                              strain_library["name"][i],
                                                              strain_library["tax_id"][i],
                                                              strain_library["resource"][i]))
                            find = True
                            break
                        
                        elif atcc_id == strain_library[key][i]:
                            clust_cand["home_strain"].append((cand,
                                                              strain_library["name"][i],
                                                              strain_library["tax_id"][i],
                                                              strain_library["resource"][i]))
                            find = True
                            break
            
                elif key == "name":
                    if os is not None and re.search(p, os):
                        break
                    
                    elif os is not None:
                        os_split = os.split()
                        genus = os_split[0]
                        specie = os_split[1]

                        name = f"{genus}\\s{specie}"
                        for i in range(len(strain_library[key])):

                            if os == strain_library[key][i]:
                                clust_cand["species"].append((cand,
                                                            strain_library["name"][i],
                                                            strain_library["tax_id"][i],
                                                            strain_library["resource"][i]))
                                find = True
                                break
                            
                            elif re.search(name, strain_library[key][i]):
                                clust_cand["species"].append((cand,
                                                                strain_library["name"][i],
                                                                strain_library["tax_id"][i],
                                                                strain_library["resource"][i]))
                                find = True
                                break
                                    
                        if find is False:
                            if dsm_id is not None:
                                clust_cand["order"].append((cand,os,ox,dsm_id))
                                find = True
                            
                            elif atcc_id is not None:
                                clust_cand["order"].append((cand,os,ox,atcc_id))
                                        
                                find = True
                    
                    elif len(alt_db) != 0:
                        clust_cand[sources[cand]].append((cand,None,None,None))
                        find = True
                    
                if find is True:
                    finded[sources[cand]].add(cand)
                    break
            
            if find is False:
                not_finded.add(cand)
            
        presel[clust] = clust_cand
    
    return presel, finded, not_finded

def uniprot_id_mapping(query):
    
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
    
    cds_map = {}
    
    mapping = result.text.split("\n")[1:-1]
    for line in mapping:
        split_line = line.split()
        cds_map[split_line[1]] = split_line[0]
        
    return cds_map

def efecth_fasta_cds_na(query, mail):
    
    Entrez.email = mail
    
    if len(query) > 200:
        raise ValueError("number of ids must be less or equal than 200")
    
    query_str = ",".join(query)
    
    handle = Entrez.efetch(db="protein", rettype="fasta_cds_na", retmode="text",
                           id=query_str)
    
    return handle

def read_cds_na(handle, all_seq, uniprot_cds_map):
    
    
    for fasta in handle.read().split(">")[1:]:
        try:
            protein_id = re.search("cds_([^\\s_]+)", fasta).group(1)
        except:
            continue
        
        if protein_id in uniprot_cds_map:
            uniprot_id = uniprot_cds_map[protein_id]
            protein_header = all_seq[uniprot_id].protein_fasta.split("\n")[0]
            nucleic_fasta = re.sub("^(.+?)\n", protein_header+"\n", fasta)
            all_seq[uniprot_id].update_nucleic_fasta(nucleic_fasta)
            all_seq[uniprot_id].set_cds_id(protein_id)
            try:
                gc = compute_gc(nucleic_fasta)
                all_seq[uniprot_id].set_gc(gc, 50.0)
            except ZeroDivisionError:
                continue
        else:
            protein_header = all_seq[protein_id].protein_fasta.split("\n")[0]
            nucleic_fasta = re.sub("^(.+?)\n", protein_header+"\n", fasta)
            all_seq[protein_id].update_nucleic_fasta(nucleic_fasta)
            all_seq[protein_id].set_cds_id(protein_id)
            try:
                gc = compute_gc(nucleic_fasta)
                all_seq[protein_id].set_gc(gc, 50.0)
            except ZeroDivisionError:
                continue
        
    return all_seq

def worker(data):
    query = data[0]
    mail = data[1]
    uniprot_cds_map = data[2]
    all_seq = data[3]
    
    handle = efecth_fasta_cds_na(query, mail)
    all_seq = read_cds_na(handle, all_seq, uniprot_cds_map)
    
    return all_seq


def compute_gc(fasta):
    
    seq = "".join(fasta.split("\n")[1:])
    
    count = Counter(seq)
    
    
    gc = ((count["G"] + count["C"]) /
          (count["A"] + count["T"] + count["G"] + count["C"])) * 100
        
    return round(gc,1)

def select_candidate(presel, all_cand, yml, sources_db, n=1):
    
    selected_cand_per_clust = {}
    
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
    
    for clust_name, clust_cand in presel.items():
        nb = 0
        selected_cand = []
        
        for category in category_rank:
            if category not in clust_cand:
                continue
            
            gc_sorted = [(cand[0], all_cand[cand[0]].gc_diff, cand[1], cand[2],
                          cand[3]) for cand in clust_cand[category]]
            gc_sorted = sorted(gc_sorted, key=lambda x: x[1])
            for cand in gc_sorted:
                selected_cand.append((cand[0], category, cand[2],
                                      cand[3], cand[4]))
                nb += 1
                if nb == n:
                    break
            
            if nb == n:
                break
        
        selected_cand_per_clust[clust_name] = selected_cand
        
    return selected_cand_per_clust

def write_results(selected_cand_per_clust, all_cand, outdir):
    
    results_categories = {"strain_library":{"table":"", "faa":"", "fna":""},
                          "order":{"table":"", "faa":"", "fna":""}}
    
    for cluster in selected_cand_per_clust:
        
        if len(selected_cand_per_clust[cluster]) == 0:
            continue
        
        for cand in selected_cand_per_clust[cluster]:
            name = cand[0]
            category = cand[1]
            os = all_cand[name].os
            ox = all_cand[name].ox
            gc = all_cand[name].gc
            genbank_id = all_cand[name].cds_id
            
            if all_cand[name].query is not None:
                query_name = all_cand[name].query[0]
                pident = all_cand[name].query[1]
                qcovhsp = all_cand[name].query[2]
                positives = all_cand[name].query[3]
                mismatch = all_cand[name].query[4]
                gaps = all_cand[name].query[5]
                e_value = all_cand[name].query[6]
            else:
                query_name, pident, qcovhsp, positives = None, None, None, None
                mismatch, gaps, e_value = None, None, None
            
            sl_org = cand[2]
            sl_tax_id = cand[3]
            if cand[4] is not None:
                sl_resource = cand[4].split(",")[0]
                sl_resource_id = cand[4].split(",")[1]
            else:
                sl_resource = None
                sl_resource_id = None
            
            line = f"{cluster}\t{name}\t{os}\t{ox}\t{genbank_id}\t{gc}\t"
            line += f"{query_name}\t{pident}\t{qcovhsp}\t{positives}\t"
            line += f"{mismatch}\t{gaps}\t{e_value}\t{category}\t"
            line += f"{sl_org}\t{sl_tax_id}\t{sl_resource}\t{sl_resource_id}\n"
            
            if category in ["tax_id", "home_strain", "species"]:
                results_categories["strain_library"]["table"] += line
                results_categories["strain_library"]["faa"] += all_cand[name].protein_fasta
                results_categories["strain_library"]["fna"] += all_cand[name].nucleic_fasta
                
            elif category == "order":
                results_categories["order"]["table"] += line
                results_categories["order"]["faa"] += all_cand[name].protein_fasta
                results_categories["order"]["fna"] += all_cand[name].nucleic_fasta
                
            else:
                if category not in results_categories:
                    results_categories[category] = {"table":"", "faa":"",
                                                    "fna":""}
                    
                results_categories[category]["table"] += line
                results_categories[category]["faa"] += all_cand[name].protein_fasta
                results_categories[category]["fna"] += all_cand[name].nucleic_fasta
                
    for category in results_categories:
        if len(results_categories[category]["table"]) == 0:
            continue
        
        category_dir = Path.joinpath(outdir, category)
        if not category_dir.exists():
            category_dir.mkdir()
            
        table_tsv = Path.joinpath(category_dir, "all_candidates.tsv")
        
        header = "Cluster\tCandidate\tOrganism\tTax_id\tEMBL-GenBank-DDBJ_CDS\t"
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
    
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    
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
    
    if args.data is not None:
        data_file = Path(args.data)
    
        all_seq, ref_candidate_count = get_query_information(all_seq, data_file)

    presel, finded, not_finded = preselect_candidates(all_seq, clusters,
                                                      sources, strain_library)

    for cand in not_finded:
        del all_seq[cand]
    
    del not_finded
    
    if "nr" not in finded:
        finded["nr"] = set()
        
    uniprot_cds_map = {}
    
    # Uniprot ID Mapping Tool
    if "uniprot" in finded:
        logging.info("Uniprot ID Mapping Tool")
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
    
    logging.info("done")
    
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

            data.append((query, yml["mail"], query_map_cds, query_cand))
        
        logging.info(f"NCBI efetch fasta_cds_na")
        with ThreadPoolExecutor(max_workers=3) as executor:
            i = 0
            logging.info(f"{i}/{n_seq_id}")
            
            
            # submit tasks
            # a task = fetch 200 fasta_cds_na and parse them to modify the header
            # and add the fasta in the corresponding candidate object
            futures = [executor.submit(worker, elem) for elem in data]
            
            # iterates over results as soon they are completed
            for future in as_completed(futures):
                updated_candidate = future.result()
                all_seq.update(updated_candidate)
                i += 200
                if i <= n_seq_id:
                    logging.info(f"{i}/{n_seq_id}")
                else:
                    logging.info(f"{n_seq_id}/{n_seq_id}")
        
        logging.info("done")
        del data
    
    for key in finded:
        if key in ["uniprot", "nr"]:
            continue
        
        elif "fna" not in db_path[key]:
            continue
        
        key_file_id = outdir / f"{key}_ids.txt"
        key_file_id.write_text("\n".join(finded[key]))
        
        result = run_seqkit(key_file_id, db_path[key]["fna"])
        
        if result.returncode != 0:
            logging.info("An error has occured during seqkit process to obtain",
                         f" nucleic sequences for {key}:",
                         f" {result.returncode}\n{result.stderr.decode('utf-8')}")
        else:
            for fasta in result.stdout.decode('utf-8').split(">")[1:]:
                seq_id = fasta.split()[0]
                all_seq[seq_id].update_nucleic_fasta(fasta)
                gc = compute_gc(fasta)
                all_seq[seq_id].set_gc(gc, 50.0)
                key_file_id.unlink()
    
    selected_cand_per_clust = select_candidate(presel,
                                               all_seq,
                                               yml,
                                               db_path)
    
    write_results(selected_cand_per_clust, all_seq, outdir)