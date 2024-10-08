#############
## Imports ##
#############

import argparse
import datetime
import logging
import re
import requests
import subprocess
import sys
import textwrap
import time
import yaml
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

###############
## Functions ##
###############

def read_yaml_config(config_path):
    """Reads the config file in yaml format

    Args:
        config_path (Path): path of the config file

    Returns:
        db_path (dict): path of the databases
        mail (str): mail used for the ncbi entrez
    """
    
    with open(config_path, "r") as f:
        yml = yaml.safe_load(f)
        
    db_path = {}
    mail = ""
    
    for key in yml:
        if "db" in key:
            db_name = key[:-3].lower()
            db_path[db_name] = {}
            
            for db_format in yml[key]:
                db_path[db_name].update(db_format)
        
        if key == "mail":
            mail = yml[key]
    
    return db_path, mail

def filter_sequence_properties(data_file, pid, cov, min_len, max_len):
    """Filters sequences based on length, %id and %cov

    Args:
        data_file (Path): The path of blastp output
        pid (float): %id treshold
        cov (float): %cov treshold
        min_len (int): minmum length
        max_len (int): maximum length

    Returns:
        map_seq (dict): as key the db where come from the sequences and a the
        set of sequence id as value
        
        map_blastp (dict): as key the sequence ID and the corresponding lines
        in blastp output as value
    """
    
    map_blastp = {}
    map_seq = {}
    
    with open(data_file, "r") as f:
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
            
            # for other databases, we use the name of the blastp output file to
            # gets the source database
            else:
                seq_id = ls[2]
                source = re.search("matches_(.+?)\\.tsv", data_file.name).group(1).lower()
                
                try:
                    map_seq[source].add(seq_id)
                except KeyError:
                    map_seq[source] = set()
                    map_seq[source].add(seq_id)
    
            try:
                map_blastp[seq_id] += line
            except KeyError:
                map_blastp[seq_id] = line
    
    return map_seq, map_blastp

def filter_hmm_properties(data_file, cov, min_len, max_len, treshold, db_path):
    """Filters hits based on sequence lenght and hmm coverage

    Args:
        data_file (Path): the path of hits.domtbl
        cov (float): %cov treshold
        min_len (int): minmum length
        max_len (int): maximum length
        treshold (float): score treshold
        db_path (dict): path of the databases

    Returns:
        map_seq (dict): as key the db where come from the sequences and a the
        set of sequence id as value
        
        map_hmm (dict): as key the sequence ID and the corresponding lines
        in blastp output as value
    """
    
    map_seq = {}
    map_hmm = {}
        
    with open(data_file, 'r') as f:
        
        dom_length = []  # length of each hit domain for a target
        ali_length = 0  #  sum of aligned domain length
        
        for line in f:
            if line.startswith("#"):
                continue
            
            split_line = line.split()
            name = split_line[0]
            tlen = int(split_line[2])
            query = split_line[3]
            qlen = int(split_line[5])
            e_value = float(split_line[6])
            score = float(split_line[7])
            hit_id = int(split_line[9])
            total_hit = int(split_line[10])
            hmm_from = int(split_line[15])
            hmm_to = int(split_line[16])
            
            if tlen > max_len or tlen < min_len:
                continue
            
            dom_length.append((hmm_to - hmm_from) +1)
            
            if hit_id == total_hit:
                ali_length = sum(dom_length)
                dom_length = []  # reset the list for the next target
                
                if ali_length > qlen:
                    continue
                
                hmm_cov = round(ali_length / qlen, 2) * 100
                
                if (hmm_cov >= cov) and (score >= treshold):
                    if 'sp|' in name or "tr|" in name:
                        seq_id = re.search("\\|(\\w+)\\|", name).group(1)
                        try:
                            map_seq["uniprot"].add(seq_id)
                        except KeyError:
                            map_seq["uniprot"] = set()
                            map_seq["uniprot"].add(seq_id)
                            
                        data = f"{seq_id}\t{tlen}\t{query}\t{score:.1f}\t{e_value:.2e}"
                        data += f"\t{hmm_cov}\n"
                        try:
                            map_hmm[seq_id] += data
                        except KeyError:
                            map_hmm[seq_id] = data
                        
                    else:
                        source = ""
                        for db in db_path:
                            if "faa" in db_path[db]:
                                prefix = Path(db_path[db]['faa']).stem
                                if prefix in data_file.stem:
                                    source = db
                                    break

                        try:
                            map_seq[source].add(name)
                        except KeyError:
                            map_seq[source] = set()
                            map_seq[source].add(name)
                            
                        data = f"{name}\t{tlen}\t{query}\t{score:.1f}\t{e_value:.2e}"
                        data += f"\t{hmm_cov}\n"
                        try:
                            map_hmm[name] += data
                        except KeyError:
                            map_hmm[name] = data              
                                    
    return map_seq, map_hmm

def uniprotkb_accessions(query, data_type):
    """Requests uniprotkb rest api

    Args:
        query (list): list of sequence IDs
        data_type (str): fasta to retrieve fasta or lineage to retrieve lineage

    Returns:
        results (requests.Response): the server response
    """
    
    query_str = ",".join(query)
    
    if data_type == "fasta":
        url = f"https://rest.uniprot.org/uniprotkb/accessions?accessions={query_str}&format=fasta&size=500"
    
    elif data_type == "lineage":
        url = f"https://rest.uniprot.org/uniprotkb/accessions?accessions={query_str}&fields=accession,lineage_ids&format=tsv&size=500"
        
    results = requests.get(url)
    
    time.sleep(0.3)
    
    return results

def superkindgoms_filter(lineage, taxon_pattern):
    """Filter candidates by suoerkingdom
    
    Archaea = 2157
    Bacteria = 2
    Eukaryota = 2759

    Args:
        lineage (str): lineage tsv from uniprot
        taxon_pattern (str): taxon id to search separate by | if more than one

    Returns:
        list_id: list of uniqud ids corresponding to the selected candidate
    """
    
    set_id = set()
    
    for line in lineage.split("\n"):
        superkingdom = re.search(f"({taxon_pattern})\\s.superkingdom.", line)
        if superkingdom is not None:
            set_id.add(line.split("\t")[0])
        else:
            continue
        
    list_id = list(set_id)
    
    return list_id

def run_seqkit(list_ids_file, db_path):
    """Runs seqkit grep

    Args:
        list_ids_file (Path): file containing the ids to retrieve
        db_path (Path): database path

    Returns:
        result (CompletedProcess): the completed process
    """
    
    command = f"seqkit grep -f {list_ids_file} {db_path}"

    result = subprocess.run(command.split(), check=True, capture_output=True)
    
    return result

def cloaca_filter(sequences, pattern):
    """Filters cloaca sequences based on superkingdoms and removes non-complete
    genes

    Args:
        sequences (str): all sequences
        pattern (str): selection pattern

    Returns:
        ids_checked (list): ids selected
    """
    
    ids_checked = []
    
    for line in sequences.split("\n"):
        if line.startswith(">") and pattern in line:
            ids_checked.append(line[1:].split()[0])
      
    return ids_checked

def tara_filter(sequences, pattern):
    """Filters TARA sequences, removes non-complete gene

    Args:
        sequences (str): all sequences
        pattern (): selection pattern

    Returns:
        ids_checked (list): ids selected
        fasta (str): sequences selected
    """
    
    seq_id = ""
    ids_checked = []
    fasta = ""
    
    
    for line in sequences.split("\n"):
        if line.startswith(">"):
            if pattern in line:
                seq_id = line[1:].split()[0]
                ids_checked.append(seq_id)
                fasta += line + "\n"
            else:
                seq_id = ""
        
        elif seq_id != "":
            fasta += line + "\n"
            
    return ids_checked, fasta

def efecth_ncbi_genpept(query, mail):
    """Requests ncbi Entrez to retrieve genpept

    Args:
        query (list): list of sequence IDs
        mail (str): email for Entrez

    Returns:
        handle (HTTTPResponse): the server response
    """
    
    Entrez.email = mail
    
    query_str = ",".join(query)
    
    handle = Entrez.efetch(db="protein", id=query_str, rettype="gp",
                           retmode="xml")
    
    return handle

def read_genpept(handle, taxon_pattern):
    """Filters and get sequences

    Args:
        handle (HTTTPResponse): the response
        taxon_pattern (str): taxon name separated by | if more than one

    Returns:
        ids_checked (list): ids selected
        fasta (str): sequences selected
    """
    
    ids_checked = []
    fasta = ""
    
    record = Entrez.read(handle)
    
    for elem in record:
        seq_id = ""
        single_fasta = ""
        tax_id = ""
        cds = False
        
        if re.match(taxon_pattern, elem['GBSeq_taxonomy']):
            for feature in elem['GBSeq_feature-table']:
                if feature['GBFeature_key'] == "CDS":
                    cds = True
                if feature['GBFeature_key'] == "source":
                    for quals in feature['GBFeature_quals']:
                        if quals['GBQualifier_name'] == "db_xref":
                            tax_id = quals['GBQualifier_value']
                            try:
                                tax_id = re.search("taxon:(\\d+)", tax_id).group(1)
                            except AttributeError:
                                continue
        
        if cds is True:
            seq_id = elem['GBSeq_accession-version']
            single_fasta = f">{seq_id} {elem['GBSeq_definition']}"
            single_fasta = single_fasta.replace("[", "OS=")
            single_fasta = single_fasta.replace("]", "")
            single_fasta += f" OX={tax_id}\n"
            single_fasta += textwrap.fill(elem['GBSeq_sequence'].upper(), 60)
            
            ids_checked.append(seq_id)
            fasta += single_fasta + "\n"
        
    return ids_checked, fasta

def worker(data):
    
    list_ids = data[0]
    mail = data[1]
    taxon_pattern = data[2]
    handle = efecth_ncbi_genpept(list_ids, mail)
    ids_selected, seq_selected = read_genpept(handle, taxon_pattern)
    
    return ids_selected, seq_selected
    
def write_data(fasta_file, fasta, out_file, out_lines, sources_file, sources_text, i):
    """Writes data

    Args:
        fasta_file (Path): fasta file
        fasta (str): selected sequences
        out_file (Path): file containing the output of blast for the selected 
        sequences
        out_lines (str): blast output lines of selected sequences
        sources_file (Path): sources.txt file
        sources_text (Path): text to writes in sources_file
        i (int): int
    """
    
    if i == 0:
        fasta_file.write_text(fasta)    
        out_file.write_text(out_lines)
        sources_file.write_text(sources_text)
            
    else:
        with open(fasta_file, "a") as f:
            f.write(fasta)
                
        with open(out_file, "a") as f:
            f.write(out_lines)
            
        with open(sources_file, "a") as f:
            f.write(sources_text)

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
    parser.add_argument("-d", "--data", type=str, metavar="", required=True,
                        help="directory containing blast output or a blast"
                        " output")
    parser.add_argument("-q", "--query", type=str, metavar="", required=True,
                        help="set of reference sequences")
    parser.add_argument('--id', type=float, metavar="", default=30.0,
                        help="percentage identity threshold"" [default: 30.0]")
    parser.add_argument("--cov", type=float, metavar="", default=80.0,
                        help="percentage coverage threshold"" [default: 80.0]")
    parser.add_argument("--min", type=int, metavar="", default=200,
                        help="minimum candidate sequence length [default: 200]")
    parser.add_argument("--max", type=int, metavar="", default=1000,
                        help="maximum candidate sequence length [default: 1000]")
    parser.add_argument("--tax", type=str, metavar="", default="ABE",
                        help="selected superkingdoms to filter candidate sequences"
                        ", A=Archae, B=Bacteria, E=Eukaryota [default: ABE]")
    parser.add_argument("--hmm-score", type=float, default=0.0)
    
    args = parser.parse_args()

    config_path = Path(args.config).absolute()
    db_path, mail = read_yaml_config(config_path)
    
    outdir = Path(args.outdir).absolute()
    if not outdir.exists():
        outdir.mkdir()

    logging.info(f"selected superkingdoms: {args.tax}")

    if Path(args.data).is_dir():
        data_list_file = [f for ext in ['tsv', "domtbl"]
                          for f in Path(args.data).glob("*."+ext)]
    elif Path(args.data).is_file():
        data_list_file = [Path(args.data)]
    else:
        logging.error(f"'{args.data}' was not found")
        sys.exit(1)
        
    map_seq = {}
    data_map = {}
    text = "## Blastp ##\n"
    
    start = datetime.datetime.now()
    
    # Reads blastp output and filters on %id %cov and size
    for data_file in data_list_file:
        logging.info(data_file.name)
        if data_file.suffix == ".tsv":
            seq, map_data_line = filter_sequence_properties(data_file=data_file,
                                                            pid=args.id,
                                                            cov=args.cov,
                                                            max_len=args.max,
                                                            min_len=args.min)
        
        elif data_file.suffix == ".domtbl":
            seq, map_data_line = filter_hmm_properties(data_file,
                                                      cov=args.cov,
                                                      max_len=args.max,
                                                      min_len=args.min,
                                                      treshold=args.hmm_score,
                                                      db_path=db_path)
        # Key correspond to the sources database
        if len(seq) == 0:
            logging.error("no sequence retained after the filter")
            sys.exit(1)
            
        key = list(seq.keys())[0]
        if key in map_seq:
            map_seq[key].update(seq[key])
        else:
            map_seq.update(seq)
        data_map.update(map_data_line)
    
    for i, key in enumerate(map_seq):
        fasta_file = outdir / "filtered_sequences.fasta"
        fasta = ""
        
        data_filtered_file = outdir / "filtered_data.tsv"
        data_filtered_lines = ""
        
        sources_file = outdir / "sources.txt"
        sources_text = ""
        
        if key == "uniprot":
            map_seq[key] = list(map_seq[key])
            n_seq_id = len(map_seq[key])
            
            taxon_superkingdom = {"A":"2157", "B":"2", "E":"2759", "AB":"2|2157",
                                  "BE":"2|2759", "AE":"2157|2759",
                                  "ABE":"2|2157|2759"}
            taxon_pattern = taxon_superkingdom[''.join(sorted(args.tax))]
            
            # Filters based on superkingdoms
            logging.info("taxonomy filtering for uniprot sequences")
            ids_checked = []

            # Get lineage and filter it
            for j in range(0, n_seq_id, 500):
                logging.info(f"{j}/{n_seq_id}")
                lineage = uniprotkb_accessions(map_seq[key][j:j+500], "lineage")
                    
                ids_checked.extend(superkindgoms_filter(lineage.text,
                                                        taxon_pattern))
            logging.info(f"{n_seq_id}/{n_seq_id}")
            
            text += f"{len(ids_checked)} unique ids from {key} after filtering\n"
            
            # Get fasta for sequences that pass the filter
            logging.info("retrieves fasta for uniprot sequences")
            for j in range(0, len(ids_checked), 500):
                logging.info(f"{j}/{len(ids_checked)}")
                sequences = uniprotkb_accessions(ids_checked[j:j+500], "fasta")
                fasta += sequences.text
            logging.info(f"{len(ids_checked)}/{len(ids_checked)}")
                
            # Get the blastp lines for the checked sequences and update the
            # lines to write in sources.txt                
            for ic in ids_checked:
                data_filtered_lines += data_map[ic]
                sources_text += f"{ic} {key}\n"

                # Directly get the fasta for each sequence id
                '''
                logging.info(f"retrieves fasta for uniprot sequences")
                for j in range(0, n_seq_id, 500):
                    logging.info(f"{j}/{n_seq_id}")
                    sequences = uniprotkb_accessions(map_seq[key][j:j+500], "fasta")
                    fasta += sequences.text
                logging.info(f"{n_seq_id}/{n_seq_id}")
                
                for ic in map_seq[key]:
                    data_filtered_lines += data_map[ic]
                    sources_text += f"{ic} {key}\n"
                '''
        elif re.match("nr", key):
            map_seq[key] = list(map_seq[key])
            n_seq_id = len(map_seq[key])
            
            taxon_superkingdom = {"A":"Archaea", "B":"Bacteria", "E":"Eukaryota",
                                  "AB":"Archaea|Bacteria",
                                  "BE":"Bacteria|Eukaryota",
                                  "AE":"Archaea|Eukaryota",
                                  "ABE":"Archaea|Bacteria|Eukaryota"}
            taxon_pattern = taxon_superkingdom[''.join(sorted(args.tax))]
            
            ids_checked = []
            
            logging.info("taxonomy filtering and retrieve fasta for nr sequences")
            
            # Build data structure
            data = [[map_seq[key][j:j+200], mail, taxon_pattern]
                    for j in range(0, len(map_seq[key]), 200)]

            # Use ThreadPoolExecturor with context manager
            with ThreadPoolExecutor(max_workers=3) as executor:
                j = 0
                logging.info(f"{j}/{n_seq_id}")
                
                # Submit tasks
                # a task = fetch 200 (max) genpept at xml format from ncbi Entrez
                # and parse the xml to filter by the taxonomy, the cds 
                # availability and get the sequences
                futures = [executor.submit(worker, elem) for elem in data]
                
                # Iterates over results as soon they are completed
                for future in as_completed(futures):
                    ids_selected, seq_selected = future.result()
                    ids_checked.extend(ids_selected)
                    fasta += seq_selected
                    j += 200
                    if j <= n_seq_id:
                        logging.info(f"{j}/{n_seq_id}")
                    else:
                        logging.info(f"{n_seq_id}/{n_seq_id}")
            
            text += f"{len(ids_checked)} unique ids from {key} after filtering\n"
                  
            for ic in ids_checked:
                data_filtered_lines += data_map[ic]
                sources_text += f"{ic} {key}\n"
            
        elif re.match("cloaca", key):
            # We keep only complete gene
            # adjusts the pattern to the selected superkingdom
            if ''.join(sorted(args.tax)) in ["AB", "ABE"]:
                taxon_pattern = "complete"
            elif ''.join(sorted(args.tax)) in ["A", "AE"]:
                taxon_pattern = "archaea complete"
            elif ''.join(sorted(args.tax)) in ["B", "BE"]:
                taxon_pattern = "bacteria complete"
            
            # We need the nucleic sequences (fna) to filter cloaca
            # sequence ids based on superkingdoms
            # First, we need to create a file containing all the ids
            cloaca_prefilter_file = outdir / "cloaca_prefilter.txt"
            cloaca_prefilter_file.write_text("\n".join(map_seq[key]))
            
            # Runs seqkit to obtain nucleic sequences
            logging.info("search cloaca nucleic sequence with seqkit")
            result = run_seqkit(cloaca_prefilter_file, db_path[key]["fna"])
            
            cloaca_prefilter_file.unlink()
            
            ids_checked = []
            
            if result.returncode != 0:
                logging.error("An error has occured during seqkit process to "
                              "obtain nucleic squences for Cloaca:"
                              f" {result.returncode}\n{result.stderr.decode('utf-8')}")
            else:
                # Specific filter for cloaca sequences
                logging.info("filters cloaca sequences based on superkingdoms "
                             "and removes non-complete gene")
                ids_checked.extend(cloaca_filter(result.stdout.decode("utf-8"),
                                                 taxon_pattern))
            
            # New file for the filtered ids
            cloaca_filtered_file = outdir / "cloaca_filtered.txt"
            cloaca_filtered_file.write_text("\n".join(ids_checked))
        
            logging.info("retrieves protein sequences")
            # Runs seqkit to obtain amino acid sequences
            result = run_seqkit(cloaca_filtered_file, db_path[key]["faa"])
            if result.returncode != 0:
                logging.error("An error has occured during seqkit process to "
                              "obtain nucleic squences for Cloaca:"
                              f" {result.returncode}\n{result.stderr.decode('utf-8')}")
            else:
                fasta = result.stdout.decode("utf-8")
                logging.info("done")
            
            text += f"{len(ids_checked)} unique ids from {key} after filtering\n"
            
            for ic in ids_checked:
                data_filtered_lines += data_map[ic]
                sources_text += f"{ic} {key}\n"
                
            cloaca_filtered_file.unlink()
            
        elif re.match("tara", key):
            if args.tax != "E":
                
                # We create a file containing all the ids to retrieve
                tara_prefilter_file = outdir / "tara_prefilter.txt"
                tara_prefilter_file.write_text("\n".join(map_seq[key]))
                
                # Runs seqkit to obtain the protein sequences
                logging.info("retrieves TARA protein sequences")
                result = run_seqkit(tara_prefilter_file, db_path[key]["faa"])
                
                tara_prefilter_file.unlink()
                
                seq_id = ""
                ids_checked = []
                
                if result.returncode != 0:
                    logging.error("An error has occured during seqkit process to "
                                  "obtain protein squences for TARA:"
                                  f" {result.returncode}\n{result.stderr.decode('utf-8')}")
                else:
                    # Specific filter for TARA sequences
                    logging.info("filters TARA sequences, removes non-complete gene")
                    ids_selected, fasta = tara_filter(result.stdout.decode("utf-8"),
                                                      "gene_type:complete")
                    
                    ids_checked.extend(ids_selected)
                
                text += f"{len(ids_checked)} unique ids from {key} after filtering\n"
                
                for ic in ids_checked:
                    data_filtered_lines += data_map[ic]
                    sources_text += f"{ic} {key}\n"
        
        else:
            # No specific filter for other databases
            key_prefilter_file = outdir / f"{key}_prefilter.txt"
            key_prefilter_file.write_text("\n".join(map_seq[key]))
            
            logging.info(f"search protein sequences in {key}")
            result = run_seqkit(key_prefilter_file, db_path[key]["faa"])
            
            key_prefilter_file.unlink()
            
            seq_id = ""
            ids_checked = []
                
            if result.returncode != 0:
                logging.error("An error has occured during seqkit process to "
                              f"obtain protein squences of {key}:"
                              f" {result.returncode}\n{result.stderr.decode('utf-8')}")
            else:
                fasta = result.stdout.decode("utf-8")
                ids_selected = [s.split()[0][1:] for s in fasta.split("\n") if s.startswith(">")]
                ids_checked.extend(ids_selected)
                logging.info("done")
            
            text += f"{len(ids_checked)} unique ids from {key} after filtering\n"
            
            for ic in ids_checked:
                    data_filtered_lines += data_map[ic]
                    sources_text += f"{ic} {key}\n"
        
        write_data(fasta_file,fasta,data_filtered_file,data_filtered_lines,
                   sources_file, sources_text, i)
    
    text += "\n"
    summary_file = outdir.parent / "summary.out"
    
    if summary_file.exists():
        with summary_file.open("a") as f:
            f.write(text)
    else:
        summary_file.write_text(text)
    
    logging.info(f"elapsed time: {datetime.datetime.now() - start}")