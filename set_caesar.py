import argparse
import difflib
import logging
import os
import psutil
import sys
import yaml
from pathlib import Path
from psutil._common import bytes2human


###############
## Functions ##
###############

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
            if yml[key] is None:
                logging.error(f"The '{key}' field in the config file: "
                              f"'{config_path}' seems to be empty")
                sys.exit(1)
            for db_format in yml[key]:
                db_path[db_name].update(db_format)
    
    return db_path, yml

def check_db_path(db_path, config_file):
    
    logging.info("Checking the validity of the database paths indicated in the "
                 "configuration file")
    
    if len(db_path) == 0:
        logging.error(f"No database has been found in your config file: "
                      f"'{config_file}'")
        sys.exit(1)
    
    for db_name in db_path:
        dmnd = False
        for db_type in db_path[db_name]:
            if db_type == "dmnd":
                dmnd = True
            if db_path[db_name][db_type] is None:
                logging.error(f"The key '{db_type}' of the '{db_name}_db' field "
                              f"in the config file: '{config_file}' is empty")
                sys.exit(1)
                
            db_file = Path(db_path[db_name][db_type])
            if not db_file.exists():
                logging.error(f"'{db_file}' was not found, please check your "
                              f"config file: '{config_file}'")
                sys.exit(1)

        if dmnd is False:
            logging.error(f"No diamond database provided for '{db_name}_db' field")
            sys.exit(1)
                
    logging.info("All paths are valid")
    
def check_config_options(yml, db_path):
    
    slurm = yml.get("slurm")
    parallel = yml.get("parallel")
    module = yml.get("module")
    sel_priority = yml.get("candidate_selection")
    
    if slurm == 1:
        logging.info("slurm option is set to '1': adapt the output "
                     "script to the slurm jobs scheduler")
        slurm = True
    else:
        logging.info("slurm option is set to '0': no adaption to the slurm"
                     " jobs scheduler")
        slurm = False
        
    if parallel == 1:
        logging.info("parallel option is set to '1': uses gnu parallel to run "
                     "two blastp in parallel")
        parallel = True
    else:
        logging.info("parallel option is set to '0': runs blastp in sequence")
        parallel = False
    
    if module is None:
        logging.info("no module required")
    elif type(module) == list:
        logging.info(f"The list of module to load: {module}")
    else:
        logging.error(f"The value: '{module}' of the key: 'module' is incorrect"
                      ", it must be a list of module to load")
        sys.exit(1)
    
    if sel_priority is None:
        logging.info()
        
    expected = ["strain_library", "order"]
    for db in db_path:
        if db not in ["swissprot", "trembl", "uniprot", "nr"]:
            expected.append(db)
    
    for elem in sel_priority:
        if elem not in expected:
            close = difflib.get_close_matches(elem.lower(), expected)
            if len(close) != 0:
                logging.warning(f"'{elem}' isn't expected in candidate_selection"
                                f" field. Perharps you mean one of these: {close}")
            else:
                logging.warning(f"'{elem}' isn't expected in candidate_selection"
                                f" field. Expected words are: {expected}")
                
    return module, slurm, parallel

def check_required_inputs(**kwargs):
    
    logging.info("Checking the required inputs")
    
    list_path = []
    
    for arg in kwargs:
        p = Path(kwargs[arg])
        
        if not p.exists():
            logging.error(f"'{p}' was not found")
            sys.exit(1)
        
        list_path.append(p)
    
    logging.info("All files are valid")
    
    return list_path

def check_general_options(slurm, threads, mem, outdir):
    
    if slurm is False:
        n = os.cpu_count()
        if threads > n:
            logging.error(f"-t, --threads value: '{threads}' is superior to the "
                          f"number of logical CPU cores: '{n}'")
            sys.exit(1)

        mem_tot = bytes2human(psutil.virtual_memory().total)
        if float(mem[:-1]) > float(mem_tot[:-1]):
            logging.error(f"-m, --mem value: '{mem}' is superior to the "
                          f"RAM of the system: '{mem_tot}'")
            sys.exit(1)
        

    outdir = Path(outdir).absolute()
    if not outdir.exists():
        outdir.mkdir()
        logging.info(f"'{outdir}' was created")
                
    return outdir

def set_blastp(slurm, parallel, args, db_path):
    
    # Path to the blastp.sh script
    main_path = Path(__file__).absolute()
    parent_path = main_path.parent
    src_path = parent_path / "src" / "blastp.sh"
    
    # Set the output directory of blastp.sh
    blastp_dir = Path(args.outdir).absolute() / "blastp"
    
    # Get all diamonad database paths
    dmnd = []
    for db in db_path:
        dmnd.append(db_path[db]["dmnd"])
        
    dmnd = " ".join(dmnd)
    
    text = "# Diamond blastp\n"
    if slurm == 1:
        pass
    else:
        text += f"bash {src_path} -t {args.threads} -q {args.query} -o {blastp_dir}"
        text += f" -i {args.blast_id} -c {args.blast_cov} "
        
        if parallel is True:
            text += f"-p {db}\n\n"
        else:
            text += db + "\n\n"
            
    return text

def check_blastp_options(pid, cov, min_len, max_len, tax):
    
    if pid < 0:
        logging.error(f"--blast-id must be greater than 0 but '{pid}' is given")
        sys.exit(1)
    elif pid < 1:
        logging.warning(f"--blast-id value: '{pid}' is ambiguous, multiply it "
                        "by 100 if you don't want a percentage less than 1%")
    elif pid > 100:
        logging.error(f"--blast-id must be lower than 100 but '{pid}' is given")
        sys.exit(1)
    
    if cov < 0:
        logging.error(f"--blast-cov must be greater than 0 but '{cov}' is given")
        sys.exit(1)
    elif cov < 1:
        logging.warning(f"--blast-cov value: '{cov}' is ambiguous, multiply it "
                        "by 100 if you don't want a percentage less than 1%")
    elif cov > 100:
        logging.error(f"--blast-cov must be lower than 100 but '{cov}' is given")
        sys.exit(1)
        
    if min_len < 0:
        logging.error(f"--min-len must be a positive number but '{min_len}' "
                      "is given")
        sys.exit(1)
    elif min_len > max_len:
        min_len, max_len = max_len, min_len
        logging.warning(f"--min-len value: '{min_len}' is greater than --max-len"
                        f" value: '{max_len}', so the values are swapped")
    if max_len < 0:
        logging.error(f"--max-len must be a positive number but '{max_len}' "
                      "is given")
        sys.exit(1)
        
    if "".join(sorted(tax)) not in ["A", "B", "E", "AB", "AE", "BE", "ABE"]:
        logging.error(f"--tax value contains invalid characters: '{tax}', the"
                      f" allowed characters are, A, B and E")

def set_filter(slurm, args):
    
    # Path to the filter.py script
    main_path = Path(__file__).absolute()
    parent_path = main_path.parent
    src_path = parent_path / "src" / "filter.py"
    
    # Set the output directory of filter.py and blastp.sh
    blastp_dir = Path(args.outdir).absolute() / "blastp"
    filtered_dir = Path(args.outdir).absolute() / "filtered"
    
    # Checks the blastp options
    pid = args.blast_id
    cov = args.blast_cov
    min_len = args.min_len
    max_len = args.max_len
    tax = args.tax
    
    check_blastp_options(pid, cov, min_len, max_len, tax)
    
    text = "# Filter\n"
    
    if slurm is True:
        pass
    else:
        text += f"python {src_path} -o {filtered_dir} -c {args.config} "
        text += f"-q {args.query} -d {blastp_dir} --id {pid} --cov "
        text += f"{cov} --min {min_len} --max {max_len} --tax {tax}\n\n"

    return text

def set_clustering(slurm, args):
    
    # Path to the clustering.sh script
    main_path = Path(__file__).absolute()
    parent_path = main_path.parent
    src_path = parent_path / "src" / "clustering.sh"
    
    # Set the directory paths
    blastp_dir = Path(args.outdir).absolute() / "blastp"
    filtered_dir = Path(args.outdir).absolute() / "filtered"
    clusters_dir = Path(args.outdir).absolute() / "clusters"
    
    # Checks the blastp options
    pid = args.cluster_id
    cov = args.cluster_cov
    
    check_clustering_options(pid, cov)
    
    text = "# Diamond clustering\n"
    
    if slurm is True:
        pass
    else:
        text += f"bash {src_path} -o {clusters_dir} -f {filtered_dir}/filtered_"
        text += f"sequences.fasta -i {pid} -c {cov} -t {args.threads} -m "
        text += f"{args.mem}\n\n"
        
    return text
    
def check_clustering_options(pid, cov):
    
    if pid < 0:
        logging.error(f"--cluster-id must be greater than 0 but '{pid}' is given")
        sys.exit(1)
    elif pid < 1:
        logging.warning(f"--cluster-id value: '{pid}' is ambiguous, multiply it "
                        "by 100 if you don't want a percentage less than 1%")
    elif pid > 100:
        logging.error(f"--cluster-id must be lower than 100 but '{pid}' is given")
        sys.exit(1)
    
    if cov < 0:
        logging.error(f"--cluster-cov must be greater than 0 but '{cov}' is given")
        sys.exit(1)
    elif cov < 1:
        logging.warning(f"--cluster-cov value: '{cov}' is ambiguous, multiply it "
                        "by 100 if you don't want a percentage less than 1%")
    elif cov > 100:
        logging.error(f"--cluster-cov must be lower than 100 but '{cov}' is given")
        sys.exit(1)

##########
## MAIN ##
##########

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s - %(asctime)s - %(message)s ")
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", '--outdir', type=str, default=Path.cwd(), metavar="",
                        help="Output directory [default: './']")
    parser.add_argument("-s", "--start", type=str, metavar="", default="blastp",
                        choices=["blastp", "filter", "clustering",
                                 "selection"],
                        help="Selects the inital step: 'blastp', 'filter',"
                        " 'clustering' or 'selection' [default: 'blastp']")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=4,
                        help="Number of CPU threads for diamond processes"
                        " [default: 4]")
    parser.add_argument("-m", "--mem", type=str, default="4G", metavar="",
                        help="Memory limit for diamond process [default: '4G']")
    
    required_input = parser.add_argument_group("Mandatory inputs")
    required_input.add_argument("-c", "--config", type=str, metavar="",
                                required=True, help="Configuration file")
    required_input.add_argument("-q", "--query", type=str, metavar="",
                                required=True, help="Set of reference sequences")
    
    blast_opt = parser.add_argument_group("Blastp options")
    blast_opt.add_argument("--blast-id", default=30.0, type=float, metavar="",
                           help="Retains only candidates above the specified"
                           " percentage of sequence identity [default: 30.0]")
    blast_opt.add_argument("--blast-cov", default=80.0, type=float, metavar="",
                           help="Retains only candidates above the specified"
                           " percentage of query cover [default: 80.0]")
    blast_opt.add_argument("--min-len", default=200, type=int, metavar="",
                           help="Retains only candidates above the specified"
                           " sequence length [default: 200]")
    blast_opt.add_argument("--max-len", default=1000, type=int, metavar="",
                           help="Retains only candidates below the specified"
                           " sequence length [default: 1000]")
    blast_opt.add_argument("--tax", default="ABE", type=str, metavar="",
                           help="Superkingdom filter, A: Archaea, B: Bacteria"
                           " and E: Eukaryota [default: 'ABE']")
    
    clustering_opt = parser.add_argument_group("Clustering options")
    clustering_opt.add_argument("--cluster-id", default=80.0, type=float,
                                metavar="",
                                help="Identity cutoff for the clustering "
                                "[default: 80.0]")
    clustering_opt.add_argument("--cluster-cov", default=80.0, type=float,
                                metavar="",
                                help="Minimum coverage of cluster member sequence")
    
    selection_opt = parser.add_argument_group("Candidates selection options")
    selection_opt.add_argument("--gc", default=50.0, type=float, metavar="",
                               help="Target GC percentage to decide between"
                               " candidates [default: 50.0]")
    ncand_group = selection_opt.add_mutually_exclusive_group()
    ncand_group.add_argument("-n", "--nb-cand", type=int, metavar="", default=1,
                        help="maximum number of candidate per cluster [default: 1]")
    ncand_group.add_argument("--cov-per-cluster", type=float, metavar="",
                             help="Uses a percentage of each cluster as "
                             "maximum number of candidates rather than a given"
                             " number")
    
    data_opt = parser.add_argument_group("Required if --start equals to"
                                         " 'filter', 'clustering' or 'selection'",
                                         "Could be a directory containing "
                                         "blastp tsv outputs if 'filter'")
    data_opt.add_argument("-d", "--data", type=str, metavar="",
                          help="Directory of blastp outputs, a blastp output "
                          "or a filtered blastp ourput")
    
    fasta_source = parser.add_argument_group("Required if --start equals to"
                                             " 'clustering' or 'selection'")
    fasta_source.add_argument("-f", "--fasta-cand", type=str, metavar="",
                              help="Multi fasta of protein sequences")
    fasta_source.add_argument("--sources", type=str, metavar="",
                              help="Sources file indicating the sources"
                              " database of each sequences")
    
    exlusion_opt = parser.add_argument_group("Exclude some protein ids",
                                            "Can be used to re-run pipeline"
                                            " and try to get other candidates")
    exlusion_opt.add_argument("-u", "--update", type=str, metavar="",
                              help="File containing the list"
                              " of proteins ids not to be selected as candidates")
    
    args = parser.parse_args()
    
    config_file, query_file = check_required_inputs(config=args.config,
                                                    query=args.query)
    
    db_path,yml = read_yaml_config(Path(args.config))
    check_db_path(db_path, args.config)
    module, slurm, parallel = check_config_options(yml, db_path)
    outdir = check_general_options(slurm, threads=args.threads, mem=args.mem,
                                   outdir=args.outdir)
    
    caesar_text = "#!/bin/bash\n\n"
    if module is not None:
        for elem in module:
            caesar_text += f"module load {elem}\n"
        
        caesar_text += "\n"
        
    if args.start == "blastp":
        caesar_text += set_blastp(slurm, parallel, args, db_path)
        
    if args.start in ["blastp", "filter"]:
        caesar_text += set_filter(slurm, args)
        
    if args.start in  ["blastp", "filter", "clustering"]:
        caesar_text += set_clustering(slurm, args)
    
    print(caesar_text)