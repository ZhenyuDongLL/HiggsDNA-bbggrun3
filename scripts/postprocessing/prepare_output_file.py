#!/usr/bin/env python
# Author Tiziano Bevilacqua (03/03/2023) and Nico Haerringer (16/07/2024)
import os
import subprocess
from optparse import OptionParser
import json, glob
from higgs_dna.utils.logger_utils import setup_logger

import logging
from concurrent.futures import ThreadPoolExecutor

# ---------------------- A few helping functions  ----------------------


def MKDIRP(dirpath, verbose=False, dry_run=False):
    if verbose:
        print("\033[1m" + ">" + "\033[0m" + ' os.mkdirs("' + dirpath + '")')
    if dry_run:
        return
    try:
        os.makedirs(dirpath)
    except OSError:
        if not os.path.isdir(dirpath):
            raise
    return


# function to activate the (pre-existing) FlashggFinalFit environment and use Tree2WS scripts
# full Power Ranger style :D
def activate_final_fit(path, command):
    current_path = os.getcwd()
    os.chdir(path)
    os.system(
        f"eval `scram runtime -sh` && source {path}/flashggFinalFit/setup.sh && cd {path}/flashggFinalFit/Trees2WS && {command} "
    )
    os.chdir(current_path)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# - EXAMPLE USAGE: ----------------------------------------------------------------------------------------------------------------------------------------------------------------#
# - python3 prepare_output_file.py --input <dir_to_HiggsDNA_dump> --merge --varDict <path_to_varDict> --root --syst --cats --catDict <path_to_catDict> --output <path_to_output_dir>
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Read options from command line
usage = "Usage: python %prog filelists [options]"
parser = OptionParser(usage=usage)
parser.add_option("--input", dest="input", type="string", default="", help="input dir")
parser.add_option(
    "--merge",
    dest="merge",
    action="store_true",
    default=False,
    help="Do merging of the .parquet files",
)
parser.add_option(
    "--varDict",
    dest="varDict",
    default=None,
    help="Path to JSON that holds dictionary that encodes the mapping of the systematic variation branches (includes nominal and object-based systematics, up/down). If not provided, use only nominal.",
)
parser.add_option(
    "--root",
    dest="root",
    action="store_true",
    default=False,
    help="Do root conversion step",
)
parser.add_option(
    "--ws",
    dest="ws",
    action="store_true",
    default=False,
    help="Do root to workspace conversion step",
)
parser.add_option(
    "--ws_config",
    dest="config",
    type="string",
    default="config_simple.py",
    help="configuration file for Tree2WS, as it is now it must be stored in Tree2WS directory in FinalFit",
)
parser.add_option(
    "--final_fit",
    dest="final_fit",
    type="string",
    default="/afs/cern.ch/user/n/niharrin/cernbox/PhD/Higgs/CMSSW_10_2_13/src/",
    help="FlashggFinalFit path",
)  # the default is just for me, it should be changed but I don't see a way to make this generally valid
parser.add_option(
    "--syst",
    dest="syst",
    action="store_true",
    default=False,
    help="Do systematics variation treatment",
)
parser.add_option(
    "--cats",
    dest="cats",
    action="store_true",
    default=False,
    help="Split into categories.",
)
parser.add_option(
    "--catDict",
    dest="catDict",
    default=None,
    help="Path to JSON that defines the conditions for splitting into multiple categories. For well-defined statistical analyses in final fits, the categories should be mutually exclusive (it is your job to ensure this!). If not provided, use only one inclusive untagged category with no conditions.",
)
parser.add_option(
    "--genBinning",
    dest="genBinning",
    default="",
    help="Optional: Path to the JSON containing the binning at gen-level.",
)
parser.add_option(
    "--skip-normalisation",
    dest="skip_normalisation",
    action="store_true",
    default=False,
    help="Independent of file type, skip normalisation step",
)
parser.add_option(
    "--args",
    dest="args",
    type="string",
    default="",
    help="additional options for root converter: --notag",
)
parser.add_option(
    "--verbose",
    dest="verbose",
    type="string",
    default="INFO",
    help="Verbosity level for the logger: INFO (default), DEBUG",
)
parser.add_option(
    "--output",
    dest="output",
    type="string",
    default="",
    help="Output path for the merged and ROOT files.",
)
parser.add_option(
    "--folder-structure",
    dest="folder_structure",
    type="string",
    default="",
    help="Uses the given folder structure for dirlist.",
)
parser.add_option(
    "--merge-data-only",
    dest="merge_data",
    action="store_true",
    default=False,
    help="Flag for merging data to an allData file.",
)
parser.add_option(
    "--make-condor-logs",
    dest="make_condor_logs",
    action="store_true",
    default=False,
    help="Condor log files activated.",
)
parser.add_option(
    "--condor-logs",
    dest="condor_logs",
    type="string",
    default="",
    help="Output path of the Condor log files.",
)
parser.add_option(
    "--apptainer",
    dest="apptainer",
    action="store_true",
    default=False,
    help="Run HTCondor with Docker image of HiggsDNA's current master branch.",
)
(opt, args) = parser.parse_args()

if (opt.verbose != "INFO") and (opt.verbose != "DEBUG"):
    opt.verbose = "INFO"
logger = setup_logger(level=opt.verbose)

folder_for_dirlist = opt.input
if opt.apptainer:
    if opt.root and not opt.merge:
        folder_for_dirlist = opt.input + "/merged"
    elif opt.folder_structure != "":
        folder_for_dirlist = opt.folder_structure
else:
    if opt.root and not opt.merge:
        folder_for_dirlist = opt.input + "/merged"
    if opt.folder_structure != "":
        folder_for_dirlist = opt.folder_structure
# os.system(
#    f"ls -l {folder_for_dirlist} | tail -n +2 | grep -v .coffea | grep -v merged | grep -v root |"
#     + "awk '{print $NF}' > dirlist.txt"
# )

os.system(
    f"find {folder_for_dirlist} -mindepth 1 -maxdepth 1 -type d | grep -v '^.$' | grep -v .coffea | grep -v '/merged$' | grep -v '/root$' |"
    + "awk -F'/' '{print $NF}' > dirlist.txt"
    )

process_dict = {
    "GluGluHtoGG": "ggh",
    "VBFHtoGG": "vbf",
    "VHtoGG": "vh",
    "ttHtoGG": "tth",
    "GluGluHtoGG_M-120_preEE": "ggh_120",
    "GluGluHtoGG_M-120_postEE": "ggh_120",
    "GluGluHtoGG_M-125_preEE": "ggh_125",
    "GluGluHtoGG_M-125_postEE": "ggh_125",
    "GluGluHtoGG_M-130_preEE": "ggh_130",
    "GluGluHtoGG_M-130_postEE": "ggh_130",
    "VBFHtoGG_M-120_preEE": "vbf_120",
    "VBFHtoGG_M-120_postEE": "vbf_120",
    "VBFHtoGG_M-125_preEE": "vbf_125",
    "VBFHtoGG_M-125_postEE": "vbf_125",
    "VBFHtoGG_M-130_preEE": "vbf_130",
    "VBFHtoGG_M-130_postEE": "vbf_130",
    "VHtoGG_M-120_preEE": "vh_120",
    "VHtoGG_M-120_postEE": "vh_120",
    "VHtoGG_M-125_preEE": "vh_125",
    "VHtoGG_M-125_postEE": "vh_125",
    "VHtoGG_M-130_preEE": "vh_130",
    "VHtoGG_M-130_postEE": "vh_130",
    "ttHtoGG_M-120_preEE": "tth_120",
    "ttHtoGG_M-120_postEE": "tth_120",
    "ttHtoGG_M-125_preEE": "tth_125",
    "ttHtoGG_M-125_postEE": "tth_125",
    "ttHtoGG_M-130_preEE": "tth_130",
    "ttHtoGG_M-130_postEE": "tth_130",
    # Shorter conventions
    "ggh_M-120_preEE": "ggh_120",
    "ggh_M-120_postEE": "ggh_120",
    "ggh_M-125_preEE": "ggh_125",
    "ggh_M-125_postEE": "ggh_125",
    "ggh_M-130_preEE": "ggh_130",
    "ggh_M-130_postEE": "ggh_130",
    "vbf_M-120_preEE": "vbf_120",
    "vbf_M-120_postEE": "vbf_120",
    "vbf_M-125_preEE": "vbf_125",
    "vbf_M-125_postEE": "vbf_125",
    "vbf_M-130_preEE": "vbf_130",
    "vbf_M-130_postEE": "vbf_130",
    "vh_M-120_preEE": "vh_120",
    "vh_M-120_postEE": "vh_120",
    "vh_M-125_preEE": "vh_125",
    "vh_M-125_postEE": "vh_125",
    "vh_M-130_preEE": "vh_130",
    "vh_M-130_postEE": "vh_130",
    "tth_M-120_preEE": "tth_120",
    "tth_M-120_postEE": "tth_120",
    "tth_M-125_preEE": "tth_125",
    "tth_M-125_postEE": "tth_125",
    "tth_M-130_preEE": "tth_130",
    "tth_M-130_postEE": "tth_130",
    "DYto2L_2Jets": "dy",
    "GG-Box-3Jets_MGG-80_postEE": "ggbox",
    "GG-Box-3Jets_MGG-80_preEE": "ggbox",
    "GJet_PT-20to40_DoubleEMEnriched_MGG-80_postEE": "gjet",
    "GJet_PT-20to40_DoubleEMEnriched_MGG-80_preEE": "gjet",
    "GJet_PT-40_DoubleEMEnriched_MGG-80_postEE": "gjet",
    "GJet_PT-40_DoubleEMEnriched_MGG-80_preEE": "gjet"
}

# the key of the var_dict entries is also used as a key for the related root tree branch
# to be consistent with FinalFit naming scheme you shoud use SystNameUp and SystNameDown, 
# e.g. "FNUFUp": "FNUF_up", "FNUFDown": "FNUF_down" 
if opt.varDict is None: # If not given in the command line
    logger.info("You did not specify the path to a variation dictionary JSON, so we will only use the nominal input trees.")
    var_dict = {
        "NOMINAL": "nominal",
    }
else:
    with open(opt.varDict, "r") as jf:
        var_dict = json.load(jf)['var_dict']
    

# Here we prepare to split the output into categories, in the dictionary are defined the cuts to be applyed by pyarrow.ParquetDataset
# when reading the data, the variable obviously has to be in the dumped .parquet
# This can be improved by passing the configuration via json loading
if opt.cats and opt.catDict is not None:
    with open(opt.catDict, "r") as jf:
        cat_dict = json.load(jf)['cat_dict']
else:
    logger.info("You chose to run without cats or you did not specify the path to a categorisation dictionary JSON, so we will only use one inclusive NOTAG category.")
    cat_dict = {"NOTAG": {"cat_filter": [("pt", ">", -1.0)]}}

# Now, after loading the JSONs from possibly relative paths, we can change the directory appropriately to get to work
EXEC_PATH = os.path.realpath(os.getcwd())
os.chdir(opt.input)
IN_PATH = os.path.realpath(os.getcwd())
SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)  # script directory

CONDOR_PATH = os.path.realpath(os.path.abspath(opt.condor_logs)) # need real absolute path. otherwise problems can arise on lxplus (between afs and eos)
# I create a dictionary and save it to a temporary json so that this can be shared between the two scripts
# and then gets deleted to not leave trash around. We have to care for the environment :P.
# Not super elegant, open for suggestions
# with open("category.json", "w") as file:
#     file.write(json.dumps(cat_dict))

# Same for variation dictionary, which is shared between merging and ROOTing steps.
# with open("variation.json", "w") as file:
#     file.write(json.dumps(var_dict))

# Using OUT_PATH for the location of the output if different from the input path
if not opt.apptainer:
    if opt.output == "":
        OUT_PATH = IN_PATH
        # os.system(f"mv category.json {SCRIPT_DIR}/../../higgs_dna/category.json")
        # os.system(f"mv variation.json {SCRIPT_DIR}/../../higgs_dna/variation.json")
        with open(f"{SCRIPT_DIR}/../../higgs_dna/category.json", "w") as file:
            file.write(json.dumps(cat_dict))
        with open(f"{SCRIPT_DIR}/../../higgs_dna/variation.json", "w") as file:
            file.write(json.dumps(var_dict))
        # if opt.folder_structure != "":
        #     dirlist_path = folder_for_dirlist+"/dirlist.txt"
        # else:
        dirlist_path = f"{EXEC_PATH}/dirlist.txt"
    else:
        OUT_PATH = opt.output
        # os.system(f"mv category.json {OUT_PATH}/category.json")
        # os.system(f"mv variation.json {OUT_PATH}/variation.json")
        # os.system(f"mv category.json {SCRIPT_DIR}/../../higgs_dna/category.json")
        # os.system(f"mv variation.json {SCRIPT_DIR}/../../higgs_dna/variation.json")
        with open(f"{SCRIPT_DIR}/../../higgs_dna/category.json", "w") as file:
            file.write(json.dumps(cat_dict))
        with open(f"{SCRIPT_DIR}/../../higgs_dna/variation.json", "w") as file:
            file.write(json.dumps(var_dict))
        # if opt.folder_structure != "":
        #     dirlist_path = folder_for_dirlist+"/dirlist.txt"
        #     os.system(f"mv {dirlist_path} {OUT_PATH}/dirlist.txt")
        # else:
        dirlist_path = f"{OUT_PATH}/dirlist.txt"
        os.system(f"mv {EXEC_PATH}/dirlist.txt {OUT_PATH}/dirlist.txt")
    cat_dict = "category.json"
else:
    if opt.output == "":
        OUT_PATH = IN_PATH
        cat_dict_loc = f"{SCRIPT_DIR}/../../higgs_dna/category.json"
        var_dict_loc = f"{SCRIPT_DIR}/../../higgs_dna/variation.json"
        # os.system(f"mv category.json {cat_dict_loc}")
        # os.system(f"mv variation.json {var_dict_loc}")
        with open(cat_dict_loc, "w") as file:
            file.write(json.dumps(cat_dict))
        with open(var_dict_loc, "w") as file:
            file.write(json.dumps(var_dict))

        dirlist_path = f"{EXEC_PATH}/dirlist.txt"
    else:
        OUT_PATH = os.path.realpath(opt.output)
        cat_dict_loc = f"{OUT_PATH}/category.json"
        var_dict_loc = f"{OUT_PATH}/variation.json"
        # os.system(f"mv category.json {cat_dict_loc}")
        # os.system(f"mv variation.json {var_dict_loc}")
        with open(cat_dict_loc, "w") as file:
            file.write(json.dumps(cat_dict))
        with open(var_dict_loc, "w") as file:
            file.write(json.dumps(var_dict))

        dirlist_path = f"{OUT_PATH}/dirlist.txt"
        os.system(f"mv {EXEC_PATH}/dirlist.txt {OUT_PATH}/dirlist.txt")
    
def submit_jobs(directory, suffix=""):
    if suffix != "":
        sub_files = glob.glob(f"{directory}/*{suffix}.sub")
    else:
        sub_files = glob.glob(f"{directory}/*.sub")
    for current_file in sub_files:
        subprocess.run(["condor_submit", "-spool", current_file])

if opt.genBinning != "":
    genBinning_str = f"--genBinning {opt.genBinning}"
else:
    genBinning_str = ""
# Define string if normalisation to be skipped
skip_normalisation_str = "--skip-normalisation" if opt.skip_normalisation else ""

# The process var below is the function that will be executed in parallel for each systematic variation. It substitutes the old loop of the systematics to speed up the process.
# Paths now must be ABSOLUTE!! - CD while multi thread is not a good idea! 
def process_var(var, var_dict, IN_PATH, OUT_PATH, SCRIPT_DIR, file, cat_dict, skip_normalisation_str):
    target_dir = f"{OUT_PATH}/merged/{file}/{var_dict[var]}"
    MKDIRP(target_dir)

    command = f"python3 merge_parquet.py --source {IN_PATH}/{file}/{var_dict[var]} --target {target_dir}/ --cats {cat_dict} {skip_normalisation_str} {genBinning_str}"
    logger.info(command)

    # Execute the command using subprocess.run
    subprocess.run(command, shell=True, cwd=SCRIPT_DIR, check=True)

# Loop to paralelize the loop over the "files", which are the ttH_125_preEE, etc. datasets
def process_file(file, IN_PATH, OUT_PATH, SCRIPT_DIR, var_dict, cat_dict, skip_normalisation_str, opt):
    file = file.strip()  # Removes newline characters and leading/trailing whitespace
    if "data" not in file.lower():
        target_path = f"{OUT_PATH}/merged/{file}"
        if os.path.exists(target_path):
            raise Exception(f"The selected target path: {target_path} already exists")
        MKDIRP(target_path)

        if opt.syst:
            # Systematic variations processing
            with ThreadPoolExecutor(max_workers=7) as executor:
                futures = [executor.submit(process_var, var, var_dict, IN_PATH, OUT_PATH, SCRIPT_DIR, file, cat_dict, skip_normalisation_str) for var in var_dict]

            for future in futures:
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Error processing variable: {e}")
        else:
            # Single nominal processing for MC
            command = f"python3 merge_parquet.py --source {IN_PATH}/{file}/nominal --target {target_path}/ --cats {cat_dict} {skip_normalisation_str} {genBinning_str}"
            subprocess.run(command, shell=True, cwd=SCRIPT_DIR, check=True)
    else:
        # Data processing
        merged_target_path = f"{OUT_PATH}/merged/{file}/{file}_merged.parquet"
        data_dir_path = f'{OUT_PATH}/merged/Data_{file.split("_")[-1]}'
        if os.path.exists(merged_target_path):
            raise Exception(f"The selected target path: {merged_target_path} already exists")
        if not os.path.exists(data_dir_path):
            MKDIRP(data_dir_path)
        command = f'python3 merge_parquet.py --source {IN_PATH}/{file}/nominal --target {data_dir_path}/{file}_ --cats {cat_dict} --is-data {genBinning_str}'
        subprocess.run(command, shell=True, cwd=SCRIPT_DIR, check=True)

if not opt.apptainer:
    if opt.merge:
        with open(dirlist_path) as fl:
            files = fl.readlines()
            
            # No more loop over the files, we will use the ThreadPoolExecutor to parallelize the process!
            with ThreadPoolExecutor(max_workers=8) as executor:
                futures = [executor.submit(process_file, file, IN_PATH, OUT_PATH, SCRIPT_DIR, var_dict, cat_dict, skip_normalisation_str, opt) for file in files]

            # Optionally, wait for all futures to complete and check for exceptions
            for future in futures:
                try:
                    future.result()
                except Exception as e:
                    # Log file-level exceptions
                    logger.error(f"Error processing file: {e}")


            # at this point Data will be split in eras if any Data dataset is present, here we merge them again in one allData file to rule them all
            # we also skip this step if there is no Data
            for file in files:
                file = file.split("\n")[0]  # otherwise it contains an end of line and messes up the os.walk() call
                if "data" in file.lower() or "DoubleEG" in file:
                    dirpath, dirnames, filenames = next(os.walk(f'{OUT_PATH}/merged/Data_{file.split("_")[-1]}'))
                    if len(filenames) > 0:
                        command = f'python3 merge_parquet.py --source {OUT_PATH}/merged/Data_{file.split("_")[-1]} --target {OUT_PATH}/merged/Data_{file.split("_")[-1]}/allData_ --cats {cat_dict} --is-data {genBinning_str}'
                        subprocess.run(command, shell=True, cwd=SCRIPT_DIR, check=True)
                        break
                    else:
                        logger.info(f'No merged parquet found for {file} in the directory: {OUT_PATH}/merged/Data_{file.split("_")[-1]}')

    if opt.root:
        logger.info("Starting root step")
        if opt.syst:
            logger.info("you've selected the run with systematics")
            args = "--do_syst"
        else:
            logger.info("you've selected the run without systematics")
            args = ""

        if opt.merge:
            IN_PATH = OUT_PATH
        # Note, in my version of HiggsDNA I run the analysis splitting data per Era in different datasets
        # the treatment of data here is tested just with that structure
        with open(dirlist_path) as fl:
            files = fl.readlines()
            for file in files:
                file = file.split("\n")[0]
                if "data" not in file.lower() and file in process_dict:
                    if os.path.exists(f"{OUT_PATH}/root/{file}"):
                        raise Exception(
                            f"The selected target path: {OUT_PATH}/root/{file} already exists"
                        )

                    if os.listdir(f"{IN_PATH}/merged/{file}/"):
                        logger.info(f"Found merged files {IN_PATH}/merged/{file}/")
                    else:
                        raise Exception(f"Merged parquet not found at {IN_PATH}/merged/")
                    MKDIRP(f"{OUT_PATH}/root/{file}")
                    os.chdir(SCRIPT_DIR)
                    os.system(
                        f"python3 convert_parquet_to_root.py {IN_PATH}/merged/{file}/merged.parquet {OUT_PATH}/root/{file}/merged.root mc --process {process_dict[file]} {args} --cats {cat_dict} --vars variation.json {genBinning_str}"
                    )
                elif "data" in file.lower():
                    if os.listdir(f'{IN_PATH}/merged/Data_{file.split("_")[-1]}/'):
                        logger.info(
                            f'Found merged data files in: {IN_PATH}/merged/Data_{file.split("_")[-1]}/'
                        )
                    else:
                        raise Exception(
                            f'Merged parquet not found at: {IN_PATH}/merged/Data_{file.split("_")[-1]}/'
                        )

                    if os.path.exists(
                        f'{OUT_PATH}/root/Data/allData_{file.split("_")[-1]}.root'
                    ):
                        logger.info(
                            f'Data already converted: {OUT_PATH}/root/Data/allData_{file.split("_")[-1]}.root'
                        )
                        continue
                    elif not os.path.exists(f"{OUT_PATH}/root/Data/"):
                        MKDIRP(f"{OUT_PATH}/root/Data")
                        os.chdir(SCRIPT_DIR)
                        os.system(
                            f'python3 convert_parquet_to_root.py {IN_PATH}/merged/Data_{file.split("_")[-1]}/allData_merged.parquet {OUT_PATH}/root/Data/allData_{file.split("_")[-1]}.root data --cats {cat_dict} --vars variation.json {genBinning_str}'
                        )
                    else:
                        os.chdir(SCRIPT_DIR)
                        os.system(
                            f'python3 convert_parquet_to_root.py {IN_PATH}/merged/Data_{file.split("_")[-1]}/allData_merged.parquet {OUT_PATH}/root/Data/allData_{file.split("_")[-1]}.root data --cats {cat_dict} --vars variation.json {genBinning_str}'
                        )

    if opt.ws:
        if not os.listdir(opt.final_fit):
            raise Exception(
                f"The selected FlashggFinalFit path: {opt.final_fit} is invalid"
            )
            
        if os.path.exists(f"{IN_PATH}/root/Data"):
            os.system(f"echo Data >> {dirlist_path}")

        data_done = False

        with open(dirlist_path) as fl:
            files = fl.readlines()
            if opt.syst:
                doSystematics = "--doSystematics"
            else:
                doSystematics = ""
            with open(f"{SCRIPT_DIR}/../../higgs_dna/category.json") as f:
                cat_file = json.load(f)
            for dir in files:
                dir = dir.split("\n")[0]
                # if MC
                if "data" not in dir.lower() and dir in process_dict:
                    if os.listdir(f"{IN_PATH}/root/{dir}/"):
                        filename = subprocess.check_output(
                            f"find {IN_PATH}/root/{dir} -name *.root -type f",
                            shell=True,
                            universal_newlines=True,
                        )
                    else:
                        raise Exception(
                            f"The selected target path: {IN_PATH}/root/{dir} it's empty"
                        )
                    doNOTAG = ""
                    if ("NOTAG" in cat_file.keys()):
                        doNOTAG = "--doNOTAG"
                    command = f"python trees2ws.py {doNOTAG} --inputConfig {opt.config} --productionMode {process_dict[dir]} --year 2017 {doSystematics} --inputTreeFile {filename}"
                    activate_final_fit(opt.final_fit, command)
                elif "data" in dir.lower() and not data_done:
                    if os.listdir(f"{IN_PATH}/root/Data/"):
                        filename = subprocess.check_output(
                            f"find {IN_PATH}/root/Data -name *.root -type f",
                            shell=True,
                            universal_newlines=True,
                        )
                    else:
                        raise Exception(
                            f"The selected target path: {IN_PATH}/root/{dir} it's empty"
                        )
                    doNOTAG = ""
                    if ("NOTAG" in cat_file.keys()):
                        doNOTAG = "--doNOTAG"
                    command = f"python trees2ws_data.py {doNOTAG} --inputConfig {opt.config} --inputTreeFile {filename}"
                    activate_final_fit(opt.final_fit, command)
                    data_done = True
        os.chdir(EXEC_PATH)

else:
    if opt.merge:
        with open(dirlist_path) as fl:
            files = fl.readlines()
            if not opt.merge_data:
                for file in files:
                    file = file.split("\n")[0]
                    # parent_id = 0
                    # MC dataset are identified as everythingthat does not contain "data" or "Data" in the name.
                    if "data" not in file.lower():
                        if opt.condor_logs != "":
                            job_file_executable = os.path.join(CONDOR_PATH, f"{file}.sh")
                            job_file_submit = os.path.join(CONDOR_PATH, f"{file}.sub")
                        else:
                            job_file_executable = os.path.join(OUT_PATH, f"{file}.sh")
                            job_file_submit = os.path.join(OUT_PATH, f"{file}.sub")

                        if not opt.make_condor_logs:
                            job_file_out = "/dev/null"
                            job_file_err = "/dev/null"
                            job_file_log = "/dev/null"
                        elif opt.condor_logs != "":
                            job_file_out = os.path.join(CONDOR_PATH, f"{file}.$(ClusterId).$(ProcId).out")
                            job_file_err = os.path.join(CONDOR_PATH, f"{file}.$(ClusterId).$(ProcId).err")
                            job_file_log = os.path.join(CONDOR_PATH, f"{file}.$(ClusterId).log")
                        else:
                            job_file_out = os.path.join(OUT_PATH, f"{file}.$(ClusterId).$(ProcId).out")
                            job_file_err = os.path.join(OUT_PATH, f"{file}.$(ClusterId).$(ProcId).err")
                            job_file_log = os.path.join(OUT_PATH, f"{file}.$(ClusterId).log")
                        
                        with open(job_file_executable, "w") as executable_file:
                            executable_file.write("#!/bin/sh\n")
                            if os.path.exists(f"{OUT_PATH}/merged/{file}"):
                                raise Exception(
                                    f"The selected target path: {OUT_PATH}/merged/{file} already exists"
                                )

                            MKDIRP(f"{OUT_PATH}/merged/{file}")
                            if opt.syst:
                                # if we have systematic variations in different files we have to split them in different directories
                                # otherwise they will be all merged at once in the same output file
                                i = 0
                                for var in var_dict:
                                    os.chdir(OUT_PATH)

                                    MKDIRP(f"{OUT_PATH}/merged/{file}/{var_dict[var]}")

                                    os.chdir(SCRIPT_DIR)
                                    logger.info(f"merge_parquet.py --source {IN_PATH}/{file}/{var_dict[var]} --target {OUT_PATH}/merged/{file}/{var_dict[var]}/ --cats {cat_dict_loc} {skip_normalisation_str} --abs {genBinning_str}")
                                    executable_file.write(f"if [ $1 -eq {i} ]; then\n")
                                    executable_file.write(f"    merge_parquet.py --source {IN_PATH}/{file}/{var_dict[var]} --target {OUT_PATH}/merged/{file}/{var_dict[var]}/ --cats {cat_dict_loc} {skip_normalisation_str} --abs {genBinning_str} || exit 107\n")
                                    executable_file.write("exit 0\n")
                                    executable_file.write("fi\n")
                                    i += 1

                            else:
                                i = 1
                                os.chdir(SCRIPT_DIR)
                                print(f"merge_parquet.py --source {IN_PATH}/{file}/nominal --target {OUT_PATH}/merged/{file}/ --cats {cat_dict_loc} {skip_normalisation_str} --abs {genBinning_str}")
                                executable_file.write(f"if [ $1 -eq 0 ]; then\n")
                                executable_file.write(f"    merge_parquet.py --source {IN_PATH}/{file}/nominal --target {OUT_PATH}/merged/{file}/ --cats {cat_dict_loc} {skip_normalisation_str} --abs {genBinning_str} || exit 107\n")
                                executable_file.write("exit 0\n")
                                executable_file.write("fi\n")
                                
                        os.system(f"chmod 775 {job_file_executable}")
                        with open(job_file_submit, "w") as submit_file:
                            if opt.condor_logs != "": submit_file.write(f"initialdir = {CONDOR_PATH}\n")
                            submit_file.write(f"executable = {job_file_executable}\n")
                            submit_file.write("arguments = $(ProcId)\n")
                            submit_file.write(f"output = {job_file_out}\n")
                            submit_file.write(f"error = {job_file_err}\n")
                            submit_file.write(f"log = {job_file_log}\n")
                            submit_file.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
                            submit_file.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
                            # if opt.max_materialize != "": submit_file.write(f"max_materialize = {opt.max_materialize}\n")
                            if opt.apptainer:
                                submit_file.write("MY.XRDCP_CREATE_DIR     = True\n")
                                submit_file.write("""MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/higgsdna:lxplus-el9-latest"\n""")
                                submit_file.write("""MY.SINGULARITY_EXTRA_ARGUMENTS = "-B /afs -B /cvmfs/cms.cern.ch -B /tmp -B /etc/sysconfig/ngbauth-submit -B ${XDG_RUNTIME_DIR} -B /eos --env KRB5CCNAME='FILE:${XDG_RUNTIME_DIR}/krb5cc'"\n""")
                            submit_file.write("max_retries = 3\n")
                            submit_file.write("requirements = Machine =!= LastRemoteHost\n")
                            submit_file.write(f'+JobFlavour = "microcentury"\n')
                            submit_file.write(f"queue {i}\n")

                    else:
                        if opt.condor_logs != "":
                            job_file_executable = os.path.join(CONDOR_PATH, f"{file}.sh")
                            job_file_submit = os.path.join(CONDOR_PATH, f"{file}.sub")
                        else:
                            job_file_executable = os.path.join(OUT_PATH, f"{file}.sh")
                            job_file_submit = os.path.join(OUT_PATH, f"{file}.sub")

                        if not opt.make_condor_logs:
                            job_file_out = "/dev/null"
                            job_file_err = "/dev/null"
                            job_file_log = "/dev/null"
                        elif opt.condor_logs != "":
                            job_file_out = os.path.join(CONDOR_PATH, f"{file}.$(ClusterId).$(ProcId).out")
                            job_file_err = os.path.join(CONDOR_PATH, f"{file}.$(ClusterId).$(ProcId).err")
                            job_file_log = os.path.join(CONDOR_PATH, f"{file}.$(ClusterId).log")
                        else:
                            job_file_out = os.path.join(OUT_PATH, f"{file}.$(ClusterId).$(ProcId).out")
                            job_file_err = os.path.join(OUT_PATH, f"{file}.$(ClusterId).$(ProcId).err")
                            job_file_log = os.path.join(OUT_PATH, f"{file}.$(ClusterId).log")

                        with open(job_file_executable, "w") as executable_file:
                            executable_file.write("#!/bin/sh\n")
                            if os.path.exists(f"{OUT_PATH}/merged/{file}/{file}_merged.parquet"):
                                raise Exception(
                                    f"The selected target path: {OUT_PATH}/merged/{file}/{file}_merged.parquet already exists"
                                )
                            if not os.path.exists(f'{OUT_PATH}/merged/Data_{file.split("_")[-1]}'):
                                MKDIRP(f'{OUT_PATH}/merged/Data_{file.split("_")[-1]}')
                            os.chdir(SCRIPT_DIR)
                            print(f'merge_parquet.py --source {IN_PATH}/{file}/nominal --target {OUT_PATH}/merged/Data_{file.split("_")[-1]}/{file}_ --cats {cat_dict_loc} --is-data --abs {genBinning_str}')
                            executable_file.write(f"if [ $1 -eq 0 ]; then\n")
                            executable_file.write(f"    merge_parquet.py --source {IN_PATH}/{file}/nominal --target {OUT_PATH}/merged/Data_{file.split('_')[-1]}/{file}_ --cats {cat_dict_loc} --is-data --abs {genBinning_str} || exit 107\n")
                            executable_file.write("exit 0\n")
                            executable_file.write("fi\n")
                            
                        os.system(f"chmod 775 {job_file_executable}")
                        with open(job_file_submit, "w") as submit_file:
                            if opt.condor_logs != "": submit_file.write(f"initialdir = {CONDOR_PATH}\n")
                            submit_file.write(f"executable = {job_file_executable}\n")
                            submit_file.write("arguments = $(ProcId)\n")
                            submit_file.write(f"output = {job_file_out}\n")
                            submit_file.write(f"error = {job_file_err}\n")
                            submit_file.write(f"log = {job_file_log}\n")
                            submit_file.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
                            submit_file.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
                            # if opt.max_materialize != "": submit_file.write(f"max_materialize = {opt.max_materialize}\n")
                            if opt.apptainer:
                                submit_file.write("MY.XRDCP_CREATE_DIR     = True\n")
                                submit_file.write("""MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/higgsdna:lxplus-el9-latest"\n""")
                                submit_file.write("""MY.SINGULARITY_EXTRA_ARGUMENTS = "-B /afs -B /cvmfs/cms.cern.ch -B /tmp -B /etc/sysconfig/ngbauth-submit -B ${XDG_RUNTIME_DIR} -B /eos --env KRB5CCNAME='FILE:${XDG_RUNTIME_DIR}/krb5cc'"\n""")
                            submit_file.write("max_retries = 3\n")
                            submit_file.write("requirements = Machine =!= LastRemoteHost\n")
                            submit_file.write(f'+JobFlavour = "microcentury"\n')
                            submit_file.write(f"queue\n")
                if opt.condor_logs != "":
                    submit_jobs(CONDOR_PATH)
                else:
                    submit_jobs(OUT_PATH)

            # at this point Data will be split in eras if any Data dataset is present, here we merge them again in one allData file to rule them all
            # we also skip this step if there is no Data
            if opt.merge_data:
                j = 0
                for file in files:
                    if j != 0: continue
                    file = file.split("\n")[0]  # otherwise it contains an end of line and messes up the os.walk() call
                    if opt.condor_logs != "":
                        job_file_executable = os.path.join(CONDOR_PATH, f"{file}_merge_data.sh")
                        job_file_submit = os.path.join(CONDOR_PATH, f"{file}_merge_data.sub")
                    else:
                        job_file_executable = os.path.join(OUT_PATH, f"{file}_merge_data.sh")
                        job_file_submit = os.path.join(OUT_PATH, f"{file}_merge_data.sub")

                    if not opt.make_condor_logs:
                        job_file_out = "/dev/null"
                        job_file_err = "/dev/null"
                        job_file_log = "/dev/null"
                    elif opt.condor_logs != "":
                        job_file_out = os.path.join(CONDOR_PATH, f"{file}_merge_data.$(ClusterId).$(ProcId).out")
                        job_file_err = os.path.join(CONDOR_PATH, f"{file}_merge_data.$(ClusterId).$(ProcId).err")
                        job_file_log = os.path.join(CONDOR_PATH, f"{file}_merge_data.$(ClusterId).log")
                    else:
                        job_file_out = os.path.join(OUT_PATH, f"{file}_merge_data.$(ClusterId).$(ProcId).out")
                        job_file_err = os.path.join(OUT_PATH, f"{file}_merge_data.$(ClusterId).$(ProcId).err")
                        job_file_log = os.path.join(OUT_PATH, f"{file}_merge_data.$(ClusterId).log")
                    if "data" in file.lower() or "DoubleEG" in file:
                        with open(job_file_executable, "w") as executable_file:
                            executable_file.write("#!/bin/sh\n")
                            dirpath, dirnames, filenames = next(os.walk(f'{OUT_PATH}/merged/Data_{file.split("_")[-1]}'))
                            if len(filenames) > 0:
                                print(f'merge_parquet.py --source {OUT_PATH}/merged/Data_{file.split("_")[-1]} --target {OUT_PATH}/merged/Data_{file.split("_")[-1]}/allData_ --cats {cat_dict_loc} --is-data --abs {genBinning_str}')
                                executable_file.write(f"if [ $1 -eq 0 ]; then\n")
                                executable_file.write(f"    merge_parquet.py --source {OUT_PATH}/merged/Data_{file.split('_')[-1]} --target {OUT_PATH}/merged/Data_{file.split('_')[-1]}/allData_ --cats {cat_dict_loc} --is-data --abs {genBinning_str} || exit 107\n")
                                executable_file.write("exit 0\n")
                                executable_file.write("fi\n")
                                #break
                            else:
                                logger.info(f'No merged parquet found for {file} in the directory: {OUT_PATH}/merged/Data_{file.split("_")[-1]}')
                        with open(job_file_submit, "w") as submit_file:
                            if opt.condor_logs != "": submit_file.write(f"initialdir = {CONDOR_PATH}\n")
                            submit_file.write(f"executable = {job_file_executable}\n")
                            submit_file.write("arguments = $(ProcId)\n")
                            submit_file.write(f"output = {job_file_out}\n")
                            submit_file.write(f"error = {job_file_err}\n")
                            submit_file.write(f"log = {job_file_log}\n")
                            submit_file.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
                            submit_file.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
                            # if opt.max_materialize != "": submit_file.write(f"max_materialize = {opt.max_materialize}\n")
                            if opt.apptainer:
                                submit_file.write("MY.XRDCP_CREATE_DIR     = True\n")
                                submit_file.write("""MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/higgsdna:lxplus-el9-latest"\n""")
                                submit_file.write("""MY.SINGULARITY_EXTRA_ARGUMENTS = "-B /afs -B /cvmfs/cms.cern.ch -B /tmp -B /etc/sysconfig/ngbauth-submit -B ${XDG_RUNTIME_DIR} -B /eos --env KRB5CCNAME='FILE:${XDG_RUNTIME_DIR}/krb5cc'"\n""")
                            submit_file.write("max_retries = 3\n")
                            submit_file.write("requirements = Machine =!= LastRemoteHost\n")
                            submit_file.write(f'+JobFlavour = "microcentury"\n')
                            submit_file.write(f"queue\n")
                    os.system(f"chmod 775 {job_file_executable}")
                    j += 1
                if opt.condor_logs != "":
                    submit_jobs(CONDOR_PATH, "merge_data")
                else:
                    submit_jobs(OUT_PATH, "merge_data")

    if opt.root:
        logger.info("Starting root step")
        if opt.syst:
            logger.info("you've selected the run with systematics")
            args = "--do_syst"
        else:
            logger.info("you've selected the run without systematics")
            args = ""

        if opt.merge:
            IN_PATH = OUT_PATH
        # Note, in my version of HiggsDNA I run the analysis splitting data per Era in different datasets
        # the treatment of data here is tested just with that structure
        with open(dirlist_path) as fl:
            files = fl.readlines()
            for file in files:
                file = file.split("\n")[0]
                if "data" not in file.lower() and file in process_dict:
                    if opt.condor_logs != "":
                        job_file_executable = os.path.join(CONDOR_PATH, f"{file}_root.sh")
                    else:
                        job_file_executable = os.path.join(OUT_PATH, f"{file}_root.sh")

                    if not opt.merge:
                        if opt.condor_logs != "":
                            job_file_submit = os.path.join(CONDOR_PATH, f"{file}_root.sub")
                        else:
                            job_file_submit = os.path.join(OUT_PATH, f"{file}_root.sub")
                        if not opt.make_condor_logs:
                            job_file_out = "/dev/null"
                            job_file_err = "/dev/null"
                            job_file_log = "/dev/null"
                        elif opt.condor_logs != "":
                            job_file_out = os.path.join(CONDOR_PATH, f"{file}_root.$(ClusterId).$(ProcId).out")
                            job_file_err = os.path.join(CONDOR_PATH, f"{file}_root.$(ClusterId).$(ProcId).err")
                            job_file_log = os.path.join(CONDOR_PATH, f"{file}_root.$(ClusterId).log")
                        else:
                            job_file_out = os.path.join(OUT_PATH, f"{file}_root.$(ClusterId).$(ProcId).out")
                            job_file_err = os.path.join(OUT_PATH, f"{file}_root.$(ClusterId).$(ProcId).err")
                            job_file_log = os.path.join(OUT_PATH, f"{file}_root.$(ClusterId).log")
                    with open(job_file_executable, "w") as executable_file:
                        executable_file.write("#!/bin/sh\n")
                        if os.path.exists(f"{OUT_PATH}/root/{file}"):
                            raise Exception(
                                f"The selected target path: {OUT_PATH}/root/{file} already exists"
                            )

                        if os.listdir(f"{IN_PATH}/merged/{file}/"):
                            logger.info(f"Found merged files {IN_PATH}/merged/{file}/")
                        else:
                            raise Exception(f"Merged parquet not found at {IN_PATH}/merged/")
                        MKDIRP(f"{OUT_PATH}/root/{file}")
                        os.chdir(SCRIPT_DIR)
                        executable_file.write(f"if [ $1 -eq 0 ]; then\n")
                        executable_file.write(f"    convert_parquet_to_root.py {IN_PATH}/merged/{file}/merged.parquet {OUT_PATH}/root/{file}/merged.root mc --process {process_dict[file]} {args} --cats {cat_dict_loc} --vars {var_dict_loc} --abs {genBinning_str} || exit 107\n")
                        executable_file.write("exit 0\n")
                        executable_file.write("fi\n")
                    os.system(f"chmod 775 {job_file_executable}")
                    with open(job_file_submit, "w") as submit_file:
                        if opt.condor_logs != "": submit_file.write(f"initialdir = {CONDOR_PATH}\n")
                        submit_file.write(f"executable = {job_file_executable}\n")
                        submit_file.write("arguments = $(ProcId)\n")
                        submit_file.write(f"output = {job_file_out}\n")
                        submit_file.write(f"error = {job_file_err}\n")
                        submit_file.write(f"log = {job_file_log}\n")
                        submit_file.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
                        submit_file.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
                        # if opt.max_materialize != "": submit_file.write(f"max_materialize = {opt.max_materialize}\n")
                        if opt.apptainer:
                            submit_file.write("MY.XRDCP_CREATE_DIR     = True\n")
                            submit_file.write("""MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/higgsdna:lxplus-el9-latest"\n""")
                            submit_file.write("""MY.SINGULARITY_EXTRA_ARGUMENTS = "-B /afs -B /cvmfs/cms.cern.ch -B /tmp -B /etc/sysconfig/ngbauth-submit -B ${XDG_RUNTIME_DIR} -B /eos --env KRB5CCNAME='FILE:${XDG_RUNTIME_DIR}/krb5cc'"\n""")
                        submit_file.write("max_retries = 3\n")
                        submit_file.write("requirements = Machine =!= LastRemoteHost\n")
                        submit_file.write(f'+JobFlavour = "microcentury"\n')
                        submit_file.write(f"queue\n")
                elif "data" in file.lower():
                    if opt.condor_logs != "":
                        job_file_executable = os.path.join(CONDOR_PATH, f"{file}_root.sh")
                    else:
                        job_file_executable = os.path.join(OUT_PATH, f"{file}_root.sh")

                    if not opt.merge:
                        if opt.condor_logs != "":
                            job_file_submit = os.path.join(CONDOR_PATH, f"{file}_root.sub")
                        else:
                            job_file_submit = os.path.join(OUT_PATH, f"{file}_root.sub")

                        if not opt.make_condor_logs:
                            job_file_out = "/dev/null"
                            job_file_err = "/dev/null"
                            job_file_log = "/dev/null"
                        elif opt.condor_logs != "":
                            job_file_out = os.path.join(CONDOR_PATH, f"{file}_root.$(ClusterId).$(ProcId).out")
                            job_file_err = os.path.join(CONDOR_PATH, f"{file}_root.$(ClusterId).$(ProcId).err")
                            job_file_log = os.path.join(CONDOR_PATH, f"{file}_root.$(ClusterId).log")
                        else:
                            job_file_out = os.path.join(OUT_PATH, f"{file}_root.$(ClusterId).$(ProcId).out")
                            job_file_err = os.path.join(OUT_PATH, f"{file}_root.$(ClusterId).$(ProcId).err")
                            job_file_log = os.path.join(OUT_PATH, f"{file}_root.$(ClusterId).log")
                    
                    with open(job_file_executable, "w") as executable_file:
                        executable_file.write("#!/bin/sh\n")
                        if os.listdir(f'{IN_PATH}/merged/Data_{file.split("_")[-1]}/'):
                            logger.info(
                                f'Found merged data files in: {IN_PATH}/merged/Data_{file.split("_")[-1]}/'
                            )
                        else:
                            raise Exception(
                                f'Merged parquet not found at: {IN_PATH}/merged/Data_{file.split("_")[-1]}/'
                            )

                        if os.path.exists(
                            f'{OUT_PATH}/root/Data/allData_{file.split("_")[-1]}.root'
                        ):
                            logger.info(
                                f'Data already converted: {OUT_PATH}/root/Data/allData_{file.split("_")[-1]}.root'
                            )
                            continue
                        elif not os.path.exists(f"{OUT_PATH}/root/Data/"):
                            MKDIRP(f"{OUT_PATH}/root/Data")
                            os.chdir(SCRIPT_DIR)
                            executable_file.write(f"if [ $1 -eq 0 ]; then\n")
                            executable_file.write(f"    convert_parquet_to_root.py {IN_PATH}/merged/Data_{file.split('_')[-1]}/allData_merged.parquet {OUT_PATH}/root/Data/allData_{file.split('_')[-1]}.root data --cats {cat_dict_loc} --vars {var_dict_loc} --abs {genBinning_str} || exit 107\n")
                            executable_file.write("exit 0\n")
                            executable_file.write("fi\n")
                        else:
                            os.chdir(SCRIPT_DIR)
                            executable_file.write(f"if [ $1 -eq 0 ]; then\n")
                            executable_file.write(f"    convert_parquet_to_root.py {IN_PATH}/merged/Data_{file.split('_')[-1]}/allData_merged.parquet {OUT_PATH}/root/Data/allData_{file.split('_')[-1]}.root data --cats {cat_dict_loc} --vars {var_dict_loc} --abs {genBinning_str} || exit 107\n")
                            executable_file.write("exit 0\n")
                            executable_file.write("fi\n")
                os.system(f"chmod 775 {job_file_executable}")
                with open(job_file_submit, "w") as submit_file:
                    if opt.condor_logs != "": submit_file.write(f"initialdir = {CONDOR_PATH}\n")
                    submit_file.write(f"executable = {job_file_executable}\n")
                    submit_file.write("arguments = $(ProcId)\n")
                    submit_file.write(f"output = {job_file_out}\n")
                    submit_file.write(f"error = {job_file_err}\n")
                    submit_file.write(f"log = {job_file_log}\n")
                    submit_file.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
                    submit_file.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
                    # if opt.max_materialize != "": submit_file.write(f"max_materialize = {opt.max_materialize}\n")
                    if opt.apptainer:
                        submit_file.write("MY.XRDCP_CREATE_DIR     = True\n")
                        submit_file.write("""MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/higgsdna:lxplus-el9-latest"\n""")
                        submit_file.write("""MY.SINGULARITY_EXTRA_ARGUMENTS = "-B /afs -B /cvmfs/cms.cern.ch -B /tmp -B /etc/sysconfig/ngbauth-submit -B ${XDG_RUNTIME_DIR} -B /eos --env KRB5CCNAME='FILE:${XDG_RUNTIME_DIR}/krb5cc'"\n""")
                    submit_file.write("max_retries = 3\n")
                    submit_file.write("requirements = Machine =!= LastRemoteHost\n")
                    submit_file.write(f'+JobFlavour = "microcentury"\n')
                    submit_file.write(f"queue\n")
        if opt.condor_logs != "":
            submit_jobs(CONDOR_PATH, "root")
        else:
            submit_jobs(OUT_PATH, "root")

    if not opt.apptainer:
    # We don't want to leave trash around
        if os.path.exists(dirlist_path):
            os.system(f"rm {dirlist_path}")
        if opt.output == "":
            if os.path.exists(f"{SCRIPT_DIR}/../../higgs_dna/category.json"):
                os.system(f"rm {SCRIPT_DIR}/../../higgs_dna/category.json")
            if os.path.exists(f"{SCRIPT_DIR}/../../higgs_dna/variation.json"):
                os.system(f"rm {SCRIPT_DIR}/../../higgs_dna/variation.json")

    # elif (opt.apptainer) and (opt.root):
    #     if os.path.exists(f"{OUT_PATH}/category.json"):
    #         os.system(f"rm {OUT_PATH}/category.json")
    #     if os.path.exists(f"{OUT_PATH}/variation.json"):
    #         os.system(f"rm {OUT_PATH}/variation.json")