#!/usr/bin/env python
from higgs_dna.utils.logger_utils import setup_logger
from higgs_dna.utils.runner_utils import get_proxy
from XRootD import client
from rich.progress import track
import argparse
import json
import os


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# - This script creates a json file containing the unprocessed samples, based on a provided sample.json and parquet directory. ----------------------------------------------------#
# - It can work with the DAS UUID naming convention or the Legacy one, where the UUID is contained in the ROOT file header. -------------------------------------------------------#
# - It will check for missing UUIDs in the parquet file names & unprocessed event in individual parquets, and will produce a json with the updated list of unprocessed samples. ---#
# - EXAMPLE USAGE: ----------------------------------------------------------------------------------------------------------------------------------------------------------------#
# - python3 get_unprocessed_files.py --convention <naming_convention> --source <dir_to_HiggsDNA_dump> --json <sample.json> --output <some_path/unprocessed_samples.json>
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

def get_fetcher_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Obtain the list of unprocessed root files from the associated samples list."
    )
    parser.add_argument(
        "--log",
        dest="log",
        type=str,
        default="INFO",
        help="Logger info level"
    )
    parser.add_argument(
        "-s",
        "--source",
        help="Directory containing the datasets.",
        required=True,
    )
    parser.add_argument(
        "-j",
        "--json",
        help="Json containing the lists of samples.",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--convention",
        help="Parquet files naming convention: DAS or Legacy.",
        required=True,
        choices=["DAS", "Legacy"]
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path in which to save the unprocessed samples, eg 'Thisdir/myoutput.json'.",
        default="unprocessed_samples.json",
    )

    return parser.parse_args()


# This function takes a list of event chunks associated with one parquet, eg [0-127,127-255,255-511]
# It will return the total number of events if the chunks are continuous (eg 511) and 0 otherwise.
def check_range(range_list):

    # Create list of event chunks and sort it
    chunks = [int(val) for rg in range_list for val in rg.split('-')]
    chunks.sort()

    total_range = chunks[-1] - chunks[0]
    # Verify that chunks in the list are successive.
    for i in range(1,len(chunks)-1,2):
        if chunks[i] != chunks[i+1]:
            total_range = 0
            break

    return total_range


# Use xrootd to read the header of a root file.
# Retrieve the associated UUID from the file header.
def get_root_uuid_using_xrootd(fpath):

    logger.debug(f"Opening {fpath} with xrootd")

    with client.File() as f:
        # Try to open the file
        status, _ = f.open(fpath)
        if not status.ok:
            raise RuntimeError(f"Failed to open file {fpath}: {status.message}")

        # Read version bytes (offset 4, size 4)
        status, version_bytes = f.read(offset=4, size=4)
        if not status.ok:
            raise RuntimeError(f"Failed to read version from file {fpath}: {status.message}")
        # Convert version to an integer
        version = int.from_bytes(version_bytes, "big")

        # Determine the offset for UUID based on the version
        uuid_offset = 59 if version >= 1000000 else 47
        # Read the UUID bytes
        status, uuid_bytes = f.read(offset=uuid_offset, size=16)
        if not status.ok:
            raise RuntimeError(f"Failed to read UUID from file {fpath}: {status.message}")
        # Convert UUID bytes to standard UUID format
        root_uuid = '-'.join([
            uuid_bytes[0:4].hex(),
            uuid_bytes[4:6].hex(),
            uuid_bytes[6:8].hex(),
            uuid_bytes[8:10].hex(),
            uuid_bytes[10:16].hex()
        ])

        return root_uuid


# Create a dict of form {'dataset':{'uuid':nevent}} from the source directory.
def create_pq_dict(path, root_dict):
    tree =  {}
    parquetinfo = {}
    source_dict = {}

    dataset_from_json = list(root_dict.keys())

    # List all subdirectories in provided path. Subdirs name are expected to be the same as datasets in sample.json
    subdir = [
            d for d in os.listdir(path) if (os.path.isdir(os.path.join(path, d)) and d in dataset_from_json)
        ]

    logger.info(f"Starting inspection of source directory {path}.")

    if not subdir:
        logger.warning(f"No datasets from the original JSON were found in the source directory.")

    # Create the dict containing the parquet files: 'dataset':(parquet1,parquet2,...)
    for dataset in subdir:
        logger.info(f"Starting inspection of directory {dataset}")
        # We look into the "nominal" subdirectory
        pq_path = os.path.join(path, dataset+"/nominal/")
        if not os.path.exists(pq_path):
            logger.info(f"{pq_path} does not exist, continue.")
            continue

        file_list = [
                f for f in os.listdir(pq_path) if os.path.isfile(os.path.join(pq_path, f))
            ]
        tree[dataset] = file_list
        logger.debug(f"Successfully read parquet files for dataset {dataset}")

    # Retrieve the uuid and event chunks from the parquet name.
    for dataset, file_list in tree.items():
        for filename in file_list:
            index = 1 if filename.startswith("_") else 0
            uuid = filename.split('_')[index]
            ev_range = (filename.split('_')[index+2]).replace(".parquet", "").replace(".txt", "")

            if dataset not in parquetinfo:
                parquetinfo[dataset] = {}

            if uuid not in parquetinfo[dataset]:
                parquetinfo[dataset][uuid] = [ev_range]
            else:
                parquetinfo[dataset][uuid].extend([ev_range])

    # Finally, we create the dict 'dataset':{'uuid':nevent}
    for dataset, uuidinfo in parquetinfo.items():
        source_dict[dataset] = {}
        for uuid, ranges in uuidinfo.items():
            total_range = check_range(ranges)
            source_dict[dataset][uuid] = total_range

        logger.debug(f"Successfully retrieved parquet information for dataset {dataset}")

    return(source_dict)


# Create a dict of form {'dataset':{'uuid':(nevent,physical_location)}} from the sample.json file.
def parse_sample_json(samples_json: str, convention: str):

    root_dict = {}
    rootf_uuid = []

    f = open(samples_json)
    samples = json.load(f)

    logger.info(f"Retrieving information on root files from {samples_json}.")
    for name in samples:
        logger.info(f"Retrieving file information for dataset {name}")

        # Get list of root files in dataset
        rootf_location = samples[name]
        # Remove redirector, eg ""root://xrootd-cms.infn.it/"
        rootf_name = ["/store"+ file.split("store")[-1] for file in rootf_location]

        # Here we retrieve the original dataset name
        dataset =  os.popen(
                # use the cvmfs source for dasgoclient because it works for everyone
                # Both local infrastructures with cvmfs and lxplus!
                ("/cvmfs/cms.cern.ch/common/dasgoclient -query='dataset file={}'").format(
                    rootf_name[0]
                )
            ).read()

        # From the dataset name, we can now retrieve all the root files contained in it
        file_list = (os.popen(
                ("/cvmfs/cms.cern.ch/common/dasgoclient -query='file dataset={} | grep file.name | grep file.nevents'").format(
                    dataset.strip()
                )
            ).read()).splitlines()

        if convention == "DAS":
            # Retrieve DAS uuid from root file name
            rootf_uuid = [file.split("/")[-1].replace(".root","") for file in rootf_name]
        elif convention == "Legacy":
            # Retrieve ROOT uuid from root file location using xrootd
            rootf_uuid = [get_root_uuid_using_xrootd(file) for file in track(rootf_location, description=f"[blue]Processing files in {name}...")]

        # Construct the output dict for each dataset
        root_dict[name] = {}
        for f in file_list:
            fname, nev = f.split()
            nev = int(nev)
            # We only want the information for the root files specified in the json
            if fname in rootf_name:
                # Find the associated information from the list: uuid and physical file location
                index = rootf_name.index(fname)
                associated_rootuuid = rootf_uuid[index]
                associated_rootf_location = rootf_location[index]
                root_dict[name][associated_rootuuid] = (nev, associated_rootf_location)

        logger.debug(f"Successfully retrieved file information for dataset {name}")

    return root_dict


# Compare two dataset (uuid and nevents) and return the list of unprocessed samples.
def compare_data(root_dict, source_dict):

    unprocessed_files = []

    root_uuid = list(root_dict.keys())
    parquet_uuid = list(source_dict.keys())

    for uuid in root_uuid:
        # Retrieve number of events and physical file location
        rt_nevent, file_location = root_dict.get(uuid, 0)

        # Check if uuid is found in parquet list
        if uuid not in parquet_uuid:
            logger.debug(f"Missing file with uuid: {uuid}")
            unprocessed_files.append(file_location)
            continue

        # Check if all events have been processed
        pq_nevent = source_dict.get(uuid, 0)

        if not rt_nevent == pq_nevent:
            logger.debug(f"Missing event for uuid: {uuid}. Expected {rt_nevent} and got {pq_nevent}")
            unprocessed_files.append(file_location)
            continue

    return unprocessed_files


if __name__ == "__main__":

    output_dict = {}

    args = get_fetcher_args()

    logger = setup_logger(level=args.log)

    if ".json" not in args.output:
        raise Exception("Output file must have '.json' extension and be a json file!")
    if ".json" not in args.json:
        raise Exception("Input json must have '.json' extension and be a json file!")
    if not os.path.isdir(args.source):
        raise Exception("Source directory does not exist. Make sure you provided a valid path.")

    # Check for valid proxy (required to use dasgoclient)
    get_proxy()

    # Create dicts from sample.json and parquet directory
    root_dict = parse_sample_json(args.json, args.convention)
    pq_dict = create_pq_dict(args.source, root_dict)

    logger.info("Starting creation of output file.")
    for dataset_name in root_dict:
        rt = root_dict[dataset_name]
        # Check if dataset in sample.json exists in parquet dir.
        if dataset_name not in pq_dict:
            logger.info(f"Dataset {dataset_name} could not be found in source directory and will be marked as unprocessed.")
            pq = {}
        else :
            pq = pq_dict[dataset_name]

        # Construct the dict of unprocessed root files
        unprocessed_files = compare_data(rt,pq)

        if unprocessed_files:
            output_dict[dataset_name] = unprocessed_files

        expected = len(rt)
        fully_processed = len(rt)-len(unprocessed_files)
        logger.info(f"Out of the {expected} files specified for dataset {dataset_name} in json file, {fully_processed} were fully processed. ")

    # If dict is empty (all samples fully processed), do not create output json
    if not output_dict:
        logger.info("All samples were fully processed. Output file will not be created.")

    else:
        logger.info(f"Output file will be saved in {args.output}.")
        with open(args.output, 'w') as f:
            json.dump(output_dict, f, indent=4)
