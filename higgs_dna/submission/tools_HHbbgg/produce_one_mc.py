import os
import argparse
import json


def write_sample_txt(keyword, year, cmsdas):
    filename = f"samples_mc_{year}_{keyword}.txt"
    with open(filename, "w") as f:
        f.write(keyword + " " + cmsdas)
    return filename


def fetch_datasets(sample_file, dbs_instance='prod/global', region='Yolo'):
    command = f"python scripts/samples/fetch_datasets.py -i {sample_file} -w {region} --dbs-instance {dbs_instance}"
    os.system(command)


def run_analysis(nano_version, parent_dir, keyword, year, memory):
    memoryLine = f"--memory {memory} " if memory is not None else ""
    command = (
        f"python scripts/run_analysis.py "
        f"--json-analysis runner_mc_{year}_{keyword}.json "
        f"--dump {parent_dir} "
        f"--doFlow-corrections "
        f"--fiducialCuts store_flag "
        f"--Smear-sigma-m "
        f"--doDeco "
        f"--executor vanilla_lxplus "
        # f"--queue longlunch "
        f"--queue workday "
        f"{memoryLine}"
        f"--debug "
        f"--nano-version {nano_version} "
        f"--timeout 200"
    )
    print(command)
    os.system(command)


def update_json_config(keyword, year):
    # Using the preliminary JSON - to be updated for final results
    with open("submission/tools_HHbbgg/prelim_runner_mc_template.json", "r") as f:
        config = json.load(f)
    config["samplejson"] = f"samples_mc_{year}_{keyword}.json"
    if "year" in config:
        config["year"].pop("GluGluToHH", None)
        config["year"][keyword] = [f"{year}"]
    if "corrections" in config:
        config["corrections"][keyword] = config["corrections"].pop("GluGluToHH", [])
        if "GluGluHtoGG" in keyword:
            config["corrections"][keyword].append("NNLOPS")
        if "2024" in year:
            if any("jet" in corr for corr in config["corrections"][keyword]):
                jerc_idxs = [
                    idx for idx, corr in enumerate(config["corrections"][keyword])
                    if "jet" in corr
                ]
                for jerc_idx in jerc_idxs:
                    # Remove '_syst' from jerc as we don't have for 2024/2025 yet
                    config["corrections"][keyword][jerc_idx] = config["corrections"][keyword][jerc_idx].replace("_syst", "")
            # Remove bTag SF correction for now as we don't have SFs
            if any("PNet_bTagShapeSF" == corr for corr in config["corrections"][keyword]):
                config["corrections"][keyword].remove("PNet_bTagShapeSF")
    if "systematics" in config:
        config["systematics"][keyword] = config["systematics"].pop("GluGluToHH", [])
        if "2024" in year:
            if any("PNet_bTagShapeSF" == syst for syst in config["systematics"][keyword]):
                config["systematics"][keyword].remove("PNet_bTagShapeSF")
    new_filename = f"runner_mc_{year}_{keyword}.json"
    with open(new_filename, "w") as f:
        json.dump(config, f, indent=4)

    return new_filename


def main():
    parser = argparse.ArgumentParser(description="Run MC production example pipeline.")
    parser.add_argument("-k", "--keyword", required=True, help="Keyword for dataset filtering")
    parser.add_argument("-c", "--cmsdas", required=True, help="Keyword for cmsdas filtering")
    parser.add_argument("-p", "--parent-dir", required=True, help="Directory to store output parquets")
    parser.add_argument("-y", "--year", required=True, choices=["2022postEE","2022preEE","2023postBPix","2023preBPix", "2024"], help="year")
    parser.add_argument("-n", "--nano", required=True, help="nano-version")
    parser.add_argument("-m", "--memory", help="condor job memory")
    parser.add_argument(
        "-w",
        "--where",
        help="Specify the region for xrootd prefix (only for grid mode).",
        default="Eurasia",
        choices=["Americas", "Eurasia", "Yolo"],
    )
    parser.add_argument(
        "--dbs-instance",
        dest="instance",
        help="The DBS instance to use for querying datasets (only for grid mode).",
        type=str,
        default="prod/global",
        choices=["prod/global", "prod/phys01", "prod/phys02", "prod/phys03"],
    )

    args = parser.parse_args()

    # Write dataset sample file
    sample_file = write_sample_txt(args.keyword, args.year, args.cmsdas)

    # Fetch datasets
    fetch_datasets(sample_file, dbs_instance=args.instance, region=args.where)

    # Update and save JSON configuration
    update_json_config(args.keyword, args.year)

    # Launch jobs
    run_analysis(args.nano, args.parent_dir, args.keyword, args.year, args.memory)


if __name__ == "__main__":
    main()
