#!/usr/bin/env python
import argparse
import json
import ast
import os
import glob
import awkward as ak
from higgs_dna.utils.logger_utils import setup_logger
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import numpy as np
from importlib import resources
from higgs_dna.scripts.postprocessing.tools.Btag_WeightSum_Calculation import Get_WeightSum_Btag, Renormalize_BTag_Weights
from higgs_dna.scripts.postprocessing.tools.postprocessing_tools import filter_and_set_diff_variable
from coffea.processor.accumulator import iadd

def process_custom_accumulator(source_path, logger):
    # Get the custom accumulator from all files in the source path
    accumulator = None
    source_files = glob.glob("%s/*.parquet" % source_path)
    for f in source_files:
        try:
            file_accumulator = pq.read_schema(f).metadata[b'custom_accumulator']
        except KeyError:
            logger.warning(f"Custom accumulator requested but not found in file {f}")
        file_accumulator = json.loads(file_accumulator)
        if accumulator is None:
            accumulator = file_accumulator
        else:
            accumulator = iadd(accumulator, file_accumulator)
    return accumulator


def main():
    parser = argparse.ArgumentParser(
        description="Simple utility script to merge all parquet files in one folder."
    )
    parser.add_argument(
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Debugging verbosity for logger.",
    )
    parser.add_argument(
        "--source",
        type=str,
        default="",
        help="Comma separated paths (with trailing slash) to folder where multiple parquet files are located. Careful: Folder should ONLY contain parquet files!",
    )
    parser.add_argument(
        "--target",
        type=str,
        default="",
        help="Comma separated paths (with trailing slash) to desired folder. Resulting merged file is placed there.",
    )
    parser.add_argument(
        "--cats",
        type=str,
        dest="cats_dict",
        default="",
        help="Dictionary containing category selections.",
    )
    parser.add_argument(
        "--is-data",
        default=False,
        action="store_true",
        help="Files to be merged are data and therefore do not require normalisation.",
    )
    parser.add_argument(
        "--merge-all-data",
        default=False,
        action="store_true",
        help="To be used for merging 'Data*' parquets to a single 'allData' parquet",
    )
    parser.add_argument(
        "--skip-normalisation",
        default=False,
        action="store_true",
        help="Independent of file type, skip normalisation step",
    )
    parser.add_argument(
        "--abs",
        dest="abs",
        action="store_true",
        default=False,
        help="Uses absolute path for the dictionary files.",
    )
    parser.add_argument(
        "--genBinning",
        type=str,
        dest="genBinning",
        default="",
        help="Optional: Path to the JSON containing the binning at gen-level.",
    )
    parser.add_argument(
    "--do-b-weight-normalisation",
    default=False,
    action="store_true",
    help="Perform the bweight normalization to make sure the number of event remain the same before and after apling the b tagging weights",
   )
    parser.add_argument(
        "--custom-accumulator",
        default=False,
        action="store_true",
        dest="custom_accumulator",
        help="If set, the script will process the custom accumulator from the parquet files.",
    )

    args = parser.parse_args()
    source_paths = args.source.split(",")
    target_paths = args.target.split(",")

    BASEDIR = resources.files("higgs_dna").joinpath("")

    if args.genBinning != "":
        if args.abs:
            genBinning_path = os.path.realpath(args.genBinning)
        else:
            genBinning_path = os.path.join(BASEDIR, "scripts/postprocessing/sample_gen_binning.json")
        with open(genBinning_path, 'r') as json_file:
            gen_binning = json.load(json_file)
    else:
        gen_binning = None

    logger_verbosity = "DEBUG" if args.verbose else "INFO"

    logger = setup_logger(level=logger_verbosity)

    if (
        (len(source_paths) != len(target_paths))
        or (args.source == "")
        or (args.target == "")
    ):
        logger.info("You gave a different number of sources and targets")
        exit


    if args.cats_dict != "":
        if args.abs:
            cats_path = os.path.realpath(args.cats_dict)
        else:
            cats_path = os.path.join(BASEDIR, "category.json")
        with open(cats_path) as pf:
            cat_dict = json.load(pf)
        for cat in cat_dict:
            logger.info(f"Found category: {cat}")
    else:
        logger.info(
            "You provided an invalid dictionary containing categories information, have a look at your version of prepare_output_file.py"
        )
        logger.info(
            "An inclusive NOTAG category is used as default"
        )
        cat_dict = {"NOTAG": {"cat_filter": [("pt", ">", -1.0)]}}

    if (not args.is_data) & (not args.skip_normalisation):
        logger.info(
            "Extracting sum of gen weights (before selection) from metadata of files to be merged."
        )
        if(args.do_b_weight_normalisation):
            IsBtagNorm_sys_arr,WeightSum_preBTag_arr,WeightSum_postBTag_arr,WeightSum_postBTag_sys_arr = Get_WeightSum_Btag(source_paths,logger)

        sum_genw_beforesel_arr = []
        for i, source_path in enumerate(source_paths):
            source_files = glob.glob("%s/*.parquet" % source_path)
            sum_genw_beforesel = 0
            for f in source_files:
                sum_genw_beforesel += float(pq.read_schema(f).metadata[b'sum_genw_presel'])
            sum_genw_beforesel_arr.append(sum_genw_beforesel)
        logger.info(
            "Successfully extracted sum of gen weights (before selection)"
        )

    for i, source_path in enumerate(source_paths):
        # Process custom accumulator
        if args.custom_accumulator:
            logger.info(f"Processing custom accumulator for {source_path}")
            custom_accumulator = process_custom_accumulator(source_path, logger)

        for cat in cat_dict:
            logger.info("-" * 125)
            logger.info(
                f"INFO: Starting parquet file merging. Attempting to read parquet dataset from {source_path}, for category: {cat}"
            )
            dataset = ds.dataset(glob.glob(source_path+"/Data*.parquet") if args.merge_all_data else source_path)
            logger.info("Parquet dataset read successfully.")
            logger.info(
                f"Attempting to merge parquet dataset and save to {target_paths[i]}."
            )
            if "Data" in target_paths[i]:
                os.makedirs("/".join(target_paths[i].split("/")[:-1]), exist_ok=True)
            else:
                os.makedirs(target_paths[i], exist_ok=True) # Create target directory if it does not exist

            # Process in batches
            output_file = target_paths[i] + cat + "_merged.parquet"

            # Process in batches
            writer = None

            for batch in dataset.to_batches(filter=pq.filters_to_expression(cat_dict[cat]["cat_filter"])):
                batch_arr = ak.from_arrow(batch)

                if (not args.is_data) & (not args.skip_normalisation):
                    if gen_binning != None:
                        for keys in gen_binning:
                            var_dict = {ast.literal_eval(key): value for key, value in gen_binning[keys].items()}
                            if len(list(var_dict.keys())[0]) == 5:
                                selectionVariableName = list(var_dict.keys())[0][0]
                                var_dict = {k[1:]: v for k, v in var_dict.items()}
                            else:
                                selectionVariableName = keys
                            batch_arr = filter_and_set_diff_variable(batch_arr, var_dict, selectionVariableName, "diffVariable_" + keys)

                    batch_arr['weight_nominal'] = batch_arr['weight']
                    syst_weight_fields = [field for field in batch_arr.fields if (("weight_" in field) and ("Up" in field or "Down" in field))]
                    for weight_field in ["weight"] + syst_weight_fields:
                        batch_arr[weight_field] = batch_arr[weight_field] / sum_genw_beforesel_arr[i]
                    logger.info("Successfully added normalised weight column")

                    if args.do_b_weight_normalisation:
                        if((WeightSum_preBTag_arr[i]/WeightSum_postBTag_arr[i])!=1):
                            batch_arr = Renormalize_BTag_Weights(batch_arr, target_paths[i], cat, WeightSum_preBTag_arr[i], WeightSum_postBTag_arr[i], WeightSum_postBTag_sys_arr[i], IsBtagNorm_sys_arr[i], logger)

                table = ak.to_arrow_table(batch_arr, extensionarray=False)
                if args.custom_accumulator:
                    logger.info("Adding custom accumulator")
                    table = table.replace_schema_metadata({b'custom_accumulator': json.dumps(custom_accumulator).encode("utf-8")})
                    logger.info("Custom accumulator added successfully")
                else:
                    table = table.replace_schema_metadata()

                if writer is None:
                    schema = table.schema
                    writer = pq.ParquetWriter(output_file, schema)
                writer.write_table(table)

            if writer:
                writer.close()

            logger.info(f"Success! Merged parquet file is located in {output_file}")
            logger.info("-" * 125)

if __name__ == "__main__":
    main()
