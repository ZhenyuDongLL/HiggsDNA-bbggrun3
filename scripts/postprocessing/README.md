# Run `faster_prepare_output_file.py` 

Depending on the nature of the dataset, different steps have to be followed.

## MC files

Merging the parquet files from HiggsDNA. Execute:

```
python faster_prepare_output_file.py --input absolute/path/input_for_merging --cats --catDict absolute/path/to/catDict.json --varDict absolute/path/to/varDict.json --syst --merge --output absolute/path/to/output_of_merging
```

Then, respective HTCondor jobs which merge the parquet files located in `absolute/path/input_for_merging` are created. Once they are merged, run:

```
python faster_prepare_output_file.py --input absolute/path/to/output_of_merging --cats --catDict absolute/path/to/catDict.json --varDict absolute/path/to/varDict.json --syst --root --output absolute/path/to/output_of_merging
```

This will create ROOT files in `absolute/path/to/output_of_merging` which contain all variations found in the `absolute/path/input_for_merging/considered_process` subfolder.

---

## Data files

For merging the parquet files of the data, one has to execute:

```
python faster_prepare_output_file.py --input absolute/path/input_for_merging_data --cats --catDict absolute/path/to/catDict_data.json --varDict absolute/path/to/varDict_data.json --syst --merge --output absolute/path/to/output_of_merging_data
```

This will create for each era (A-F) a parquet file. In order to merge the eras into an `allData.root` file, one has to run:

```
python faster_prepare_output_file.py --input absolute/path/input_for_merging_data --cats --catDict absolute/path/to/catDict_data.json --varDict absolute/path/to/varDict_data.json --syst --merge --output absolute/path/to/output_of_merging_data --merge-data-only
```

Then to get the ROOT file, run:

```
python faster_prepare_output_file.py --input absolute/path/to/output_of_merging_data --cats --catDict absolute/path/to/catDict_data.json --varDict absolute/path/to/varDict_data.json --syst --root --output absolute/path/to/output_of_merging_data
```
