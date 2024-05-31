=============
Main Concepts
=============


.. _def-cltool:

-----------------
Command Line Tool
-----------------
If you want to run an analysis with processors and taggers that have already been developed, the suggested way is to use the command line tool ``run_analysis.py``.

The main part that define an analysis are the following:

* ``datasets``:
  path to a JSON file in the form ``{"dataset_name": [list_of_files]}`` (like the one dumped by dasgoclient)
* ``workflow``:
  the coffea processor you want to use to process you data, can be found in the modules located inside the subpackage ``higgs_dna.workflows``
* ``metaconditions``:
  the name (without ``.json`` extension) of one of the JSON files inherited from FLASHgg and located inside ``higgs_dna.metaconditions``
* ``year``: the year condition for each sample, let us use ``2016preVFP``, ``2016postVFP``, ``2017``, ``2018``, ``2022preEE``, ``2022postEE``, ``2023``
* ``taggers``:
  the set of taggers you want to use, can be found in the modules located inside the subpackage ``higgs_dna/workflows/taggers``
* ``systematics``: the set of systematics you want to use for each sample
* ``corections``: the set of corrections you want to use for each sample

These parameters are specified in a JSON file and passed to the command line with the flag ``--json-analysis``

.. code-block:: bash

        run_analysis.py --json-analysis simple_analysis.json

where ``simple_analysis.json`` looks like this:

.. code-block:: json

      {
          "samplejson": "/work/gallim/devel/HiggsDNA/tmp/DY-data-test.json",
          "workflow": "tagandprobe",
          "metaconditions": "Era2017_legacy_xgb_v1",
          "taggers": [
              "DummyTagger1"
          ],
          "year":[
            "SampleName1": ["2022preEE"],
            "SampleName2": ["2017"]            
          ]
          "systematics": {
              "SampleName1": [
                  "SystematicA", "SystematicB"
              ],
              "SampleName2": [
                  "SystematicC"
              ]
          },
          "corrections": {
                  "SampleName1": [
                      "CorrectionA", "CorrectionB"
                  ],
                  "SampleName2": [
                      "CorrectionC"
                  ]
              }
      }

where the ``taggers`` list and the ``systematics`` or ``corrections`` dictionaries can be left empty if no taggers or systematics are applied. For most of the ``systematics``, a ``year`` condition is necessary.


The next two flags that you will want to specify are ``dump`` and ``executor``: the former receives the path to a directory where the parquet output files will be stored, while the latter specifies the Coffea executor used to process the chunks of data. It can be one of the following:

* ``iterative``
* ``futures``
* ``dask/condor``
* ``dask/slurm``
* ``dask/lpc``
* ``dask/lxplus``
* ``dask/casa``
* ``parsl/slurm``
* ``parsl/condor``

There are then a few other options that depend on the backend. When running with Dask, for instance, you may want to use the command line to change the following parameters:

* ``workers``:
  number processes/threads used **on the same node** (this means that every job will use this amount of cores)
* ``scaleout``:
  minimum number of nodes to scale out to (i.e. minimum number of jobs submitted)
* ``max-scaleout``:
  maximum number of nodes to adapt your cluster to (i.e. maximum number of jobs submitted)

As usual, a description of all the options is printed when running::

        run_analysis.py --help

--------------
Postprocessing
--------------

In order to get ROOT files from the output of `run_analysis.py`, we have to post-process the data. For that, `prepare_output_file.py` comes in handy. This script has two functionalities. First, merging the parquet files with respect to the desired categories and second, creating a ROOT file out of that. This ROOT file can be of further use (for example as input in FinalFit). One can choose between two modes: Local or HTCondor. The latter comes is very useful when you have a lot of categories.

By default, `prepare_output_file.py` uses the local execution to process files. If one wants to process the files via HTCondor, the `--condor` flag is to be used.

Data processing with local
--------------------------

To process the data locally, we have to know some things. First, we need to specify the absolute input path (`--input`) which leads to your output of `run_analysis.py` (unmerged parquet files). The output folder in which the merged parquet files are stored needs to be specified with `--output`. If one wants to categorize the files, the `--cats` keyword is used in conjunction with `--catsDict` which points to the `category.json` to be considered. Are systematics desired, they have to be activated with `--syst`.

In order to merge the parquet files according to the categories and produce the ROOT files in the same step, the following command is to be used:

.. code-block:: python

    python prepare_output_file.py --input /absolute/input/path --cats --catDict /absolute/path/to/cat_data.json --varDict /absolute/path/to/varDict_data.json --syst --merge --root --output /absolute/output/path

Using the condor-way, one has to pay attention when processing data as an additional step wrt. the local-way is required, and the merge and ROOT-production step have to be separated:

Data processing with condor
---------------------------

The first step is to merge the data parquet files according to the chosen categories. Since the data come in so-called eras (era `A`, era `B`, etc.), they have to be merged, such that we have per era and category a parquet file. This is the purpose of the following command, which has to be executed first:

.. code-block:: python

    python prepare_output_file.py --input /absolute/input/path --cats --catDict /absolute/path/to/cat_data.json --varDict /absolute/path/to/varDict_data.json --syst --merge --output /absolute/output/path --condor

Studies in the past showed that for 2022 data there is not much of a difference significance-wise between splitting `preEE` and `postEE` datasets (referencing to the ECAL Endcap water leak in 2022) and merging them. For this reason, it was merged to one big dataset for HIG-23-014. The following command merges the era datasets to an `allData.parquet` file according to the categories. One needs in addition the flag `--merge-data-only`:

.. code-block:: python

    python prepare_output_file.py --input /absolute/input/path --cats --catDict /absolute/path/to/cat_data.json --varDict /absolute/path/to/varDict_data.json --syst --merge --output /absolute/output/path --merge-data-only --condor

Finally, we convert the parquet files to ROOT:

.. code-block:: python

    python prepare_output_file.py --input /absolute/input_path/to_folder_with_merged --cats --catDict /absolute/path/to/cat_data.json --varDict /absolute/path/to/varDict_data.json --syst --root --output /absolute/input_path/to_folder_with_merged --condor

Whenever the parquet files are merged (after the first step), a folder `merged` in the `/absolute/output/path` is created. For getting the ROOT files, one has to use the folder `/absolute/output/path` (which is now containing the `merged` subfolders) as the new input folder. The file processing for MC samples functions in a similar way:

MC processing with condor
-------------------------

Similar to data, the MC samples can be processed with HTCondor. Here we only have two steps. The first consists of merging the parquet files according to the categories just like in the data case:

.. code-block:: python

    python prepare_output_file.py --input /absolute/input/path --cats --catDict /absolute/path/to/cat_mc.json --varDict /absolute/path/to/varDict_mc.json --syst --merge --output /absolute/output/path --condor

In order to convert the parquet files to ROOT, one executes:

.. code-block:: python

    python prepare_output_file.py --input /absolute/input_path/to_folder_with_merged --cats --catDict /absolute/path/to/cat_mc.json --varDict /absolute/path/to/varDict_mc.json --syst --root --output /absolute/input_path/to_folder_with_merged --condor

It can happen that if the samples are stored on or the conda environment is installed on `afs`, that the postprocessing can trigger the AFS throttling of your user account which leads to a much slower postprocessing of the files. Usually the throttle is revoked after a few hours.

To avoid that, one has to throttle condor to such an extent, that the AFS throttling is not triggered. This can be done with the `max_materialize` function of HTCondor. This option limits the number of simultaneous processed condor jobs.

An **important** point is that whenever one uses `max_materialize` while using files on `eos`, all condor log, err, out, sub, and sh files have to be put to `afs`. The reason for that is, that the option switches the submission host from lxplus directly to the schedd, which doesn't have access to `eos`. One can specify a separate path which is hosting all the sub and sh files with `--condor-logs`. If the condor log, err, and out files are desired (e.g. for debugging purposes) they can be explicitly produced with `--make-condor-logs`.

A valid command would for example be:

.. code-block:: python

    python prepare_output_file.py --input /absolute/input/path --cats --catDict /absolute/path/to/cat_mc.json --varDict /absolute/path/to/varDict_mc.json --syst --merge --output /absolute/output/path --max-materialize 5 --condor-logs /absolute/path/to/condor/logs --make-condor-logs --condor


.. _def-processor:

----------
Processors
----------
Processors are items defined within Coffea where the analysis workflow is described. While a general overview is available in the `Coffea documentation <https://coffeateam.github.io/coffea/concepts.html#coffea-processor>`_, here we will focus on the aspects that are important for HiggsDNA.

Since in Higgs to diphoton analysis there are some operations that are common to every analysis workflow, we wrote a base processor `HggBaseProcessor <https://higgs-dna.readthedocs.io/en/latest/modules/higgs_dna.workflows.html#higgs_dna.workflows.base.HggBaseProcessor>`_ which can be used in many basic analyses. If more complex operations are needed, one can still write a processor that inherits from the base class and redefines the function ``process``. The operations that one can find within ``HggBaseprocessor.process`` are the following:

* application of filters and triggers
* Chained Quantile Regression to correct shower shapes and isolation variables
* photon IdMVA
* diphoton IdMVA
* photon preselection
* event tagging
* application of systematic uncertainties

Write a New Processor
---------------------

There are cases in which the workflows implemented in HiggsDNA are not enough for your studies. In these cases you might need to **write your own processor**. Depending on the scenario, there are different guidelines to do this.

1. **Hgg-like workflow**. In this case your analysis is similar to the one implemented in the Hgg basic processor, but you need to perform other operations on top (e.g. additional cuts, application of NNs, etc.). In order to reduce the amount **repeated code**, what you can do is write a processor that inherits from ``HggBaseProcessor`` and redefine the function ``process_extra``. You can find an example of this in `DYStudiesProcessor <https://higgs-dna.readthedocs.io/en/latest/modules/higgs_dna.workflows.html#higgs_dna.workflows.dystudies.DYStudiesProcessor>`_.

2. **Non Hgg-like workflow**. This is the case in which the operations you need to perform are different from the ones performed in the ``process`` function of ``HggBaseProcess``. In this kind of scenario you can still inherit from ``HggBaseProcessor`` in order to have access to the same attributes, but you also need to rewrite the ``process`` function. An example of this is the `TagAndProbeProcessor <https://higgs-dna.readthedocs.io/en/latest/_modules/higgs_dna/workflows/dystudies.html#TagAndProbeProcessor>`_. In this case, we cannot use the standard workflow since we manipulate objects in a different way (for instance, we have *tag* and *probe* photons instead of lead and sublead and since each item of a pair can be either tag or probe we need to double the number of candidates - this is an operation that we would never do in a standard workflow).

-------
Taggers
-------

------------------------
Systematic Uncertainties
------------------------
