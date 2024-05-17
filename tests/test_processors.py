from higgs_dna.utils.logger_utils import setup_logger
from higgs_dna.workflows import DYStudiesProcessor, TagAndProbeProcessor, ZmmyProcessor
import os
import subprocess
from coffea import processor
import json
from importlib import resources


def test_processors():
    """
    Check that both the standard base and the TagAndProbeProcessor can run over a basic nanoAOD data v11 file and no errors occur.
    """
    # Here, we also pull the golden JSON, probably this should be done in a separate unit test in given time...
    # But okay, we need to reevaluate this after the switch to coffea 2023 anyhow
    os.chdir("scripts")
    subprocess.run("python pull_files.py --target GoldenJSON", shell=True)
    os.chdir("..")

    fileset = {
        "Data_2022F": [
            "./tests/samples/skimmed_nano/EGamma_F_Skim.root"
        ]
    }

    with resources.open_text("higgs_dna.metaconditions", "Era2017_legacy_v1.json") as f:
        metaconditions = json.load(f)

    processor_instance = DYStudiesProcessor(
        year={"Data_2022F": ["2022postEE"]},
        metaconditions=metaconditions,
        apply_trigger=True,
        skipCQR=True,
        skipJetVetoMap=True,
        output_location="output/basics"
    )

    iterative_run = processor.Runner(
        executor = processor.IterativeExecutor(compression=None),
        schema=processor.NanoAODSchema,
    )

    out = iterative_run(
        fileset,
        treename="Events",
        processor_instance=processor_instance,
    )

    processor_instance = TagAndProbeProcessor(
        year={"Data_2022F": ["2022postEE"]},
        metaconditions=metaconditions,
        apply_trigger=True,
        skipCQR=True,
        skipJetVetoMap=True,
        output_location="output/basics"
    )

    iterative_run = processor.Runner(
        executor = processor.IterativeExecutor(compression=None),
        schema=processor.NanoAODSchema,
    )

    out = iterative_run(
        fileset,
        treename="Events",
        processor_instance=processor_instance,
    )
    
    processor_instance = ZmmyProcessor(
        year={"Data_2022F": ["2022postEE"]},
        metaconditions=metaconditions,
        apply_trigger=True,
        skipCQR=True,
        skipJetVetoMap=True,
        output_location="output/basics"
    )

    iterative_run = processor.Runner(
        executor = processor.IterativeExecutor(compression=None),
        schema=processor.NanoAODSchema,
    )

    out = iterative_run(
        fileset,
        treename="Events",
        processor_instance=processor_instance,
    )
    
    # clean up
    subprocess.run("rm -r output", shell=True)
