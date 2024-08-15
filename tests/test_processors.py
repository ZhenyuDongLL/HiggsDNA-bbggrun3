import os
import subprocess
import json
import pytest
from importlib import resources
from higgs_dna.workflows import DYStudiesProcessor, TagAndProbeProcessor, HHbbggProcessor, HplusCharmProcessor, lowmassProcessor, ParticleLevelProcessor, TopProcessor, ZeeProcessor, ZmmyProcessor
from coffea import processor


# Not tested by default since does not start with "test_"
def run_processor(processor_instance, fileset):
    """
    Helper function to run a given processor instance on the provided fileset.
    """
    iterative_run = processor.Runner(
        executor=processor.IterativeExecutor(compression=None),
        schema=processor.NanoAODSchema,
    )
    out = iterative_run(
        fileset,
        treename="Events",
        processor_instance=processor_instance,
    )
    return out


@pytest.mark.parametrize("processor_class", [
    DYStudiesProcessor,
    TagAndProbeProcessor,
    HHbbggProcessor,
    # Hpc Cannot be included in a simple way here since the arguments are not defaulted
    #HplusCharmProcessor, 
    # Unclear to me why low mass does not work here, unit test should also be designed for this processor
    #lowmassProcessor,
    ParticleLevelProcessor,
    TopProcessor,
    ZeeProcessor,
    #ZmmyProcessor
])
def test_processors(processor_class):
    """
    Test that each processor can run over a basic nanoAOD data v11 file without errors.
    """
    # Pull the golden JSON file
    os.chdir("scripts")
    subprocess.run("python pull_files.py --target GoldenJSON", shell=True)
    os.chdir("..")
    
    # Choose datasets to run over appropriately
    # In the future, should specify datasets on eos instead of local files
    # These should be appropriate for the processor being tested (e.g. muon for Zmmy or DY for T&P)
    # The skeleton below should be adjusted

    if processor_class == DYStudiesProcessor:
        MC = "./tests/samples/skimmed_nano/ggH_M125_amcatnlo_v13.root"
        Data = "./tests/samples/skimmed_nano/EGamma_2022E_v13.root"
    elif processor_class == TagAndProbeProcessor:
        MC = None
        Data = "./tests/samples/skimmed_nano/EGamma_2022E_v13.root"
    elif processor_class == HHbbggProcessor:
        MC = "./tests/samples/skimmed_nano/ggH_M125_amcatnlo_v13.root"
        Data = "./tests/samples/skimmed_nano/EGamma_2022E_v13.root"
    elif processor_class == HplusCharmProcessor:
        MC = "./tests/samples/skimmed_nano/ggH_M125_amcatnlo_v13.root"
        Data = "./tests/samples/skimmed_nano/EGamma_2022E_v13.root"
    elif processor_class == lowmassProcessor:
        MC = "./tests/samples/skimmed_nano/ggH_M125_amcatnlo_v13.root"
        Data = "./tests/samples/skimmed_nano/EGamma_2022E_v13.root"
    elif processor_class == ParticleLevelProcessor:
        MC = "./tests/samples/skimmed_nano/ggH_M125_amcatnlo_v13.root"
        Data = None
    elif processor_class == TopProcessor:
        MC = "./tests/samples/skimmed_nano/ggH_M125_amcatnlo_v13.root"
        Data = "./tests/samples/skimmed_nano/EGamma_2022E_v13.root"
    elif processor_class == ZeeProcessor:
        MC = "./tests/samples/skimmed_nano/ggH_M125_amcatnlo_v13.root"
        Data = "./tests/samples/skimmed_nano/EGamma_2022E_v13.root"
    elif processor_class == ZmmyProcessor:
        MC = None
        Data = None

    fileset = {}

    if Data is not None:
        fileset["Data"] = [Data]

    if MC is not None:
        fileset["MC"] = [MC]

    with resources.open_text("higgs_dna.metaconditions", "Era2017_legacy_v1.json") as f:
        metaconditions = json.load(f)

    processor_instance = processor_class(
        year={"Data": ["2022postEE"], "MC": ["2022postEE"]},
        metaconditions=metaconditions,
        apply_trigger=True,
        skipCQR=True,
        skipJetVetoMap=True,
        output_location="output/basics"
    )
    
    # Run the processor and verify output
    run_processor(processor_instance, fileset)
    
    # Clean up
    subprocess.run("rm -r output", shell=True)
