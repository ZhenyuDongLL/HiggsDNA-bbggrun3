import os
import yaml
import awkward as ak
import argparse
import json
from higgs_dna.utils.logger_utils import setup_logger
logger = setup_logger(level='INFO')

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", required=True, type=str, help="path to json mapping of samples")
parser.add_argument('-o', '--output', type=str, dest='output', default="./plots/", help='output folder to plots (created if not existent).')
parser.add_argument("--gjets", action = 'store', nargs = '*', default = [], help = 'GJet Samples to use')
parser.add_argument("--qcd", action = 'store', nargs = '*', default = [], help = 'QCD Samples to use')
parser.add_argument("--ggjets", action = 'store', nargs = '*', default = ['GGJets'], help = 'GGJet Samples to use')
parser.add_argument("--other", action = 'store', nargs = '*', default = [], help = 'Single Higgs and Prompt Photon (non-GGJets) Samples to use')
parser.add_argument("--data", action = 'store', nargs = '*', default = [], help = 'Data Samples to use')
parser.add_argument("--ddqcdgjets", action = 'store', nargs = '*', default = [], help = 'Data Driven Nonprompt Photon Samples to use')
parser.add_argument("--modify", action = 'store', default = 'weight_tot', help = 'Modify weight col (or new weight col) or yaml xs')
parser.add_argument("--new-sample", action = 'store_true', help = 'Generate new DD sample')
parser.add_argument(
        "--run-era",
        action = 'store',
        default = 'Run3_2022postEE',
        choices = ['Run2_2016APV', 'Run2_2016nonAPV', 'Run2_2017', 'Run2_2018', 'Run3_2022preEE','Run3_2022postEE', 'Run3_2023preBPix', 'Run3_2023postBPix', 'Run3_2024'],
        help="Choose PDF from run era"
        )
args = parser.parse_args()

os.makedirs(f'{args.output}/InitialFiles/', exist_ok = True) 
if args.new_sample:
    command  = f'python scripts/DataDrivenBkgEstForNonPrompt/DD_Estimate.py '
    command += f'-i {args.input} '
    command += f'-o "{args.output}/InitialFiles/" '
    command += f'--run-era {args.run_era} '
    os.system(command)
    args.ddqcdgjets = ['DDQCDGJets']
else:
    logger.info('Will not generate a new DD file. Using input yaml file as is')
    os.system(f'cp ' + args.input + ' ' + args.output + '/InitialFiles/' + args.input.split('/')[-1].replace('.yaml', '_updatedDDQCDGJets.yaml'))

with open(args.output + '/InitialFiles/' + args.input.split('/')[-1].replace('.yaml', '_updatedDDQCDGJets.yaml'), 'r') as f:
    files = yaml.load(f, Loader = yaml.Loader)
if len(args.data) == 0:
    args.data = [d for d in files if 'Data' in d]

os.makedirs(f'{args.output}/TFPrepFiles/', exist_ok = True)
logger.info(f'Created {args.output}/TFPrepFiles/')
samples = {'data':args.data, 'gjets': args.gjets, 'qcd':args.qcd, 'ggjets':args.ggjets, 'other':args.other, 'ddqcdgjets':args.ddqcdgjets}
mapping = {}
command  = f'python scripts/DataDrivenBkgEstForNonPrompt/DD_TemplateFitPrep.py '
command += '-i ' + args.output + '/InitialFiles/' + args.input.split('/')[-1].replace('.yaml', '_updatedDDQCDGJets.yaml') + ' '
command += f'-o "{args.output}/TFPrepFiles/" '
for cat in samples:
    if len(samples[cat]) > 0:
        command += f'--{cat} '
        for s in samples[cat]:
            command += s + ' '
os.system(command)

command  = f'python scripts/DataDrivenBkgEstForNonPrompt/DD_TemplateFit.py '
if len(args.qcd) > 0 and len(args.gjets) > 0:
    command += '--useMC '
command += f'--inputDir {args.output}/TFPrepFiles/files/ '
command += f'--output_json --output_latex '
os.system(command)

with open(args.output+'/TFPrepFiles/files/scales.json', 'r') as f:
    scales = json.load(f)

newfiles = dict(files)

os.makedirs(f'{args.output}/reweighted_samples/', exist_ok = True)
for s in files:
    if s in args.ggjets:
        sf = scales['GGJets']['final_sf']
    elif s in args.gjets:
        sf = scales['GJets']['final_sf']
    elif s in args.qcd:
        sf = scales['QCD']['final_sf']
    elif s in args.ddqcdgjets:
        sf = scales['DDQCDGJets']['final_sf']
    else:
        continue
    if args.modify in ['weight', 'weight_tot']:
        logger.info(f'Adjusting weights of {s}...')
        df = ak.from_parquet(files[s]['path'])
        df[args.modify] = df.weight*sf
        logger.info('Saving...')
        ak.to_parquet(df, f'{args.output}/reweighted_samples/{s}_Rescaled.parquet')
        newfiles[f'{s}_Unscaled'] = {'path':files[s]['path'], 'xs':files[s]['xs']}
        newfiles[f'{s}'] = {'path':f'{args.output}/reweighted_samples/{s}_Rescaled.parquet', 'xs':files[s]['xs']}
        del df
    else:
        logger.info(f'Scales applied to xs in sample yaml, not weights of {s}...')
        newfiles[f'{s}_Unscaled'] = {'path':files[s]['path'], 'xs':files[s]['xs']}
        newfiles[f'{s}'] = {'path':files[s]['path'], 'xs':files[s]['xs']*sf}

with open(f'{args.output}/reweighted_samples/{args.input.split("/")[-1].replace(".yaml", "_Rescaled.yaml")}', 'w') as f:
    yaml.dump(newfiles, f, Dumper = yaml.Dumper) 
