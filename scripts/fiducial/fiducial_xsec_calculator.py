import argparse
import awkward as ak
import numpy as np
from scipy import interpolate

available_processes = ['ggH', 'VBFH', 'VH', 'ttH', 'all', 'xH']
available_years = ['2022']
available_eras = ['preEE', 'postEE', 'all']

# Setup command-line argument parsing
parser = argparse.ArgumentParser(description = "Calculate the inclusive fiducial cross section of pp->H(yy)+X process(es) based on processed samples without detector-level selections. ")
parser.add_argument('path', type = str, help = "Path to the top-level folder containing the different directories. Please only run this on the output of the ParticleLevelProcessor.")
parser.add_argument('--process', type = str, choices = available_processes, default = 'ggH', help = "Please specify the process(es) for which you want to calculate the inclusive fiducial xsec.")
parser.add_argument('--year', type = str, choices = available_years, default = '2022', help = 'Please specify the desired year if you want to combine samples from multiple eras.')
parser.add_argument('--era', type = str, choices = available_eras, default = 'postEE', help = "Please specify the era(s) that you want to run over. If you specify 'all', an inverse variance weighting is performed to increase the precision.")
parser.add_argument('--bin', type = str, default = '|0|5000|', help = "Bin boundaries of the differential XS.")
parser.add_argument('--obs', type = str, default = 'PTH', help = "Name of the differential observable.")
parser.add_argument('--powheg', action="store_true", help="To process powheg sample.")

args = parser.parse_args()

path_folder = args.path # Use the specified folder path
# Pepare the processes array appropriately
if args.process == 'all':
    processes = available_processes
    processes.remove('all')
    processes.remove('xH')
elif args.process == 'xH':
    processes = ['VBFH', 'VH', 'ttH']
else:
    processes = [args.process]
year = args.year
# Pepare the eras array appropriately
if args.era == 'all':
    eras = available_eras
    eras.remove('all')
else:
    eras = [args.era]

# Convert bin boundaries in a list
obs_bins = [float(num) for num in args.bin.strip("|").split("|")]

# See also the following pages (note numbers always in picobarn)
# 13: for 125
# 13p6: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap, for 125.38
# 14: for 125
XS_map = {'13':   {'ggH': 48.58, 'VBFH': 3.782, 'VH': 2.2569, 'ttH': 0.5071},
         '13p6': {'ggH': 51.96, 'VBFH': 4.067, 'VH': 2.3781, 'ttH': 0.5638},
         '14':   {'ggH': 54.67, 'VBFH': 4.278, 'VH': 2.4991, 'ttH': 0.6137},}


# This depends on how you named your samples in HiggsDNA
processMap = {'ggH':  'GluGluHtoGG',
              'VBFH': 'VBFHtoGG',
              'VH':   'VHtoGG',
              'ttH':   'ttHtoGG',}

BR = 0.2270/100 # SM value for mH close to 125: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
mass_points = [120,125,130] # The fiducial acceptance should be extrapolated at 125.38 using a spline between 120, 125, and 130
# For powheg only the 125 GeV sample is available
# [FIXME] For the time being, no extrapolation, its effects is in any case small.
# [FIXME] Idea: compute the relative variation between 125 and 125.38 GeV with Madgraph and then apply it to powheg
mass_powheg = {120:125, 125:125, 130:125}

fid_xsecs_per_bin = {}
for b in range(len(obs_bins)-1):
    fid_xsecs_per_bin_process = {}
    for process in processes:
        print(f'INFO: Now extracting fraction of in-fiducial events for process {process} ...')
        in_frac_per_mass = {}
        for mass in mass_points:
            in_frac_per_mass_era = {}
            sumw2_tmp = []
            for era in eras:
                print(f'INFO: Now extracting numbers for era {era}, process {process}, and bin {b} ...')
                # Extract the events
                process_string = path_folder + processMap[process] + '_M-' + str(mass) + '_' + era
                if args.powheg: process_string = path_folder + processMap[process] + '_M-' + str(mass_powheg[mass]) + '_powheg'
                arr = ak.from_parquet(process_string)
                # Calculating the relevant fractions
                inFiducialFlag = (arr.fiducialGeometricFlag == True) & (abs(arr[args.obs]) > obs_bins[b]) & (abs(arr[args.obs]) <= obs_bins[b+1]) # Only for this type of tagger right now, can be customised in the future
                sumwAll = ak.sum(arr.weight)

                sumwIn = ak.sum(arr.weight[(inFiducialFlag)])
                in_frac = sumwIn/sumwAll

                print(f"INFO: Fraction of in-fiducial events: {in_frac} ...")

                sumw2 = ak.sum(arr.weight[(inFiducialFlag)]**2) # This is the MC stat variance
                sumw2_tmp.append(sumw2)

                in_frac_per_mass_era[era] = in_frac * 1/sumw2

            sumw2_tmp = np.asarray(sumw2_tmp)
            result = np.sum(np.asarray([in_frac_per_mass_era[era] for era in eras]))
            in_frac_per_mass[mass] = result / np.sum(1/sumw2_tmp)

        points = [in_frac_per_mass[p] for p in mass_points]
        spline = interpolate.splrep(mass_points, points, k=2)
        fid_xsecs_per_bin_process[process] = float(interpolate.splev(125.38, spline)) * XS_map['13p6'][process] * 1000 * BR

    fid_xsecs_per_bin[b] = np.sum(np.asarray([fid_xsecs_per_bin_process[process] for process in processes]))
    print(f"The fiducial cross section for {args.obs} in [{obs_bins[b]},{obs_bins[b+1]}] is: {fid_xsecs_per_bin[b]} fb")


final_fid_xsec = np.sum(np.asarray([fid_xsecs_per_bin[b] for b in range(len(obs_bins)-1)]))
print(f"The inclusive fiducial cross section is given by: {final_fid_xsec} fb")

output = 'fidXS_'+args.obs+'_'+args.process
if args.powheg: output += '_powheg'

with open(output+'.py', 'w') as f:
        f.write('Boundaries = '+str(obs_bins)+' \n')
        f.write('fidXS = '+str(list(fid_xsecs_per_bin.values()))+' \n')
