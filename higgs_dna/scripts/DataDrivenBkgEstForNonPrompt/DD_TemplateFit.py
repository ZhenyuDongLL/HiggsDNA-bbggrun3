import os, sys, copy
import ROOT
import numpy
import ctypes
import argparse
import json
from higgs_dna.utils.logger_utils import setup_logger
logger = setup_logger(level='INFO')

parser = argparse.ArgumentParser()
parser.add_argument("--inputDir", required=True, help = "input directory containing root files with histograms to fit", type = str)
parser.add_argument("--useMC", help = "perform fits using MC nonprompt photon samples", action = 'store_true', default = False)
parser.add_argument("--output_latex", help = "Store test file containing Template Fit results in LATEX format", action = 'store_true', default = False)
parser.add_argument("--output_json", help = "Store JSON file with scale factors", action = 'store_true', default = False)
args = parser.parse_args()

def find_jet_bin(jet_bin):
  return jet_bin + 1

def combine_hists(hist1, hist2, name):
  target_hist = ROOT.TH1D(name, "", hist1.GetNbinsX() + hist2.GetNbinsX(), 0, 1)
  jet_bin_number = 3
  for i in range(hist1.GetNbinsX()):
    target_hist.SetBinContent(i+1, hist1.GetBinContent(i+1, jet_bin_number))
    for j in range(jet_bin_number, hist1.GetNbinsY()):
        target_hist.AddBinContent(i+1, hist1.GetBinContent(i+1, j+1))
  for i in range(hist2.GetNbinsX()):
    target_hist.SetBinContent(hist1.GetNbinsX()+i+1, hist2.GetBinContent(i+1, jet_bin_number))
    for j in range(jet_bin_number, hist2.GetNbinsY()):
        target_hist.AddBinContent(hist1.GetNbinsX()+i+1, hist2.GetBinContent(i+1, j+1))
  for i in range(hist1.GetNbinsX() + hist2.GetNbinsX()):
    if target_hist.GetBinContent(i+1) < 0:
      target_hist.SetBinContent(i+1, 0)

  return target_hist


latex_dict = { 
        "GGJets"     : "$ \\gamma \\gamma $ + jets",
        "GJets"        : "$ \\gamma $ + jets",
        "QCD"          : "QCD",
        "DDQCDGJets"   : "QCD, $ \\gamma $ + jets",
        "Other_Bkgs"   : "Other Backgrounds"
        }

file_names = ["n_jets:Min_mvaID", "n_jets:Max_mvaID"]

templates_dict = { 
        "data" : "Data",
        "pp"   : "GGJets",
        "fp"   : "GJets",
        "ff"   : "QCD",
		"fffp" : "DDQCDGJets",
        "bkg"  : "Other_Bkgs"
        }

hists = {}
hists_unweighted = {}; hists_unweighted_minMVAID = {}; hists_unweighted_maxMVAID = {}
hists_weighted = {}; hists_weighted_minMVAID = {}; hists_weighted_maxMVAID = {}
hists_weights = {}

for file_name in file_names:
  for template in templates_dict:
    name = templates_dict[template]
    try:
        logger.info('Loading ' + args.inputDir+"/"+name+'_'+file_name+"_unweighted.root")
        f_unweighted = ROOT.TFile(args.inputDir+"/"+name+'_'+file_name+"_unweighted.root")
    except:
        logger.info('File for ' + name + ' not found. Continuing...')
        continue
    if not f_unweighted or f_unweighted.IsZombie():
        logger.info('Issue loading file for ' + name + ' not found. Continuing...')
        continue

    if "min" in file_name.lower(): hists_unweighted_minMVAID[template] = copy.deepcopy(f_unweighted.Get("totalSM"))
    elif "max" in file_name.lower(): hists_unweighted_maxMVAID[template] = copy.deepcopy(f_unweighted.Get("totalSM"))

    f_weighted = ROOT.TFile(args.inputDir+"/"+name+'_'+file_name+"_weighted.root")
    if "min" in file_name.lower(): hists_weighted_minMVAID[template] = copy.deepcopy(f_weighted.Get("totalSM"))
    elif "max" in file_name.lower(): hists_weighted_maxMVAID[template] = copy.deepcopy(f_weighted.Get("totalSM"))

for template in hists_weighted_minMVAID:
  name = templates_dict[template]
  hists_unweighted[template] = combine_hists(hists_unweighted_minMVAID[template], hists_unweighted_maxMVAID[template], template + "unweighted")
  hists_weighted[template] = combine_hists(hists_weighted_minMVAID[template], hists_weighted_maxMVAID[template], template + "weighted")
  hists_weights[template] = hists_weighted[template].Clone(template + "weights")
  hists_weights[template].Divide(hists_unweighted[template])
  for i in range(hists_weights[template].GetNbinsX()):
    if hists_weights[template].GetBinContent(i+1) <= 0:
       hists_weights[template].SetBinContent(i+1, 0.000000001)
  hists[template] = [ hists_unweighted[template], hists_weighted[template], hists_weights[template] ]

DDavail = ("fffp" in hists) and ("pp" in hists)
MCavail = ("ff" in hists) and ("fp" in hists) and ("pp" in hists)
if MCavail and args.useMC:
    logger.info('useMC requested and all samples available. Performing template fit on MC...')
    use = ["ff", "fp", "pp"]
elif DDavail and args.useMC:
    logger.info('useMC requested but insufficient samples available. Performing template fit using DD Estimate...')
    use = ["fffp", "pp"]
elif DDavail:
    logger.info('Performing template fit using DD Estimate...')
    use = ["fffp", "pp"]
else:
    logger.info('Insufficient Samples to perform template fit using DD Estimate. Exiting...')
    print([f for f in hists])
    asdf
if "bkg" in hists:
    use.append("bkg")

mc = ROOT.TObjArray(len(use))
hist_unweighted_list = []; hist_weighted_list = []; hist_weights_list = []

for cat in use:
    hist_unweighted_list.append(hists[cat][0])
    hist_weighted_list.append(hists[cat][1])
    hist_weights_list.append(hists[cat][2])
    mc.Add(hists[cat][0])

h_data =  hists["data"][1]
initial_fracs = []
for hist in hist_weighted_list:
  initial_fracs.append(hist.Integral() / h_data.Integral())

fit = ROOT.TFractionFitter(h_data, mc, "V")#, "Q")

if 'bkg' in use:
  i_bkg = use.index('bkg')
else:
  i_bkg = -1
i_ggjets = use.index('pp')

for i in range(len(use)):
  logger.info(f'Weighting {use[i]}')
  fit.SetWeight(i, hist_weights_list[i]) # set bin-by-bin weights for raw MC counts 
  # Other bkgs, fixed in fit
  if i == i_bkg:
    fit.Constrain(i, initial_fracs[i]-0.000000001, initial_fracs[i]+0.000000001)
    logger.info(f'Freezing {use[i]}')
  else:
    fit.Constrain(i, 0.0, 1.0)

fit.Fit()

fracs = []
frac_errs = []
scales = []
for i in range(len(use)):
  f = ctypes.c_double()
  ferr = ctypes.c_double()
  fit.GetResult(i, f, ferr)
  fracs.append(f)
  frac_errs.append(ferr)
  scales.append(f.value/initial_fracs[i])

fracs = [a.value for a in fracs]
frac_errs= [a.value for a in frac_errs]

# Print fit results
logger.info( f"Included BkGs {use}")
logger.info( f"Included BkGs {[templates_dict[f] for f in use]}")
logger.info( f"Initial fractions {initial_fracs}")
logger.info( f"Fitted fractions {fracs}")
logger.info( f"Errors {frac_errs}")
logger.info( f"Scales {scales}")


if args.output_latex:
    with open(f'{args.inputDir}/results_latex.tex', 'w' ) as f:
        # Make pretty table
        f.write( "\\begin{center} \\Fontvi \\begin{tabular}{|l|| r| r| r|} \\hline \n")
        f.write( "Template & Initial Fraction & Fitted Fraction & Scale \\\\ \\hline \n")
        for i in range(len(initial_fracs)):
            template = templates_dict[use[i]]
            f.write( "%s & %.4f & %.4f $ \\pm $ %.4f & %.4f \\\\ \n" % (latex_dict[template], initial_fracs[i], fracs[i], frac_errs[i], scales[i]))
        f.write('\\hline \n')
        f.write( "\\end{tabular} \\end{center}")

if args.output_json:
  table = {}
  for i in range(len(initial_fracs)):
    template = templates_dict[use[i]]
    table[template] = {
        'inital_frac':initial_fracs[i], 
        'final_frac':fracs[i], 
        'frac_unc':frac_errs[i], 
        'final_sf':scales[i]
      }
  with open(f'{args.inputDir}/scales.json', 'w' ) as f:
    json.dump(table, f)
  
