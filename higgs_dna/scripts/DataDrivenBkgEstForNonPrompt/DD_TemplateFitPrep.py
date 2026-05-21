import ROOT
import matplotlib.pyplot as plt

import os
import awkward
import argparse
import yaml
import numpy as np
from higgs_dna.scripts.DataDrivenBkgEstForNonPrompt.DD_Estimate import loadData
from higgs_dna.utils.logger_utils import setup_logger
logger = setup_logger(level='INFO')

plotXTitles = {}
plotXTitles["Max_mvaID"] = "max #gamma MVA ID"
plotXTitles["Min_mvaID"] = "min #gamma MVA ID"
plotXTitles["n_jets:Max_mvaID"] = "max #gamma MVA ID"
plotXTitles["n_jets:Min_mvaID"] = "min #gamma MVA ID"
plotXTitles["n_jets"] = "N_{jets}"


plotYTitles = {}
plotYTitles["n_jets"] = "Events"
plotYTitles["Max_mvaID"] = "Events"
plotYTitles["Min_mvaID"] = "Events"
plotYTitles["n_jets:Max_mvaID"] = "N_{jets}"
plotYTitles["n_jets:Min_mvaID"] = "N_{jets}"

bins = {
        "Max_mvaID":(50, -1, 1),
        "Min_mvaID":(50, -1, 1),
        "n_jets":(10, 0, 10)
        }

def process(new, sample, weighted, output):
    if weighted:
        postfix = f'_weighted.png'
        title = sample + 'Weighted'
    else:
        postfix = '_unweighted.png'
        title = sample + 'Unweighted'
    temp = {}
    
    temp["Max_mvaID"] = new.Max_mvaID.to_numpy()
    temp["Min_mvaID"] = new.Min_mvaID.to_numpy()
    temp["n_jets"] = new.n_jets.to_numpy()
    temp["weight"] = new.weight_tot.to_numpy()
    dftemp = ROOT.RDF.FromNumpy(temp)
    title = 'totalSM'

    v2 = "n_jets"
    for v1 in ["Max_mvaID", "Min_mvaID", "n_jets"]:

        mode = v1
        plotXTitle = plotXTitles[mode]
        plotYTitle = plotYTitles[mode]
        nBinsX,lowBinX,highBinX = bins[mode]
        f = ROOT.TFile(output + '/files/' + sample.replace(' ', '_') + '_' + mode + postfix.replace('png', 'root'),"RECREATE")
        c = ROOT.TCanvas('c')
        if weighted:
            hist = dftemp.Histo1D((title,";"+plotXTitle+";"+plotYTitle,nBinsX,lowBinX,highBinX), v1, 'weight')
        else:
            hist = dftemp.Histo1D((title,";"+plotXTitle+";"+plotYTitle,nBinsX,lowBinX,highBinX), v1)
        hist.Draw()
        c.Update()
        c.Print(output + sample.replace(' ', '_') + '_' + mode + postfix)
        hist.Write()
        f.Close()
        if v1 == v2:
            continue

        mode = v2 + ':' + v1
        plotXTitle = plotXTitles[mode]
        plotYTitle = plotYTitles[mode]
        nBinsX,lowBinX,highBinX = bins[v1]
        nBinsY,lowBinY,highBinY = bins[v2]
        f = ROOT.TFile(output + '/files/' + sample.replace(' ', '_') + '_' + mode + postfix.replace('png', 'root'),"RECREATE")
        c = ROOT.TCanvas('c')
        if weighted:
            c.SetGridx()
            c.SetGridy()
            hist = dftemp.Histo2D((title, ";"+plotXTitle+";"+plotYTitle, nBinsX, lowBinX, highBinX, nBinsY, lowBinY, highBinY), v1, v2, 'weight') 
        else:
            hist = dftemp.Histo2D((title, ";"+plotXTitle+";"+plotYTitle, nBinsX, lowBinX, highBinX, nBinsY, lowBinY, highBinY), v1, v2)
        hist.Draw('COLZ')
        c.Update()
        c.Print(output + sample.replace(' ', '_') + '_' + mode + postfix)
        hist.Write()
        f.Close()

def addWeightsandCuts(files, proc, df):    
    df = df[ (df.mass < 180) & (df.mass > 100) ]
    df = df[(df.lead_mvaID > -0.7)&(df.sublead_mvaID > -0.7)]

    if not 'Data' in proc and not 'DD' in proc:
        lumi_tot = 0
        for s in files:
            if 'Data' in s:
                lumi_tot += files[s]['lumi']
        md = files[proc]
        xs = md['xs']
        df['weight_tot'] = df.weight*xs*lumi_tot
    else:
        df['weight_tot'] = df.weight
    return df

def buildArray(files, procs):
    for i, proc in enumerate(procs):
        if i == 0:
            array = loadData(files, proc, min_cols = True)
            array = addWeightsandCuts(files, proc, array)
        else:
            temp = loadData(files, proc, min_cols = True)
            temp = addWeightsandCuts(files, proc, temp)
            array = awkward.concatenate([array, temp])
    return array

def makeHist(array, weights, xlabel, ylabel, proc, var):
    fig = plt.figure(figsize = (10, 8))
    currbins = [bins[var][1], bins[var][2], bins[var][0] + 1]
    plt.hist(array, weights = weights, bins = np.linspace(*currbins))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(args.output + f'plt_{proc.replace(" ", "_")}_{var}.png')
    plt.close()


def plot(df, proc, args):
    makeHist(
        df.Min_mvaID,
        df.weight,
        r'min[$\gamma$ MVA ID]',
        'Events',
        proc,
        'Min_mvaID',
    )
    makeHist(
        df.Max_mvaID,
        df.weight,
        r'max[$\gamma$ MVA ID]',
        'Events',
        proc,
        'Max_mvaID',
    )

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", required=True, type=str, help="path to json mapping of samples")
parser.add_argument('-o', '--output', type=str, dest='output', default="./plots/", help='output folder to plots (created if not existent).')
parser.add_argument("--gjets", action = 'store', nargs = '*', default = [], help = 'Which GJets samples to prepare')
parser.add_argument("--qcd", action = 'store', nargs = '*', default = [], help = 'Which QCD samples to prepare')
parser.add_argument("--ggjets", action = 'store', nargs = '*', default = [], help = 'Which GGJets samples to prepare')
parser.add_argument("--other", action = 'store', nargs = '*', default = [], help = 'Which Single Higgs and other samples to prepare')
parser.add_argument("--data", action = 'store', nargs = '*', default = [], help = 'Which Data files to prepare')
parser.add_argument("--ddqcdgjets", action = 'store', nargs = '*', default = [], help = 'Which DD Nonprompt Background files to prepare')
args = parser.parse_args()

os.makedirs(args.output + '/files/', exist_ok = True)

with open(args.input, 'r') as f:
    files = yaml.load(f, Loader = yaml.Loader)

bkgs = {
    'GJets':args.gjets,
    'QCD': args.qcd,
    'GGJets': args.ggjets,
    'Data': args.data,
    'Other Bkgs': args.other,
    'DDQCDGJets': args.ddqcdgjets,
}
for proc in bkgs:
    if len(bkgs[proc]) == 0:
        continue
    events = buildArray(files, bkgs[proc])
    plot(events, proc, args)
    process(events, proc, True, args.output)
    process(events, proc, False, args.output)
