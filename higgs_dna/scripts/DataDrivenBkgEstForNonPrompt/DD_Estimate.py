import awkward
import argparse
import numpy as np
from tqdm import tqdm
import yaml

from higgs_dna.utils.logger_utils import setup_logger
logger = setup_logger(level='INFO')

#############################################################################################################
#  A parquet file containing the events from the data-driven method can be produced with the command:       #
#                                                                                                           #
#       ipython3 scripts/DD_Estimate.py -- -i ${INPUTFILE} -o ${OUTPATH} --run_era ${RUN_ERA}               #
#                                                                                                           #
#  For the run era, choose a run era for which there is a PDF available (Currently 2022preEE and            #
#  2022postEE). For the input file, please use a yaml file formatted like the one linked below:             #
#                                                                                                           #
#       https://gitlab.cern.ch/hhbbgg/docs/-/blob/master/v1/samples_v1_Run3_2022postEE.yaml?ref_type=heads  #
#                                                                                                           #
#############################################################################################################

pdf_coeffs = {
    "Run2_2016APV": {
        "EE": [0.30180157393047863, 0.06046058022684358, 0.6425271397345449, -0.80071245021483, 0.20661911373221137, -0.9207559783395367, 0.5712345931355396, 0.06522299737945454, 0.001716813801638239], 
        "EB": [0.3060731194049686, -0.5354497800785627, 0.10835630252503482, 1.1547922666765857, 1.5619532920135284, -2.4234902779928773]
    },
    "Run2_2016nonAPV": {
        "EE": [0.26333419890100884, -0.3851785964002431, 0.6055964735852113, 0.33896010224283535, -0.04450967234136684, -0.7254013780862512, 1.7014149792468098, -0.8442990529163116, -1.2316622242885247, -0.26627842863243334, 1.3403590467825381, -0.45493260686042625, 0.04165412020931958], 
        "EB": [0.2792311554315338, -0.10691072010119501, 0.6542792746258072, -0.6852063885764395, 0.47894360678677994, -0.6480239607339634, 0.28269124264981554]
    },
    "Run2_2017": {
        "EE": [0.3101944030130834, -0.17921149844643805, 0.27395636269569773, -0.39021570983195814, 1.2159512483679937, -1.0689452763467577], 
        "EB": [0.2848828101651926, -0.21607248878354288, 0.667406586045607, -0.5962994961736051, 0.4804150100498138, -0.39441667009162124, 0.11293366533200291]
        },
    "Run2_2018": {
        "EE": [0.2656859587253099, 0.035874183812668506, 0.6744375979819817, -0.7672881364247405, 0.37160541983332296, -1.029707620922217, 0.6215468846434865], 
        "EB": [0.34405417655668535, -0.3897587728576553, -0.11261425628320787, 0.0800512175717314, 2.457232102249345, -1.002170274216002, -1.2407724790580794]
        },
    "Run3_2022preEE":{
        "EE": [0.3745739369285325, -0.26155191192118127, 0.11398241596508048, 0.5117120026643069, 0.5610293197232438, -1.3923778516490706, 0.6717545764886738, 0.42480298623979806, -0.4740730156553807, -1.1744532213062604, 0.9311427400226765],
        "EB": [0.3902701455451318, -0.5115656696175526, -0.002384306796159553, 1.0900049177717233, 1.1879835363688696, -2.089155161196582]
        },
    "Run3_2022postEE":{
        "EE": [0.39422362871848493, -0.37434203553043915, -0.5389902158493887, 0.8262942159487295, 3.0630025721767993, -2.7858422084474195, -3.660114997786721, 3.6462932432417685, 6.799674103365295, -4.764110799261395, -6.583793780495415, 2.1130404204293103, 2.12283155450363],
        "EB": [0.39259062005763345, -0.4095276669800686, 0.1644570327734382, 1.021480044884044, 0.1399577902061723, -1.931989490230088, 1.6157883788638956, 0.7978700324002087, -1.9115037823819594, -1.443304761732576, 1.7363122773450523]
        },
    "Run3_2023preBPix": {
        "EE": [0.316730817142084, -0.102710340356223, 0.6356405298491293, -0.6222730619633555, 0.31563032692520193, -0.50263552870272, 0.1352738268236311, 0.0652375506779517, 0.004974446650253304], 
        "EB": [0.3413056589740043, -0.12493312907871527, 0.5435991598875018, -0.5668722758763954, 0.2933636923537276, -0.3660036653923755, 0.16144309930473852]
        },
    "Run3_2023postBPix": {
        "EE": [0.33095307237020616, -0.08414531973307297, 0.5640082154085979, -0.5920181994923355, 0.3019912821576091, -0.43973103711917755, 0.2050434851083392], 
        "EB": [0.2920162415649736, -0.0017904857572657112, 0.639274052047168, -0.6550990282824134, 0.3518697069814657, -0.7149142545158571, 0.3652115427244586]
        },
    "Run3_2024": {
        "EE": [0.3868910223581321, -0.2780911329648896, 0.23307364718402426, 0.6735199879641455, 0.2558582280188992, -1.2923080898533055, 0.9342460071561383, 0.1729546264373844, -1.068777613110097, -0.8318512305835933, 1.071297867351851], 
        "EB": [0.3846668584572264, -0.341115193269897, 0.2801691799353723, 0.641982165526039, 0.14760708449685261, -1.4125244150759988, 1.092402243543304, 0.04920663248185148, -1.2138529287888213, -0.8160207691141318, 1.291894362578286]
        }   
}

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", 
        "-i", 
        required=True, 
        type=str, 
        help="path to yaml mapping of samples"
        )

    parser.add_argument(
        "--output", 
        "-o", 
        default = './DDQCD_Output/', 
        type=str, 
        help="where to output files"
        )
 
    parser.add_argument(
        "--run-era",
        action = 'store',
        default = 'Run3_2022postEE',
        choices = ['Run2_2016APV', 'Run2_2016nonAPV', 'Run2_2017', 'Run2_2018', 'Run3_2022preEE','Run3_2022postEE', 'Run3_2023preBPix', 'Run3_2023postBPix', 'Run3_2024'],
        help="Choose PDF from run era"
        )

    parser.add_argument(
        "--unrestricted",
        action = 'store_true',
        help="Use alternate sampling method"
        )


    return parser.parse_args()


def sample(pdfcoeffs, maxMVAIDs, idcut, args):
    ppf_res = 10000

    # Define PPF(Inverse of CDF) for Sampling
    pdf = np.polynomial.Polynomial(pdfcoeffs)
    maxMVAIDs = np.array(maxMVAIDs)
    cdf = pdf.integ(lbnd = -0.9)
    denom = cdf(idcut)
    maxodds = cdf(maxMVAIDs)
    minodds = cdf(idcut)

    y = np.linspace(-0.9, 1, ppf_res)
    ppfdist = cdf(y)
    ppf = lambda x: np.interp(x, ppfdist, y)

    #Initial min MVA IDs Sample
    probs = np.random.random(size = len(maxMVAIDs))
    probs = ((maxodds - minodds)*probs) + minodds
    minMVAIDs = ppf(probs)
    #Resample until min MVA IDs satisfy requirements, truncate if too many attempts used
    #Should not engage, exists as a failsafe
    if args.unrestricted:
        idx = (minMVAIDs > 1) | (minMVAIDs < idcut)
    else:
        idx = (minMVAIDs > maxMVAIDs) | (minMVAIDs < idcut)
    if any(idx):
        N = 100000
        for att in tqdm(range(N)):
            if not any(idx):
                break
            probs[idx] = np.random.random(size = np.sum(idx))
            minMVAIDs = ppf(probs)
            if args.unrestricted:
                idx = (minMVAIDs > 1) | (minMVAIDs < idcut)
            else:
                idx = (minMVAIDs > maxMVAIDs) | (minMVAIDs < idcut) 
        if att + 1 == N:
            logger.info('Attempts maxed out.')
            logger.info(str(np.sum(minMVAIDs > maxMVAIDs)) + 'unweighted events failed the max MVA ID cut.')
            logger.info(str(np.sum(minMVAIDs < idcut)) + 'unweighted events failed the preselection cut.')

    
    # set weights based on max MVA ID, PDF
    if args.unrestricted:
        weights = np.ones_like(maxMVAIDs)
        newmin = np.min([minMVAIDs, maxMVAIDs], axis = 0)
        newmax = np.max([minMVAIDs, maxMVAIDs], axis = 0)
        minMVAIDs = newmin
        maxMVAIDs = newmax
    weights = (cdf(maxMVAIDs) - denom)/denom

    return minMVAIDs, maxMVAIDs, weights

def loadData(files, proc, min_cols = True):
    #Load Sample and apply correct lumi, xs, and create new useful variables
    md = files[proc]

    #If min_cols, use minimum necessary columns to save memory
    if min_cols:
        cols = [
                'weight',
                'lead_mvaID',
                'lead_isScEtaEE',
                'lead_isScEtaEB',
                'sublead_mvaID',
                'sublead_isScEtaEE',
                'sublead_isScEtaEB',
                'n_jets',
                'mass'
                ]
    else:
        cols = None

    #Load Sample File
    path = md['path']
    logger.info('Loading ' + path)
    df = awkward.from_parquet(path, columns = cols)
    #Adding necessary variables
    logger.info('Adding Min/Max MVAID Variables')
    df['Max_mvaID'] = np.max([df.lead_mvaID, df.sublead_mvaID], axis = 0)
    df['Min_mvaID'] = np.min([df.lead_mvaID, df.sublead_mvaID], axis = 0)
    idx = df.lead_mvaID == df.Max_mvaID
    for var in ['isScEtaEB', 'isScEtaEE']:
        df['Max_mvaID_'+var] = awkward.where(
            idx,
            df['lead_'+var],
            df['sublead_'+var],
        )
        df['Min_mvaID_'+var] = awkward.where(
            ~idx,
            df['lead_'+var],
            df['sublead_'+var],
        )
    return df

def main(args):
    
    #Load File Mapping
    with open(args.input, 'r') as f:
        files = yaml.load(f, Loader = yaml.Loader)

    #Load appropriate PDF Coeffs
    coeffs = pdf_coeffs[args.run_era]
    coeffs_EE = coeffs['EE']
    coeffs_EB = coeffs['EB']
    idcut = -0.7


    #Get sideband data events
    for i, proc in enumerate([f for f in files if 'Data' in f]):
        if i == 0:
            data = loadData(files, proc)
        else:
            temp = loadData(files, proc)
            data = awkward.concatenate([data, temp])
    sdbd_events = data[(data.Max_mvaID > idcut) & (data.Min_mvaID < idcut)]
    logger.info('Total SBD Yield: ' + str(np.sum(sdbd_events.weight)))
    del data

    #Generate new min MVA IDs using EE, EB PDFs
    Min_mvaID_EE, Max_mvaID_EE, Min_weights_EE = sample(coeffs_EE, sdbd_events.Max_mvaID, idcut, args)
    Min_mvaID_EB, Max_mvaID_EB, Min_weights_EB = sample(coeffs_EB, sdbd_events.Max_mvaID, idcut, args)
    
    #Collect proper mva IDs
    newmin = np.ones_like(Min_mvaID_EB)
    newweight = np.ones_like(Min_mvaID_EB)
    newmin[sdbd_events.Min_mvaID_isScEtaEB] = Min_mvaID_EB[sdbd_events.Min_mvaID_isScEtaEB]
    newmin[sdbd_events.Min_mvaID_isScEtaEE] = Min_mvaID_EE[sdbd_events.Min_mvaID_isScEtaEE]
    newweight[sdbd_events.Min_mvaID_isScEtaEB] = Min_weights_EB[sdbd_events.Min_mvaID_isScEtaEB]
    newweight[sdbd_events.Min_mvaID_isScEtaEE] = Min_weights_EE[sdbd_events.Min_mvaID_isScEtaEE]
    logger.info('DD EE Yields: ' + str(np.sum(newweight[sdbd_events.Min_mvaID_isScEtaEE])))
    logger.info('DD EB Yields: ' + str(np.sum(newweight[sdbd_events.Min_mvaID_isScEtaEB])))

    #Stitch min MVA IDs to lead, sublead photons
    leadismax = sdbd_events.lead_mvaID == sdbd_events.Max_mvaID
    sdbd_events['lead_mvaID'] = awkward.where(
        ~leadismax,
        newmin,
        sdbd_events.lead_mvaID
    )
    sdbd_events['sublead_mvaID'] = awkward.where(
        leadismax,
        newmin,
        sdbd_events.sublead_mvaID
    )
    sdbd_events['weight'] = newweight

    #Drop extra variables and save parquet file
    dropper = [f for f in sdbd_events.fields if 'Min_mvaID' in f and 'Max_mvaID' in f]
    keeper_fields = [f for f in sdbd_events.fields if not f in dropper]
    sdbd_events = sdbd_events[keeper_fields]
    filename = args.output + "/DDQCDGJET.parquet"
    logger.info('Saving ' + filename)
    awkward.to_parquet(sdbd_events, filename)
    files['DDQCDGJets'] = {'path':filename, 'xs':1}
    filename = args.output + args.input.split('/')[-1].replace('.yaml', '_updatedDDQCDGJets.yaml')
    logger.info('Saving ' + filename)
    with open(filename, 'w') as f:
        yaml.dump(files, f)
    return sdbd_events
    

if __name__ == "__main__":
    args = parse_arguments()
    df = main(args)
