
from higgs_dna.tools.decorrelator import cdfCalc
import higgs_dna.tools.decorrelator as decorr
import concurrent.futures
import awkward as ak
import pandas as pd
import numpy as np
import argparse
import glob


def read_parquet_file(file_path):
    return pd.read_parquet(file_path)


def printProgressBar(iteration,total,prefix='',suffix='',decimals=1,length=100,fill=chr(9608),printEnd="\r"):

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
    if iteration == total:
        print()


def diphoton_ak_array(diphotons: ak.Array) -> ak.Array:

    output = {}
    for field in ak.fields(diphotons):
        output[field] = diphotons[field]
    return ak.Array(output)


def getArrayBranchName(branchname, fieldname, index):
    if index != ():
        return '{}{}'.format(branchname, index[0])
    return '{}'.format(branchname)


def main(options):
    # reading the parquet files
    if not options.infilepath.endswith('/'):
        print(f"WARNING: Please make sure that {options.infilepath} is a path to .parquet files and ends with /")
        options.infilepath = options.infilepath + "/"
        print(f"INFO: To help you out, the path is changed to: {options.infilepath}")

    df = pd.DataFrame()
    files = glob.glob(str(options.infilepath) + "*.parquet")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        data = list(executor.map(read_parquet_file, files))
    events = pd.concat(data, ignore_index=True)
    print(f"INFO: found {len(events)} events")

    # helper and parquet column read in
    if options.var in events:
        df[options.var] = events[options.var].to_numpy()
    else:
        print(f"ERROR: var not in columns of parquet files from {options.infilepath}")
        print(f"Please choose e.g. from {[col for col in events.columns if col.startswith('sigma')]} or inspect columns for more information!")
        exit()
    df[options.dVar] = events[options.dVar].to_numpy()
    df["weight"] = events.weight.to_numpy()

    # handling of decorrelator
    calc = cdfCalc(df, options.var, options.dVar, np.linspace(100, 180, 161))
    calc.calcCdfs()
    cdfs = calc.cdfs
    dummyDf = pd.DataFrame({'{}'.format(options.var): [0], '{}'.format(options.dVar): [0]})
    decl = decorr.decorrelator(dummyDf, options.var, options.dVar, np.linspace(100., 180., 161))
    decl.cdfs = cdfs
    print("var, dVar:", options.var,", ", options.dVar)

    # giving variables to dataframe and resetting index
    decl.df = df.loc[:, [options.var, options.dVar]]
    decl.df.reset_index(inplace=True)

    # doing the decorr
    decorrelated_var = decl.doDecorr(options.ref)
    df['{}_decorr'.format(options.var)] = decorrelated_var

    # create pickle file and remove the old sigma_m_over_m variable
    df = df.drop(columns=options.var)
    calc = cdfCalc(df, '{}_decorr'.format(options.var), options.dVar, np.linspace(100, 180, 161))
    print(f"INFO: decorrelated CDF contains: {df.columns.tolist()}")
    if options.era:
        if options.outFile:
            calc.dumpCdfs(str(options.outFile) + "_" + str(options.era) + ".pkl.gz")
        else:
            calc.dumpCdfs(str(options.var) + "_" + str(options.era) + "_CDFs" + ".pkl.gz")
    else:
        if options.outFile:
            calc.dumpCdfs(str(options.outFile) + ".pkl.gz")
        else:
            calc.dumpCdfs(str(options.var) + "_CDFs" + ".pkl.gz")
    print("Created pickle file!")

    # a parquet file with the decorrelated variable will be created by default but takes more time
    if not options.parquetGenerationOff:
        print("INFO: new parquet will be created, this can take some time!")
        events['{}_decorr'.format(options.var)] = decorrelated_var
        if options.era:
            if options.outFile:
                events.to_parquet(options.outFile + "_" + str(options.era) + ".parquet")
            else:
                events.to_parquet(str(options.var) + "_decorr_" + str(options.era) + ".parquet")
        else:
            if options.outFile:
                events.to_parquet(options.outFile + ".parquet")
            else:
                events.to_parquet(str(options.var) + "_decorr" + ".parquet")

        print("Created parquet file!")
    else:
        print("INFO: new parquet file will not be created, only pickle!")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group('Required Arguements')
    requiredArgs.add_argument('-i', '--infilepath', action='store', type=str, required=True, help="Input file path, e.g. /net/.../samples/")
    requiredArgs.add_argument('-v','--var', default='sigma_m_over_m', action='store', type=str, help="variable you want to decorrelate (default: sigma_m_over_m)")
    requiredArgs.add_argument('-d','--dVar', default='mass', action='store', type=str, help="variable you want to correlate against, most likely mass, (default: mass)")
    requiredArgs.add_argument('-o','--outFile', action='store', type=str, help="filename and path to the decorrelated files, default .<var>(_era)_CDFs")
    optArgs = parser.add_argument_group('Optional Arguments')
    optArgs.add_argument('-r', '--ref', action='store', type=float, default=125., help="reference mass for decorrelation")
    optArgs.add_argument('--columns', nargs='+')
    optArgs.add_argument('--nomColumns', nargs='+')
    optArgs.add_argument('--vecColumns', nargs='+')
    optArgs.add_argument('-p', '--parquetGenerationOff', action='store_true', help="Set this flag to avoid generating parquet files (by default parquet files are created)")
    optArgs.add_argument('-e','--era', action='store', type=str, help="optional: choose era to give extra information in file name")
    options = parser.parse_args()
    main(options)
