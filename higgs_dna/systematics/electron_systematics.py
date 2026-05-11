from higgs_dna.systematics.EGM_SS_systematics import Scale_EGM, Smearing_EGM


def Electron_Scale_EGM(pt, events, year="2022postEE", is_correction=True, restriction=None):
    """
    Applies the photon pt scale corrections (use on data!) and corresponding uncertainties (on MC!).
    JSONs need to be pulled first with scripts/pull_files.py
    """

    return Scale_EGM(pt, events, year, is_correction, restriction, is_electron=True)


def Electron_Smearing_EGM(pt, events, year="2022postEE", is_correction=True):
    """
    Applies the photon smearing corrections and corresponding uncertainties (on MC!).
    JSON needs to be pulled first with scripts/pull_files.py
    """

    return Smearing_EGM(pt, events, year, is_correction, is_electron=True)
