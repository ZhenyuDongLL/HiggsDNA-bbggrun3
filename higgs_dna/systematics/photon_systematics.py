import numpy as np
import awkward as ak
import correctionlib
import os
import sys
import logging
from higgs_dna.systematics.EGM_SS_systematics import Scale_EGM, Smearing_EGM, Scale_ZeeZmmg, Smearing_ZeeZmmg

logger = logging.getLogger(__name__)


# first dummy, keeping it at this point as reference for even simpler implementations
def photon_pt_scale_dummy(pt, **kwargs):
    return (1.0 + np.array([0.05, -0.05], dtype=np.float32)) * pt[:, None]


# Not nice but working: if the functions are called in the base processor by Photon.add_systematic(... "what"="pt"...), the pt is passed to the function as first argument.
# I need the full events here, so I pass in addition the events. Seems to only work if it is explicitly a function of pt, but I might be missing something. Open for better solutions.
def Photon_Scale_EGM(pt, events, year="2022postEE", is_correction=True, restriction=None):
    """
    Applies the photon pt scale corrections (use on data!) and corresponding uncertainties (on MC!).
    JSONs need to be pulled first with scripts/pull_files.py
    """

    return Scale_EGM(pt, events, year, is_correction, restriction, is_electron=False)


def Photon_Smearing_EGM(pt, events, year="2022postEE", is_correction=True):
    """
    Applies the photon smearing corrections and corresponding uncertainties (on MC!).
    JSON needs to be pulled first with scripts/pull_files.py
    """

    return Smearing_EGM(pt, events, year, is_correction, is_electron=False)


def Scale(pt, events, year="2022postEE", is_correction=True, gaussians="1G", restriction=None, is_Zee=True):
    """
    Applies the IJazZ photon pt scale corrections (use on data!) and corresponding uncertainties (on MC!).
    JSONs need to be pulled first with scripts/pull_files.py.
    The IJazZ corrections are independent and detached from the Egamma corrections.
    """

    return Scale_ZeeZmmg(pt, events, year, is_correction, gaussians, restriction, is_Zee=is_Zee)


def Smearing(pt, events, year="2022postEE", is_correction=True, gaussians="1G"):
    """
    Applies the photon smearing corrections and corresponding uncertainties (on MC!).
    JSON needs to be pulled first with scripts/pull_files.py
    """

    return Smearing_ZeeZmmg(pt, events, year, is_correction, gaussians)


def energyErrShift(energyErr, events, year="2022postEE", is_correction=True, workflow="base"):
    # See also https://indico.cern.ch/event/1131803/contributions/4758593/attachments/2398621/4111806/Hgg_Differentials_Approval_080322.pdf#page=47
    # 2% with flows justified by https://indico.cern.ch/event/1495536/#20-study-of-the-sigma_mm-mismo
    if is_correction:
        return events
    else:
        _energyErr = ak.flatten(events.Photon.energyErr)
        if workflow == "lowmass":
            # For lowmass: 4% is a preliminary estimate, based on the observed data-MC discrepancies in the energy resolution.
            uncertainty_up = np.ones(len(_energyErr)) * 1.04
            uncertainty_dn = np.ones(len(_energyErr)) * 0.96
        else:
            uncertainty_up = np.ones(len(_energyErr)) * 1.02
            uncertainty_dn = np.ones(len(_energyErr)) * 0.98
        return (
            np.concatenate(
                (uncertainty_up[:, None], uncertainty_dn[:, None]), axis=1
            )
            * _energyErr[:, None]
        )


# Not nice but working: if the functions are called in the base processor by Photon.add_systematic(... "what"="pt"...), the pt is passed to the function as first argument.
# I need the full events here, so I pass in addition the events. Seems to only work if it is explicitly a function of pt, but I might be missing something. Open for better solutions.
def FNUF(pt, events, year="2017", is_correction=True):
    """
    ---This is an implementation of the FNUF uncertainty copied from flashgg,
    --- Preliminary JSON (run2 I don't know if this needs to be changed) file created with correctionlib starting from flashgg: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    Applies the photon pt and energy scale corrections and corresponding uncertainties.
    To be checked by experts
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)
    eta = ak.flatten(abs(events.Photon.ScEta))
    r9 = ak.flatten(events.Photon.r9)
    _pt = ak.flatten(events.Photon.pt)

    # era/year defined as parameter of the function
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
    if year not in avail_years:
        logger.error(f"Only FNUF corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        sys.exit(1)
    elif "2022" in year or "2023" in year or "2024" in year:
        logger.warning(f"""You selected the year_string {year}, which is a Run 3 era.
                        FNUF was not re-derived for Run 3 yet, but we fall back to the Run 2 2018 values.
                        These values only constitute up/down variations, no correction is applied.
                        The values are the averaged corrections from Run 2, turned into a systematic and inflated by 25%.
                        Please make sure that this is what you want. You have been warned.""")
        # The values have been provided by Badder for HIG-23-014 and Fabrice suggested to increase uncertainty a bit.
        year = "2022"
    elif "2016" in year:
        year = "2016"

    jsonpog_file = os.path.join(os.path.dirname(__file__), f"JSONs/FNUF/{year}/FNUF_{year}.json")
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["FNUF"]

    if is_correction:
        correction = evaluator.evaluate("nominal", eta, r9)
        corr_pt = _pt * correction

        corrected_photons = events.Photon
        corr_pt = ak.unflatten(corr_pt, counts)
        corrected_photons["pt"] = corr_pt
        events["Photon"] = corrected_photons

        return events

    else:
        correction = evaluator.evaluate("nominal", eta, r9)
        # When creating the JSON I already included added the variation to the returned value
        uncertainty_up = evaluator.evaluate("up", eta, r9) / correction
        uncertainty_dn = evaluator.evaluate("down", eta, r9) / correction
        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return (
            np.concatenate(
                (uncertainty_up[:, None], uncertainty_dn[:, None]), axis=1
            )
            * _pt[:, None]
        )


# Same old same old, just reiterated: if the functions are called in the base processor by Photon.add_systematic(... "what"="pt"...), the pt is passed to the function as first argument.
# Open for better solutions.
def ShowerShape(pt, events, year="2017", is_correction=True, workflow="base"):
    """
    ---This is an implementation of the ShowerShape uncertainty copied from flashgg,
    --- Preliminary JSON (run2 I don't know if this needs to be changed) file created with correctionlib starting from flashgg: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    Applies the photon pt and energy scale corrections and corresponding uncertainties (only on the pt because it is what is used in selection).
    To be checked by experts
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)
    eta = ak.flatten(events.Photon.ScEta)
    r9 = ak.flatten(events.Photon.r9)
    _pt = ak.flatten(events.Photon.pt)
    if workflow == "lowmass":
        logger.warning(f"""You selected the workflow {workflow}.
                        The ShowerShape systematic is not yet rederived for the Run3 lowmass analysis, but we apply the same systematic of Run2 2018.
                       """)
        year = "2018"
        # add also energy here for lowmass
        _energy = ak.flatten(events.Photon.energy)

    # era/year defined as parameter of the function
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018"]
    if year not in avail_years:
        logger.error(f"Only ShowerShape corrections for the year strings {avail_years} are already implemented! ShowerShape should not be used in run3 eras. \n Exiting. \n")
        raise ValueError(f"Only ShowerShape corrections for the year strings {avail_years} are already implemented! ShowerShape should not be used in run3 eras. \n Exiting. \n")
    elif "2016" in year:
        year = "2016"

    jsonpog_file = os.path.join(os.path.dirname(__file__), f"JSONs/ShowerShape/{year}/ShowerShape_{year}.json")
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["ShowerShape"]

    if is_correction:
        correction = evaluator.evaluate("nominal", eta, r9)
        corr_pt = _pt * correction

        corrected_photons = events.Photon
        corr_pt = ak.unflatten(corr_pt, counts)
        corrected_photons["pt"] = corr_pt
        # add also energy here for lowmass
        if workflow == "lowmass":
            corr_energy = _energy * correction
            corr_energy = ak.unflatten(corr_energy, counts)
            corrected_photons["energy"] = corr_energy
        events["Photon"] = corrected_photons

        return events

    else:
        correction = evaluator.evaluate("nominal", eta, r9)
        # When creating the JSON I already included added the variation to the returned value
        uncertainty_up = evaluator.evaluate("up", eta, r9) / correction
        uncertainty_dn = evaluator.evaluate("down", eta, r9) / correction
        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return (
            np.concatenate(
                (uncertainty_up[:, None], uncertainty_dn[:, None]), axis=1
            )
            * _pt[:, None]
        )


def Material(pt, events, year="2017", is_correction=True):
    """
    ---This is an implementation of the Material uncertainty copied from flashgg,
    --- JSON file for run2 created with correctionlib starting from flashgg: https://github.com/cms-analysis/flashgg/blob/2677dfea2f0f40980993cade55144636656a8a4f/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py
    Applies the photon pt and energy scale corrections and corresponding uncertainties.
    To be checked by experts
    """

    # for later unflattening:
    counts = ak.num(events.Photon.pt)
    eta = ak.flatten(abs(events.Photon.ScEta))
    r9 = ak.flatten(events.Photon.r9)
    _pt = ak.flatten(events.Photon.pt)

    # era/year defined as parameter of the function, only 2017 is implemented up to now
    avail_years = ["2016", "2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
    if year not in avail_years:
        logger.error(f"Only eVetoSF corrections for the year strings {avail_years} are already implemented! \n Exiting. \n")
        sys.exit(1)
    elif "2016" in year:
        year = "2016"
    # use Run 2 files also for Run 3, preliminary
    elif year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]:
        logger.warning(f"""You selected the year_string {year}, which is a Run 3 era.
                  Material was not rederived for Run 3 yet, but we fall back to the Run 2 2018 values.
                  Please make sure that this is what you want. You have been warned.""")
        year = "2018"

    jsonpog_file = os.path.join(os.path.dirname(__file__), f"JSONs/Material/{year}/Material_{year}.json")
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["Material"]

    if is_correction:
        correction = evaluator.evaluate("nominal", eta, r9)
        corr_pt = _pt * correction

        corrected_photons = events.Photon
        corr_pt = ak.unflatten(corr_pt, counts)
        corrected_photons["pt"] = corr_pt
        events["Photon"] = corrected_photons

        return events

    else:
        correction = evaluator.evaluate("nominal", eta, r9)
        # When creating the JSON I already included added the variation to the returned value
        uncertainty_up = evaluator.evaluate("up", eta, r9) / correction
        uncertainty_dn = evaluator.evaluate("down", eta, r9) / correction
        # coffea does the unflattenning step itself and sets this value as pt of the up/down variations
        return (
            np.concatenate(
                (uncertainty_up[:, None], uncertainty_dn[:, None]), axis=1
            )
            * _pt[:, None]
        )


def PhotonIDMVAShape(mvaID, events, year="2024", is_correction=True, workflow="base"):
    """
    Photon MVA ID shape systematic derived via quantile mapping of the
    data--MC mvaID discrepancy.

    There is no correction associated with this systematic -- the normalizing
    flows *are* the correction.  Calling with ``is_correction=True`` raises
    a ``RuntimeError``.

    When ``is_correction=False`` the function loads the pre-derived
    correctionlib JSON, evaluates an additive shift delta(mvaID) for every
    photon, and returns the up / down varied mvaID values as an (N, 2) array.

    https://indico.cern.ch/event/1624987/#31-shape-systematic-on-the-pho
    """
    if is_correction:
        raise RuntimeError(
            "PhotonIDMVAShape has no correction to apply -- the normalizing "
            "flow corrections fulfil that role.  Use --doFlow-corrections "
            "to apply them, and request this entry only as a systematic."
        )

    if not hasattr(events, "GenPart"):
        raise ValueError("PhotonIDMVAShape corrections should only be applied to MC!")

    if len(mvaID) > 0 and ak.all(mvaID == ak.flatten(events.Photon.mvaID)):
        raise ValueError("mvaID values are identical to those in the original NanoAOD. You must apply flow corrections to use `PhotonIDMVAShape`!")

    # For Lowmass: a linear shift of the mvaID values, with different parameters for
    # EB and EE is applied as a preliminary estimate of the systematic uncertainty
    if workflow == "lowmass":
        _mvaID = ak.flatten(events.Photon.mvaID)[:, None]
        _isScEtaEB = ak.flatten(events.Photon.isScEtaEB)[:, None]

        x_min = {"EB": 0, "EE": 0}
        x_max = {"EB": 1, "EE": 1}
        y_min = {"EB": 0.04, "EE": 0.05}
        y_max = {"EB": 0.01, "EE": 0.01}
        slope = {
            "EB": (y_min["EB"] - y_max["EB"]) / (-x_min["EB"] + x_max["EB"]),
            "EE": (y_min["EE"] - y_max["EE"]) / (-x_min["EE"] + x_max["EE"]),
        }

        Up_IDMVA = {}
        Down_IDMVA = {}

        for region in ["EB", "EE"]:
            Up_IDMVA[region] = np.where(
                _mvaID < x_min[region],
                (_mvaID + y_min[region]),
                (y_min[region] - ((_mvaID - x_min[region]) * slope[region]) + _mvaID),
            )
            Up_IDMVA[region] = np.clip(Up_IDMVA[region], None, 1)
            Down_IDMVA[region] = np.where(
                _mvaID < x_min[region],
                (_mvaID - y_min[region]),
                (-y_min[region] + ((_mvaID - x_min[region]) * slope[region]) + _mvaID),
            )
            Down_IDMVA[region] = np.clip(Down_IDMVA[region], -1, None)

        uncertainty_up = ak.where(_isScEtaEB, Up_IDMVA["EB"], Up_IDMVA["EE"])
        uncertainty_down = ak.where(_isScEtaEB, Down_IDMVA["EB"], Down_IDMVA["EE"])

        return np.concatenate((uncertainty_up, uncertainty_down), axis=1)

    jsonpog_file = os.path.join(
        os.path.dirname(__file__),
        f"JSONs/PhotonIDMVAShape/{year}/PhotonIDMVAShape.json.gz",
    )
    evaluator = correctionlib.CorrectionSet.from_file(jsonpog_file)["PhotonIDMVAShape"]

    delta = evaluator.evaluate(mvaID)

    mvaID_up = np.clip(mvaID + delta, -1.0, 1.0)
    mvaID_down = np.clip(mvaID - delta, -1.0, 1.0)

    return ak.concatenate(
        [mvaID_up[:, None], mvaID_down[:, None]],
        axis=1,
    )
