import numpy as np
import awkward as ak
import correctionlib
import os
import logging

logger = logging.getLogger(__name__)


def Tau_EnergyScale(ptmass, *, events, year, is_correction=True):
    year_mapping = {
        "2022preEE": "2022_Summer22",
        "2022postEE": "2022_Summer22EE",
        "2023preBPix": "2023_Summer23",
        "2023postBPix": "2023_Summer23BPix",
        "2024": "2024_Summer24",
    }

    folder_name = year_mapping.get(year)
    json_name = "tau_energy_scale"
    if folder_name and json_name:
        path_to_json = os.path.join(os.path.dirname(__file__), f"JSONs/POG/TAU/{folder_name}/tau.json.gz")
        evaluator = correctionlib.CorrectionSet.from_file(path_to_json)[json_name]
    else:
        raise ValueError(f"Year {year} not supported for TAU ES.")

    name = "DeepTau2018v2p5"
    wp_VSjet = "Medium"
    wp_VSe = "Tight"

    taus = events.Tau
    counts = ak.num(taus)
    taus = ak.flatten(taus)

    # mask dm 5 and 6 taus to ensure correctionlib does not crash (such taus are dropped in the tau preselection)
    mask = (taus.decayMode == 5) | (taus.decayMode == 6)
    masked_dm = ak.where(mask, 0, taus.decayMode)
    masked_dm = ak.where(taus.decayMode == 2, 1, masked_dm)  # map 2->1 for TES json (no separate bin for 2pi0)

    if is_correction:
        logger.info(f"Applying TAU energy scale correction for year {year}")
        correction = evaluator.evaluate(taus.pt, taus.eta, masked_dm, taus.genPartFlav, name, wp_VSjet, wp_VSe, "nom")

        taus["pt"] = taus.pt * correction
        taus["mass"] = ak.where((taus.decayMode != 0), taus.mass * correction, taus.mass)
        events["Tau"] = ak.unflatten(taus, counts)
        return events

    else:
        logger.info(f"Applying TAU energy scale systematic variation for year {year}")
        nominal_correction = evaluator.evaluate(taus.pt, taus.eta, masked_dm, taus.genPartFlav, name, wp_VSjet, wp_VSe, "nom")
        uncertainty_up = evaluator.evaluate(taus.pt, taus.eta, masked_dm, taus.genPartFlav, name, wp_VSjet, wp_VSe, "up")
        uncertainty_down = evaluator.evaluate(taus.pt, taus.eta, masked_dm, taus.genPartFlav, name, wp_VSjet, wp_VSe, "down")
        up_variation = uncertainty_up / nominal_correction
        down_variation = uncertainty_down / nominal_correction

        pt_var = ak.concatenate((up_variation[:, None], down_variation[:, None]), axis=1) * taus.pt[:, None]
        varied_mass = ak.concatenate((up_variation[:, None], down_variation[:, None]), axis=1) * taus.mass[:, None]
        unvaried_mass = ak.concatenate((np.ones(len(taus))[:, None], np.ones(len(taus))[:, None]), axis=1) * taus.mass[:, None]
        mass_var = ak.where((taus.decayMode != 0)[:, None], varied_mass, unvaried_mass)
        return ak.zip({"pt": pt_var, "mass": mass_var}, depth_limit=1)
