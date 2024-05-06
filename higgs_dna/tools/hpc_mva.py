from typing import List, Optional, Tuple

import awkward
import numpy
import xgboost as xgb
import logging

logger = logging.getLogger(__name__)


def calculate_ch_vs_ggh_mva(
    self,
    mva: Tuple[Tuple[Optional[xgb.Booster], Optional[xgb.Booster]], List[str]],
    diphotons: awkward.Array,
    events: awkward.Array,
) -> awkward.Array:
    """
    Calculate cH vs ggH bdt scores for events.
    """

    if mva[0] is None:
        return diphotons, events
    elif len(diphotons) == 0:
        logger.info("no events surviving event selection, adding fake ch vs ggh bdt score")
        diphotons["ch_vs_ggh_bdt_score"] = awkward.zeros_like(diphotons.mass)
        return diphotons, events

    ch_vs_ggh = []
    ch_vs_ggh.append(mva[0][0])
    ch_vs_ggh.append(mva[0][1])

    var_order = mva[1]

    events_bdt = events

    events_bdt["customLeadingPhotonIDMVA"] = diphotons.pho_lead.mvaID
    events_bdt["customSubLeadingPhotonIDMVA"] = diphotons.pho_sublead.mvaID
    events_bdt["leadingPhoton_eta"] = diphotons.pho_lead.eta
    events_bdt["subleadingPhoton_eta"] = diphotons.pho_sublead.eta
    events_bdt["leadingPhoton_relpt"] = diphotons.pho_lead.pt / diphotons.mass
    events_bdt["subleadingPhoton_relpt"] = diphotons.pho_sublead.pt / diphotons.mass
    events_bdt["leadingJet_pt"] = diphotons.first_jet_pt
    events_bdt["leadingJet_eta"] = diphotons.first_jet_eta

    lead_jets = awkward.zip(
        {
            "pt": diphotons.first_jet_pt,
            "eta": diphotons.first_jet_eta,
            "phi": diphotons.first_jet_phi,
            "mass": diphotons.first_jet_mass,
            "charge": diphotons.first_jet_charge
        }
    )
    lead_jets = awkward.with_name(lead_jets, "PtEtaPhiMCandidate")

    lpj_dphi = diphotons.pho_lead.delta_phi(lead_jets)
    spj_dphi = diphotons.pho_sublead.delta_phi(lead_jets)

    events_bdt["DeltaPhi_gamma1_cjet"] = lpj_dphi
    events_bdt["DeltaPhi_gamma2_cjet"] = spj_dphi

    events_bdt["nJets_revised"] = awkward.where(
        diphotons.n_jets > 3,
        awkward.ones_like(diphotons.n_jets) * 3,
        diphotons.n_jets
    )

    for name in var_order:
        events_bdt[name] = awkward.fill_none(events_bdt[name], -999.0)

    bdt_features = []
    for x in var_order:
        if isinstance(x, tuple):
            name_flat = "_".join(x)
            events_bdt[name_flat] = events_bdt[x]
            bdt_features.append(name_flat)
        else:
            bdt_features.append(x)

    events_bdt = awkward.values_astype(events_bdt, numpy.float64)
    features_bdt = awkward.to_numpy(events_bdt[bdt_features])

    features_bdt_matrix = xgb.DMatrix(
        features_bdt.view((float, len(features_bdt.dtype.names)))
    )

    scores = []
    for bdt in ch_vs_ggh:
        scores.append(bdt.predict(features_bdt_matrix))

    for var in bdt_features:
        if "dipho" not in var:
            diphotons[var] = events_bdt[var]

    scores_out = awkward.where(
        events.event % 4 < 2,
        scores[0],
        scores[1]
    )

    diphotons["ch_vs_ggh_bdt_score"] = awkward.ones_like(diphotons.mass)
    diphotons["ch_vs_ggh_bdt_score"] = scores_out

    return diphotons, events


def calculate_ch_vs_cb_mva(
    self,
    mva: Tuple[Tuple[Optional[xgb.Booster], Optional[xgb.Booster]], List[str]],
    diphotons: awkward.Array,
    events: awkward.Array,
) -> awkward.Array:
    """
    Calculate cH vs ggH bdt scores for events.
    """

    if mva[0] is None:
        return diphotons, events
    elif len(diphotons) == 0:
        logger.info("no events surviving event selection, adding fake ch vs cb bdt score")
        diphotons["ch_vs_cb_bdt_score"] = awkward.zeros_like(diphotons.mass)
        return diphotons, events

    ch_vs_cb = []
    ch_vs_cb.append(mva[0][0])
    ch_vs_cb.append(mva[0][1])
    var_order = mva[1]

    events_bdt = events

    events_bdt["customLeadingPhotonIDMVA"] = diphotons.pho_lead.mvaID
    events_bdt["customSubLeadingPhotonIDMVA"] = diphotons.pho_sublead.mvaID
    events_bdt["leadingPhoton_eta"] = diphotons.pho_lead.eta
    events_bdt["subleadingPhoton_eta"] = diphotons.pho_sublead.eta
    events_bdt["leadingPhoton_relpt"] = diphotons.pho_lead.pt / diphotons.mass
    events_bdt["subleadingPhoton_relpt"] = diphotons.pho_sublead.pt / diphotons.mass
    events_bdt["leadingJet_pt"] = diphotons.first_jet_pt
    events_bdt["leadingJet_eta"] = diphotons.first_jet_eta

    lead_jets = awkward.zip(
        {
            "pt": diphotons.first_jet_pt,
            "eta": diphotons.first_jet_eta,
            "phi": diphotons.first_jet_phi,
            "mass": diphotons.first_jet_mass,
            "charge": diphotons.first_jet_charge
        }
    )
    lead_jets = awkward.with_name(lead_jets, "PtEtaPhiMCandidate")

    lpj_dphi = diphotons.pho_lead.delta_phi(lead_jets)
    spj_dphi = diphotons.pho_sublead.delta_phi(lead_jets)

    events_bdt["DeltaPhi_gamma1_cjet"] = lpj_dphi
    events_bdt["DeltaPhi_gamma2_cjet"] = spj_dphi

    events_bdt["nJets_revised"] = awkward.where(
        diphotons.n_jets > 3,
        awkward.ones_like(diphotons.n_jets) * 3,
        diphotons.n_jets
    )

    for name in var_order:
        events_bdt[name] = awkward.fill_none(events_bdt[name], -999.0)

    bdt_features = []
    for x in var_order:
        if isinstance(x, tuple):
            name_flat = "_".join(x)
            events_bdt[name_flat] = events_bdt[x]
            bdt_features.append(name_flat)
        else:
            bdt_features.append(x)

    events_bdt = awkward.values_astype(events_bdt, numpy.float64)
    features_bdt = awkward.to_numpy(events_bdt[bdt_features])

    features_bdt_matrix = xgb.DMatrix(
        features_bdt.view((float, len(features_bdt.dtype.names)))
    )

    scores = []
    for bdt in ch_vs_cb:
        scores.append(bdt.predict(features_bdt_matrix))

    for var in bdt_features:
        if "dipho" not in var:
            diphotons[var] = events_bdt[var]

    scores_out = awkward.where(
        events.event % 4 < 2,
        scores[0],
        scores[1]
    )

    diphotons["ch_vs_cb_bdt_score"] = awkward.ones_like(diphotons.mass)
    diphotons["ch_vs_cb_bdt_score"] = scores_out

    return diphotons, events
