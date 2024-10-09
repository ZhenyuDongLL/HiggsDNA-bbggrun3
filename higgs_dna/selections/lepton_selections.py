import awkward
from higgs_dna.selections.object_selections import delta_r_mask


def select_electrons(
    self,
    electrons: awkward.highlevel.Array,
    diphotons: awkward.highlevel.Array,
) -> awkward.highlevel.Array:
    pt_cut = electrons.pt > self.electron_pt_threshold

    eta_cut = abs(electrons.eta) < self.electron_max_eta

    if self.el_id_wp == "WP90":
        id_cut = electrons.mvaIso_WP90
    elif self.el_id_wp == "WP80":
        id_cut = electrons.mvaIso_WP80
    elif self.el_id_wp == "loose":
        id_cut = electrons.cutBased >= 2
    # WPL is not supported anymore with the Run 3 electron ID, CMSSW 130X v12 nanoAODs only have WP80 and WP90 options
    # elif self.el_id_wp == "WPL":
    #    id_cut = electrons.mvaIso_WPL
    else:
        id_cut = electrons.pt > 0.

    dr_phoLead_cut = delta_r_mask(electrons, diphotons.pho_lead, self.electron_photon_min_dr)
    dr_phoSublead_cut = delta_r_mask(electrons, diphotons.pho_sublead, self.electron_photon_min_dr)

    return pt_cut & eta_cut & id_cut & dr_phoLead_cut & dr_phoSublead_cut


def select_muons(
    self,
    muons: awkward.highlevel.Array,
    diphotons: awkward.highlevel.Array
) -> awkward.highlevel.Array:
    pt_cut = muons.pt > self.muon_pt_threshold

    eta_cut = abs(muons.eta) < self.muon_max_eta

    if self.mu_id_wp == "tight":
        id_cut = muons.tightId
    elif self.mu_id_wp == "medium":
        id_cut = muons.mediumId
    elif self.mu_id_wp == "loose":
        id_cut = muons.looseId
    else:
        id_cut = muons.pt > 0

    # if I understand https://twiki.cern.ch/twiki/bin/view/CMS/MuonRun32022#Medium_pT_15_GeV_to_200_GeV correctly, only loose and tight are PF isos are supported with a SF (so far?)
    # also very loose, very tight and very very tight WPs are available: 1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight
    if self.mu_iso_wp == "tight":
        iso_cut = muons.pfIsoId >= 4
    elif self.mu_iso_wp == "medium":
        iso_cut = muons.pfIsoId >= 3
    elif self.mu_iso_wp == "loose":
        iso_cut = muons.pfIsoId >= 2
    else:
        iso_cut = muons.pt > 0

    if self.global_muon:
        global_cut = muons.isGlobal
    else:
        global_cut = muons.pt > 0

    dr_phoLead_cut = delta_r_mask(muons, diphotons.pho_lead, self.muon_photon_min_dr)
    dr_phoSublead_cut = delta_r_mask(muons, diphotons.pho_sublead, self.muon_photon_min_dr)

    return pt_cut & eta_cut & id_cut & iso_cut & dr_phoLead_cut & dr_phoSublead_cut & global_cut
