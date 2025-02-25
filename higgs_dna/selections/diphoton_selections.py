import awkward
import numpy
import warnings


def apply_fiducial_cut_det_level(
    self,
    diphotons: awkward.Array,
) -> awkward.Array:
    lead_rel_iso = diphotons.pho_lead.pfRelIso03_all if hasattr(diphotons.pho_lead, "pfRelIso03_all") else diphotons.pho_lead.pfRelIso03_all_quadratic  # photons.pfRelIso03_chg for v11, photons.pfRelIso03_chg_quadratic v12 and above
    sublead_rel_iso = diphotons.pho_sublead.pfRelIso03_all if hasattr(diphotons.pho_sublead, "pfRelIso03_all") else diphotons.pho_sublead.pfRelIso03_all_quadratic  # photons.pfRelIso03_chg for v11, photons.pfRelIso03_chg_quadratic v12 and above
    # Determine if event passes fiducial Hgg cuts at detector-level
    if self.fiducialCuts == 'classical':
        fid_det_passed = (diphotons.pho_lead.pt / diphotons.mass > 1 / 3) & (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & (lead_rel_iso * diphotons.pho_lead.pt < 10) & ((sublead_rel_iso * diphotons.pho_sublead.pt) < 10) & (numpy.abs(diphotons.pho_lead.eta) < 2.5) & (numpy.abs(diphotons.pho_sublead.eta) < 2.5)

    elif self.fiducialCuts == 'geometric':
        fid_det_passed = (numpy.sqrt(diphotons.pho_lead.pt * diphotons.pho_sublead.pt) / diphotons.mass > 1 / 3) & (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & (lead_rel_iso * diphotons.pho_lead.pt < 10) & (sublead_rel_iso * diphotons.pho_sublead.pt < 10) & (numpy.abs(diphotons.pho_lead.eta) < 2.5) & (numpy.abs(diphotons.pho_sublead.eta) < 2.5)

    elif self.fiducialCuts == 'store_flag':
        fid_classical = (diphotons.pho_lead.pt / diphotons.mass > 1 / 3) & (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & (lead_rel_iso * diphotons.pho_lead.pt < 10) & ((sublead_rel_iso * diphotons.pho_sublead.pt) < 10) & (numpy.abs(diphotons.pho_lead.eta) < 2.5) & (numpy.abs(diphotons.pho_sublead.eta) < 2.5)
        fid_geometric = (numpy.sqrt(diphotons.pho_lead.pt * diphotons.pho_sublead.pt) / diphotons.mass > 1 / 3) & (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & (lead_rel_iso * diphotons.pho_lead.pt < 10) & (sublead_rel_iso * diphotons.pho_sublead.pt < 10) & (numpy.abs(diphotons.pho_lead.eta) < 2.5) & (numpy.abs(diphotons.pho_sublead.eta) < 2.5)
        diphotons['pass_fiducial_classical'] = fid_classical
        diphotons['pass_fiducial_geometric'] = fid_geometric

    elif self.fiducialCuts == 'classical_noIso':
        fid_det_passed = (diphotons.pho_lead.pt / diphotons.mass > 1 / 3) & (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & (diphotons.mass > 100) & (diphotons.mass < 180)

    elif self.fiducialCuts == 'none':
        fid_det_passed = diphotons.pho_lead.pt > -10  # This is a very dummy way but I do not know how to make a true array of outer shape of diphotons

    else:
        warnings.warn("You chose %s the fiducialCuts mode, but this is currently not supported. You should check your settings. For this run, no fiducial selection at detector level is applied." % self.fiducialCuts)
        fid_det_passed = diphotons.pho_lead.pt > -1

    return diphotons[fid_det_passed] if self.fiducialCuts != 'store_flag' else diphotons
