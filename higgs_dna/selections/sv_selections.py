from higgs_dna.selections.object_selections import delta_r_mask
import awkward


def match_sv(
    self,
    jets: awkward.highlevel.Array,
    sv: awkward.highlevel.Array,
    lead_only: bool
) -> awkward.highlevel.Array:
    if lead_only:
        jets = awkward.firsts(jets)
    dr_max = 0.4
    dr_dipho_cut = delta_r_mask(sv, jets, dr_max)

    return (
        ~dr_dipho_cut
    )
