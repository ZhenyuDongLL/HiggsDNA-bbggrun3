import awkward as ak
from coffea.nanoevents.methods.nanoaod import behavior as nano_behavior
from coffea.nanoevents.methods.base import Systematic
from coffea.nanoevents.methods.systematics.UpDownSystematic import UpDownSystematic
from copy import copy


@ak.behaviors.mixins.mixin_class(nano_behavior)
class UpDownMultiSystematic(UpDownSystematic):
    """
    Up/Down systematic that supports multi-field `what` by zipping fields and
    passing a single record to the varying function, then writing the result
    back field-by-field.

    Use when `what` is a sequence (e.g. ["pt", "phi"]) such as for MET.

    Contract with the varying function:
      - Coffea will call:   varying_function(whatarg, *args, **kwargs)
      - If `what` is a list/tuple, `whatarg` is a *flattened* record with those fields,
        e.g. `whatarg.pt` and `whatarg.phi`, each shape (M,).
      - The varying function must return a record-of-arrays with the *same* fields,
        each of shape (M, 2) where [:,0] = Up and [:,1] = Down.

    For scalar `what` (str) and for "weight", behavior matches UpDownSystematic.
    """
    def _build_variations(self, name, what, varying_function, *args, **kwargs):
        if what == "weight":
            whatarg = self["__systematics__", "__ones__"]
        elif isinstance(what, (list, tuple)):
            whatarg = ak.zip({f: self[f] for f in what}, depth_limit=1)
        else:
            whatarg = self[what]

        self["__systematics__", f"__{name}__"] = varying_function(
            whatarg, *args, **kwargs
        )

    def get_variation(self, name, what, astype, updown):
        # Select up(0)/down(1) slice out of the variation payload
        varied = self["__systematics__", f"__{name}__", :, self._udmap[updown]]

        # Write results back into the base record, field-by-field if needed
        fields = ak.fields(self)
        if "__systematics__" in fields:
            fields.remove("__systematics__")

        out = {field: self[field] for field in fields}

        if what == "weight":
            out[f"weight_{name}"] = varied
        elif isinstance(what, (list, tuple)):
            for f in what:
                out[f] = varied[f]
        else:
            out[what] = varied

        params = copy(self.layout.parameters)
        params["variation"] = f"{name}-{what}-{updown}"

        return ak.zip(out, depth_limit=1, parameters=params,
                      behavior=self.behavior, with_name=astype)


# Let Coffea know this kind exists
nano_behavior[("__typestr__", "UpDownMultiSystematic")] = "UpDownMultiSystematic"
Systematic.add_kind("UpDownMultiSystematic")
