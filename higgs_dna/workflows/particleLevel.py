from higgs_dna.workflows.skeleton import HggSkeletonProcessor
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_corrections as available_weight_corrections
from higgs_dna.utils.dumping_utils import diphoton_ak_array, dump_ak_array, diphoton_list_to_pandas, dump_pandas

from higgs_dna.tools.gen_helpers import (
    get_fiducial_flag,
    get_genJets,
    get_higgs_gen_attributes,
    get_higgs_truth_attributes,
)
from higgs_dna.utils.misc_utils import choose_jet, DPhiV1V2, rapidity_from_pt_eta_mass

from typing import Any, Dict, List, Optional
import awkward as ak
import logging
import warnings
import numpy
import sys
from coffea.analysis_tools import Weights

logger = logging.getLogger(__name__)


class ParticleLevelProcessor(HggSkeletonProcessor):
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        systematics: Dict[str, List[Any]] = None,
        corrections: Dict[str, List[Any]] = None,
        apply_trigger: bool = False,
        nano_version: int = None,
        bTagEffFileName: Optional[str] = None,
        output_location: Optional[str] = None,
        taggers: Optional[List[Any]] = None,
        trigger_group=".*DoubleEG.*",
        analysis="mainAnalysis",
        applyCQR: bool = False,
        skipJetVetoMap: bool = False,
        year: Dict[str, List[str]] = None,
        fiducialCuts: str = "none",
        doDeco: bool = False,
        Smear_sigma_m: bool = False,
        doFlow_corrections: bool = False,
        validate_with_electrons: bool = False,
        output_format: str = "parquet",
    ) -> None:
        if validate_with_electrons:
            raise ValueError(f"Validation with electrons is not supported in {self.__class__.__name__}")
        super().__init__(
            metaconditions,
            systematics=systematics,
            corrections=corrections,
            apply_trigger=apply_trigger,
            nano_version=nano_version,
            bTagEffFileName=bTagEffFileName,
            output_location=output_location,
            taggers=taggers,
            trigger_group=trigger_group,
            analysis=analysis,
            applyCQR=applyCQR,
            skipJetVetoMap=skipJetVetoMap,
            year=year,
            fiducialCuts=fiducialCuts,
            doDeco=doDeco,
            Smear_sigma_m=Smear_sigma_m,
            doFlow_corrections=doFlow_corrections,
            validate_with_electrons=validate_with_electrons,
            output_format=output_format
        )

    def process_extra(self, events: ak.Array) -> ak.Array:
        return events, {}

    def process(self, events: ak.Array) -> Dict[Any, Any]:
        dataset_name = events.metadata["dataset"]

        # metadata array to append to higgsdna output
        metadata = {}

        # data or monte carlo?
        self.data_kind = "mc" if hasattr(events, "GenPart") else "data"

        if self.data_kind == "data":
            logger.info("The 'particleLevel' processor can only be run on MC. Aborting now...")
            sys.exit(0)

        # here we start recording possible coffea accumulators
        # most likely histograms, could be counters, arrays, ...
        histos_etc = {}
        histos_etc[dataset_name] = {}
        if self.data_kind == "mc":
            histos_etc[dataset_name]["nTot"] = int(
                ak.num(events.genWeight, axis=0)
            )
            histos_etc[dataset_name]["nPos"] = int(ak.sum(events.genWeight > 0))
            histos_etc[dataset_name]["nNeg"] = int(ak.sum(events.genWeight < 0))
            histos_etc[dataset_name]["nEff"] = int(
                histos_etc[dataset_name]["nPos"] - histos_etc[dataset_name]["nNeg"]
            )
            histos_etc[dataset_name]["genWeightSum"] = float(
                numpy.sum(events.genWeight.to_numpy())
            )
        else:
            histos_etc[dataset_name]["nTot"] = int(len(events))
            histos_etc[dataset_name]["nPos"] = int(histos_etc[dataset_name]["nTot"])
            histos_etc[dataset_name]["nNeg"] = int(0)
            histos_etc[dataset_name]["nEff"] = int(histos_etc[dataset_name]["nTot"])
            histos_etc[dataset_name]["genWeightSum"] = float(len(events))

        # Add sum of gen weights before selection for normalisation in postprocessing
        metadata["sum_genw_presel"] = str(numpy.sum(events.genWeight.to_numpy()))

        # read which systematics and corrections to process
        try:
            correction_names = self.corrections[dataset_name]
        except KeyError:
            correction_names = []

        for correction_name in correction_names:
            if correction_name in available_object_corrections.keys():
                logger.info(
                    f"Applying correction {correction_name} to dataset {dataset_name}"
                )
                varying_function = available_object_corrections[correction_name]
                events = varying_function(events=events, year=self.year[dataset_name][0])
            elif correction_name in available_weight_corrections:
                # event weight corrections will be applied after photon preselection / application of further taggers
                continue
            else:
                # may want to throw an error instead, needs to be discussed
                warnings.warn(f"Could not process correction {correction_name}.")
                continue

        # Filling with some dummy values
        diphotons = ak.Array({"pt": numpy.ones(len(events))})

        diphotons['fiducialClassicalFlag'] = get_fiducial_flag(events, flavour='Classical')
        diphotons['fiducialGeometricFlag'] = get_fiducial_flag(events, flavour='Geometric')
        GenPTH, GenYH, GenPhiH, GenLeadPho, GenSubleadPho = get_higgs_gen_attributes(events)
        TruthPTH, TruthYH = get_higgs_truth_attributes(events)

        genJets = get_genJets(
            events,
            pt_cut=30.,
            eta_cut=4.7,
            jet_pho_min_dr=self.jet_pho_min_dr,
            jet_ele_min_dr=self.jet_ele_min_dr,
            jet_muo_min_dr=self.jet_muo_min_dr,
            electron_pt_threshold=self.electron_pt_threshold,
            electron_max_eta=self.electron_max_eta,
            muon_pt_threshold=self.muon_pt_threshold,
            muon_max_eta=self.muon_max_eta,
        )
        genJets_absEta2p5 = genJets[numpy.abs(genJets.eta) < 2.5]

        ######################
        # Diphoton Variables #
        ######################
        GenPTH = ak.fill_none(GenPTH, -999.0)
        diphotons['GenPTH'] = GenPTH

        GenYH = ak.fill_none(GenYH, -999)
        GenYH = ak.where(numpy.isnan(GenYH), -999, GenYH)
        diphotons['GenYH'] = GenYH

        diphotons["GenLeadPt"] = ak.fill_none(GenLeadPho.pt, -999.0)
        diphotons["GenLeadEta"] = ak.fill_none(GenLeadPho.eta, -999.0)
        diphotons["GenLeadPhi"] = ak.fill_none(GenLeadPho.phi, -999.0)
        diphotons["GenSubleadPt"] = ak.fill_none(GenSubleadPho.pt, -999.0)
        diphotons["GenSubleadEta"] = ak.fill_none(GenSubleadPho.eta, -999.0)
        diphotons["GenSubleadPhi"] = ak.fill_none(GenSubleadPho.phi, -999.0)

        GenDeltaPhoPhi = GenLeadPho.phi - GenSubleadPho.phi
        GenDeltaPhoPhi_pi_array = ak.full_like(GenDeltaPhoPhi, 2 * numpy.pi)
        # Select the smallest angle
        GenDeltaPhoPhi = ak.where(
            GenDeltaPhoPhi > numpy.pi,
            GenDeltaPhoPhi - GenDeltaPhoPhi_pi_array,
            GenDeltaPhoPhi
        )
        GenDeltaPhoPhi = ak.where(
            GenDeltaPhoPhi < -numpy.pi,
            GenDeltaPhoPhi + GenDeltaPhoPhi_pi_array,
            GenDeltaPhoPhi
        )
        GenAcop = ak.full_like(GenDeltaPhoPhi, numpy.pi) - GenDeltaPhoPhi
        GenThetaEtaStar = numpy.tan(GenAcop / 2) / numpy.cosh((GenLeadPho.eta - GenSubleadPho.eta) / 2)
        GenThetaEtaStar = ak.fill_none(GenThetaEtaStar, -999.0)
        diphotons['GenThetaEtaStar'] = GenThetaEtaStar

        GenAbsDeltaPhoPhi = numpy.abs(GenDeltaPhoPhi)
        GenPhiAcop = ak.full_like(GenAbsDeltaPhoPhi, numpy.pi) - GenAbsDeltaPhoPhi
        GenPhiEtaStar = numpy.tan(GenPhiAcop / 2) / numpy.cosh((GenLeadPho.eta - GenSubleadPho.eta) / 2)
        GenAbsPhiEtaStar = numpy.abs(GenPhiEtaStar)
        GenAbsPhiEtaStar = ak.fill_none(GenAbsPhiEtaStar, -999.0)
        diphotons["GenAbsPhiEtaStar"] = GenAbsPhiEtaStar

        GenDiphoton = GenLeadPho + GenSubleadPho
        GenDiPhoMass = GenDiphoton.mass
        GenDiPhoPT = GenDiphoton.pt
        GenCosThetaStarCS = 2 * (((GenLeadPho.pz * GenSubleadPho.energy) - (GenLeadPho.energy * GenSubleadPho.pz)) / (GenDiPhoMass * numpy.sqrt(GenDiPhoMass**2 + GenDiPhoPT**2)))
        GenCosThetaStarCS = ak.fill_none(GenCosThetaStarCS, -999.0)
        diphotons['GenCosThetaStarCS'] = GenCosThetaStarCS

        #########################
        # Leading Jet Variables #
        #########################
        # Choose zero (leading) jet and pad with -999 if none
        GenPTJ0 = choose_jet(genJets.pt, 0, -999.0)
        diphotons['GenPTJ0'] = GenPTJ0
        diphotons["GenPTJ0_pt30_absEta2p5"] = choose_jet(genJets_absEta2p5.pt, 0, -999.0)

        gen_first_jet_eta = choose_jet(genJets.eta, 0, -999.0)
        gen_first_jet_mass = choose_jet(genJets.mass, 0, -999.0)
        gen_first_jet_phi = choose_jet(genJets.phi, 0, -999.0)

        diphotons['gen_first_jet_eta'] = gen_first_jet_eta
        diphotons['gen_first_jet_mass'] = gen_first_jet_mass
        diphotons['gen_first_jet_phi'] = gen_first_jet_phi

        with numpy.errstate(over='ignore', invalid='ignore'):
            gen_first_jet_pz = GenPTJ0 * numpy.sinh(gen_first_jet_eta)
            gen_first_jet_pz = ak.where(gen_first_jet_eta == -999, -999, gen_first_jet_pz)
            gen_first_jet_energy = numpy.sqrt((GenPTJ0**2 * numpy.cosh(gen_first_jet_eta)**2) + gen_first_jet_mass**2)

            GenYJ0 = 0.5 * numpy.log((gen_first_jet_energy + gen_first_jet_pz) / (gen_first_jet_energy - gen_first_jet_pz))

        GenYJ0 = ak.fill_none(GenYJ0, -999)
        GenYJ0 = ak.where(numpy.isnan(GenYJ0), -999, GenYJ0)
        diphotons['GenYJ0'] = GenYJ0

        GenPTJ0_absEta2p5 = choose_jet(genJets_absEta2p5.pt, 0, -999.0)
        gen_first_jet_eta_absEta2p5 = choose_jet(genJets_absEta2p5.eta, 0, -999.0)
        gen_first_jet_mass_absEta2p5 = choose_jet(genJets_absEta2p5.mass, 0, -999.0)
        diphotons["GenYJ0_pt30_absEta2p5"] = rapidity_from_pt_eta_mass(
            GenPTJ0_absEta2p5,
            gen_first_jet_eta_absEta2p5,
            gen_first_jet_mass_absEta2p5,
            fill_value=-999.0,
        )

        GenDYHJ0 = GenYJ0 - GenYH
        # Set all entries above 500 to -999
        GenDYHJ0 = ak.where(
            numpy.abs(GenDYHJ0) > 500,
            -999,
            GenDYHJ0
        )
        GenDYHJ0 = ak.fill_none(GenDYHJ0, -999.0)
        diphotons["GenDYHJ0"] = GenDYHJ0

        GenHPhi = ak.fill_none(GenPhiH, -999)

        GenDPhiHJ0 = gen_first_jet_phi - GenHPhi
        GenDPhiHJ0 = (GenDPhiHJ0 + numpy.pi) % (2 * numpy.pi) - numpy.pi
        GenDPhiHJ0 = ak.where(GenHPhi == -999, -999, GenDPhiHJ0)
        GenDPhiHJ0 = ak.where(gen_first_jet_phi == -999, -999, GenDPhiHJ0)
        diphotons["GenDPhiHJ0"] = GenDPhiHJ0

        #################################
        # Next-to-leading Jet Variables #
        #################################
        GenPTJ1 = choose_jet(genJets.pt, 1, -999.0)
        diphotons['GenPTJ1'] = GenPTJ1

        gen_second_jet_eta = choose_jet(genJets.eta, 1, -999.0)
        gen_second_jet_mass = choose_jet(genJets.mass, 1, -999.0)
        gen_second_jet_phi = choose_jet(genJets.phi, 1, -999.0)

        diphotons['gen_second_jet_eta'] = gen_second_jet_eta
        diphotons['gen_second_jet_mass'] = gen_second_jet_mass
        diphotons['gen_second_jet_phi'] = gen_second_jet_phi

        with numpy.errstate(over='ignore', invalid='ignore'):
            gen_second_jet_pz = GenPTJ1 * numpy.sinh(gen_second_jet_eta)
            gen_second_jet_pz = ak.where(gen_second_jet_eta == -999, -999, gen_second_jet_pz)
            gen_second_jet_energy = numpy.sqrt((GenPTJ1**2 * numpy.cosh(gen_second_jet_eta)**2) + gen_second_jet_mass**2)

            GenYJ1 = 0.5 * numpy.log((gen_second_jet_energy + gen_second_jet_pz) / (gen_second_jet_energy - gen_second_jet_pz))

        GenYJ1 = ak.fill_none(GenYJ1, -999)
        GenYJ1 = ak.where(numpy.isnan(GenYJ1), -999, GenYJ1)
        diphotons['GenYJ1'] = GenYJ1

        GenDYJ0J1 = GenYJ0 - GenYJ1
        # Set all entries above 500 to -999
        GenDYJ0J1 = ak.where(
            numpy.abs(GenDYJ0J1) > 500,
            -999,
            GenDYJ0J1
        )
        # Set all entries which are precisely 0 to -999
        GenDYJ0J1 = ak.where(
            GenDYJ0J1 == 0,
            -999,
            GenDYJ0J1
        )
        GenDYJ0J1 = ak.fill_none(GenDYJ0J1, -999.0)
        diphotons["GenDYJ0J1"] = GenDYJ0J1

        gen_first_jet_vector = ak.zip({
            "pt": GenPTJ0,
            "eta": gen_first_jet_eta,
            "phi": gen_first_jet_phi,
            "mass": gen_first_jet_mass
        }, with_name="Momentum4D")

        gen_second_jet_vector = ak.zip({
            "pt": GenPTJ1,
            "eta": gen_second_jet_eta,
            "phi": gen_second_jet_phi,
            "mass": gen_second_jet_mass
        }, with_name="Momentum4D")

        GenDPhiJ0J1 = DPhiV1V2(gen_first_jet_vector, gen_second_jet_vector)
        diphotons["GenDPhiJ0J1"] = GenDPhiJ0J1

        padded_genJets = genJets[ak.argsort(genJets.pt, ascending=False)]
        # First build the dijet system out of the leading and subleading jet (in pt)
        padded_genJets = ak.pad_none(genJets, 2)
        genDijet = padded_genJets[:, 0] + padded_genJets[:, 1]

        GenMassJ0J1 = ak.fill_none(genDijet.mass, -999.0)
        diphotons["GenMassJ0J1"] = GenMassJ0J1

        GenDijetEta = ak.fill_none(genDijet.eta, -999.0)
        GenDiphotonEta = ak.fill_none(GenDiphoton.eta, -999.0)
        GenDEtaJ0J1H = GenDijetEta - GenDiphotonEta
        # Set all entries which are precisely 0 to -999
        GenDEtaJ0J1H = ak.where(
            GenDEtaJ0J1H == 0,
            -999,
            GenDEtaJ0J1H
        )
        # Set all entries which are above 500 in absolute value to -999 (come from either no diphoton or no dijet system)
        GenDEtaJ0J1H = ak.where(
            numpy.abs(GenDEtaJ0J1H) > 500,
            -999,
            GenDEtaJ0J1H
        )
        GenDEtaJ0J1H = ak.fill_none(GenDEtaJ0J1H, -999.0)
        diphotons["GenDEtaJ0J1H"] = GenDEtaJ0J1H

        GenDijetPhi = ak.fill_none(genDijet.phi, -999)
        GenHPhi = ak.fill_none(GenPhiH, -999)

        GenDPhiHJ0J1 = GenDijetPhi - GenPhiH
        GenDPhiHJ0J1 = (GenDPhiHJ0J1 + numpy.pi) % (2 * numpy.pi) - numpy.pi
        GenDPhiHJ0J1 = ak.where(GenHPhi == -999, -999, GenDPhiHJ0J1)
        GenDPhiHJ0J1 = ak.where(GenDijetPhi == -999, -999, GenDPhiHJ0J1)
        diphotons["GenDPhiHJ0J1"] = GenDPhiHJ0J1

        GenEtaJ0J1 = gen_first_jet_eta - gen_second_jet_eta
        # Set all entries which are precisely 0 to -999
        GenEtaJ0J1 = ak.where(
            GenEtaJ0J1 == 0,
            -999,
            GenEtaJ0J1
        )
        # Set all entries which are above 500 in absolute value to -999 (come from either no diphoton or no dijet system)
        GenEtaJ0J1 = ak.where(
            numpy.abs(GenEtaJ0J1) > 500,
            -999,
            GenEtaJ0J1
        )
        GenEtaJ0J1 = ak.fill_none(GenEtaJ0J1, -999.0)
        diphotons["GenEtaJ0J1"] = GenEtaJ0J1

        ###########################
        # Event Level Observables #
        ###########################
        diphotons['GenNJ'] = ak.num(genJets)
        diphotons["GenNJ_pt30_absEta2p5"] = ak.num(genJets_absEta2p5)

        # B-Jets
        # Following the recommendations of https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools for hadronFlavour
        # and the Run 2 recommendations for the bjets
        genJetCondition = (genJets.pt > 30) & (numpy.abs(genJets.eta) < 4.7)
        genBJetCondition = genJetCondition & (genJets.hadronFlavour == 5)
        genJets = ak.with_field(genJets, genBJetCondition, "GenIsBJet")
        num_bjets = ak.sum(genJets["GenIsBJet"], axis=-1)
        diphotons["GenNBJet"] = num_bjets

        gen_first_bjet_pt = choose_jet(genJets[genJets["GenIsBJet"] == True].pt, 0, -999.0)
        diphotons["GenPTbJ0"] = gen_first_bjet_pt

        # Jet Rapidity Observable
        # Iterate over max six largest pt jets to compute tauJC
        GenTauJC_list = []
        GenTauJC_maxJets = 10
        for i in range(GenTauJC_maxJets):
            mass = choose_jet(genJets.mass, i, -999.0)
            pt = choose_jet(genJets.pt, i, -999.0)
            eta = choose_jet(genJets.eta, i, -999.0)

            with numpy.errstate(over='ignore', invalid='ignore'):
                cosh_eta = numpy.cosh(eta)
                sinh_eta = numpy.sinh(eta)

                energy = numpy.sqrt((pt**2 * cosh_eta**2) + mass**2)
                pz = pt * sinh_eta

                # If energy or pz is inf (cause of the hyperbolic functions), set them to -999
                # later set every GenTauJC to -999 which has a value of precisely 0 (corresponding to a transverse momentum of exactly 0 GeV)
                energy = ak.where(numpy.isinf(energy), -999, energy)
                pz = ak.where(numpy.isinf(pz), -999, pz)

                transverse_mass = numpy.sqrt(pt**2 + mass**2)

                y = numpy.log((numpy.sqrt((mass**2 + pt**2) * cosh_eta**2) + (pt * sinh_eta)) / transverse_mass)

                tau_jc = numpy.sqrt(energy**2 - pz**2) / (2 * numpy.cosh(y - GenYH))
            GenTauJC_list.append(tau_jc)

            logger.debug(f"GenTauJC: Jet {i}: Energy={energy}, Pz={pz}, TauJC={tau_jc}")

        # Convert to awkward array for proper axis manipulation
        GenTauJC_array = ak.Array(GenTauJC_list)
        GenTauJC = ak.max(GenTauJC_array, axis=0)
        GenTauJC = ak.where(GenTauJC == 0, -999, GenTauJC)
        GenTauJC = ak.fill_none(GenTauJC, -999.0)
        diphotons["GenTauJC"] = GenTauJC

        # workflow specific processing
        events, process_extra = self.process_extra(events)
        histos_etc.update(process_extra)

        # set diphotons as part of the event record
        events["diphotons"] = diphotons
        # annotate diphotons with event information
        diphotons["event"] = events.event
        diphotons["lumi"] = events.luminosityBlock
        diphotons["run"] = events.run
        # annotate diphotons with dZ information (difference between z position of GenVtx and PV) as required by flashggfinalfits
        diphotons["genWeight"] = events.genWeight
        diphotons["dZ"] = events.GenVtx.z - events.PV.z
        diphotons["TruthPTH"] = TruthPTH
        diphotons["TruthYH"] = TruthYH

        for i in range(9):
            diphotons["LHEScaleWeight_" + str(i)] = events.LHEScaleWeight[:,i]
        for i in range(103):
            diphotons["LHEPdfWeight_" + str(i)] = events.LHEPdfWeight[:,i]

        # return if there is no surviving events
        if len(diphotons) == 0:
            logger.info("No surviving events in this run!")

        # Retain all events
        selection_mask = numpy.ones(len(diphotons), dtype=bool)
        # initiate Weight container here, after selection, since event selection cannot easily be applied to weight container afterwards
        event_weights = Weights(size=len(events[selection_mask]))
        # set weights to generator weights
        event_weights._weight = ak.to_numpy(events["genWeight"][selection_mask])
        # corrections to event weights:
        for correction_name in correction_names:
            if correction_name in available_weight_corrections:
                logger.info(
                    f"Adding correction {correction_name} to weight collection of dataset {dataset_name}"
                )
                varying_function = available_weight_corrections[
                    correction_name
                ]
                event_weights = varying_function(
                    events=events[selection_mask],
                    photons=events["diphotons"][
                        selection_mask
                    ],
                    weights=event_weights,
                    dataset_name=dataset_name,
                    year=self.year[dataset_name][0],
                )
        diphotons["weight_central"] = event_weights.weight() / events["genWeight"][selection_mask]  # Here, if diphotons none, then also the weight is None.
        diphotons["weight"] = event_weights.weight()

        if self.output_location is not None:
            if self.output_format == "root":
                df = diphoton_list_to_pandas(self, diphotons)
            else:
                akarr = diphoton_ak_array(self, diphotons)
            fname = (
                events.attrs[
                    "@events_factory"
                ]._partition_key.replace("/", "_")
                + ".%s" % self.output_format
            )

            subdirs = []
            if "dataset" in events.metadata:
                subdirs.append(events.metadata["dataset"])
            if self.output_format == "root":
                dump_pandas(self, df, fname, self.output_location, subdirs)
            else:
                dump_ak_array(
                    self, akarr, fname, self.output_location, metadata, subdirs,
                )

        return histos_etc

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass
