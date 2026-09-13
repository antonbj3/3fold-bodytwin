#!/usr/bin/env python3
"""Build a public claims graph over the 35 physiology cells in src/bodytwin/cells/.

Step 1  writes designs/cell_designs.json: one design entry per cell module (id = module name,
        cluster = subsystem package) plus one aggregate node per subsystem package.
Step 2  calls the framework builder (src/bodytwin/framework/bodytwin_build_anchor_graph.py)
        on that catalog, which writes ANCHOR_GRAPH.json through the guarded writer. The builder
        also emits a companion inventory markdown whose section list is fixed to a different
        cluster vocabulary, so it would render only part of this catalog; it is written to a
        scratch file instead of into this directory.
Step 3  writes FOLD_LEDGER.jsonl: one entry per cell run, with the decisive number, the exact
        command that reproduces it relative to the repository root, and its consumers.

Dependency rule (the only rule; no hand-added edges):
  R1  import edge: cell A depends on cell B if A imports B's module. Census over the 35
      modules finds no intra-package import, so R1 contributes no edge.
  R2  shared-quantity edge: four coupling quantities are detected in the cell source by a
      literal regex; each has one designated producer cell (the cell whose own gates measure
      that quantity). Every other cell whose source matches the regex depends on the producer.
        ATP                  r"\\bATP\\b"                        producer mitochondrial_oxphos
        O2                   r"\\bO2\\b|oxygen"                  producer mitochondrial_oxphos
        membrane potential   r"membrane potential|\\bVm\\b|\\bmV\\b"  producer na_k_atpase
        Ca                   r"\\bCa2\\+|\\bCa\\b|calcium"         producer cardiac_cicr_ode_model
      (cortisol, regex r"cortisol", matches no module in this set and yields no edge.)
  R3  precedence, to keep the graph acyclic: the producers are ordered by supply direction
      mitochondrial_oxphos -> na_k_atpase -> cardiac_cicr_ode_model (energy supply before
      electrical state before Ca release). An edge from producer X to producer Y is kept only
      if Y precedes X in that order; self-edges and duplicates are dropped.
  R4  aggregate edge: each subsystem aggregate node depends on every cell in its package.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CELLS = REPO / "src" / "bodytwin" / "cells"

QUANTITIES = [
    ("ATP", r"\bATP\b", "mitochondrial_oxphos"),
    ("O2", r"\bO2\b|oxygen", "mitochondrial_oxphos"),
    ("membrane_potential", r"membrane potential|\bVm\b|\bmV\b", "na_k_atpase"),
    ("Ca", r"\bCa2\+|\bCa\b|calcium", "cardiac_cicr_ode_model"),
    ("cortisol", r"cortisol", None),
]
PRECEDENCE = ["mitochondrial_oxphos", "na_k_atpase", "cardiac_cicr_ode_model"]

# One row per cell module. `decisive` and `status` are transcribed from the module's own
# selftest run (see FOLD_LEDGER.jsonl for the reproducing command). Status rule:
#   PROVEN   the run prints an explicit overall-pass token (all_gates_pass true, OVERALL PASS,
#            VERDICT CONFIRMED, REPRODUCES, SELFTEST PASS)
#   ASSUMED  the run prints a mixed summary (at least one gate FAIL/MISS) and still produces
#            its decisive number
#   REFUTED  the run reports that the claimed number is not reproduced
#   OPEN     the module is a shared library with no gates of its own
CELL_ROWS = [
 dict(sub="aging", mod="telomere_attrition", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Do measured leukocyte attrition rates and the Hayflick doubling arithmetic agree on a critical telomere length?",
      mech="literature-anchored arithmetic: rate x time integration plus linear regression over primary attrition rates",
      anchors="PMID 10432279, 10739676, 11413492, 11595186, 1438199",
      decisive="5/5 required and bonus gates pass: adult rate inside the 20-40 bp/year band, Hayflick arithmetic lands the critical length at 4-5 kb",
      hidden="critical telomere length at replicative senescence", regime="adult human leukocyte and fibroblast"),
 dict(sub="cardiovascular", mod="cardiac_cicr_ode_model", type="EMPIRICAL", risk="HIGH", status="ASSUMED",
      q="Does a time-resolved RyR2/SERCA/NCX ODE reproduce the measured cytosolic Ca transient of a ventricular myocyte?",
      mech="ODE model: 4-state Shannon-Bers RyR2 Markov gating, reversible Hill SERCA, phenomenological NCX",
      anchors="PMID 15347581, 8014907, 8488088, 8928885; DOI 10.1529/biophysj.104.047449",
      decisive="6/8 pre-registered gates PASS; the disclosed MISS is a ~7x over-estimate from the non-GHK trigger amplitude",
      hidden="junctional-cleft Ca concentration", regime="adult mammalian ventricular myocyte, 37 C"),
 dict(sub="cardiovascular", mod="ventricular_ap_core_tt04style", type="ENGINEERING", risk="LOW", status="OPEN",
      q="What reduced human-ventricular ionic ODE do the action-potential runners share?",
      mech="shared library: reduced ten Tusscher-Noble-Noble-Panfilov 2004 ionic model, no gates of its own",
      anchors="PMID 14656705",
      decisive="no gate output: the module is imported as a library and prints nothing when run",
      hidden="none (library)", regime="human ventricular tissue model"),
 dict(sub="energy", mod="enzyme_eyring_ceiling", type="THEORY", risk="LOW", status="PROVEN",
      q="Does the zero-free-parameter Eyring ceiling sit 5-8 decades above the fastest known enzyme turnover numbers?",
      mech="transition-state theory: kB*T/h prefactor plus an Eyring plot fitted to measured rate-temperature pairs",
      anchors="catalase, carbonic anhydrase II and ketosteroid isomerase turnover numbers as quoted in code",
      decisive="margin 5.19 decades (catalase), 6.79 (carbonic anhydrase II), 6.97 (KSI); Eyring fit R2 = 0.9986, shuffled-pairing void floor R2 = 0.218; 7/7 gates pass",
      hidden="activation free energy of the fastest attainable catalysed step", regime="298 K aqueous"),
 dict(sub="energy", mod="glycolytic_oscillations_selkov", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does an independent rebuild of the Selkov two-variable glycolytic oscillator produce a genuine limit cycle?",
      mech="2-variable ODE with allosteric product activation; Hopf analysis at the fixed point plus a stabilised void floor",
      anchors="Selkov 1968 formulation as given in code; Wolf 2000 period cited but not imported",
      decisive="nominal period 9.6225 dimensionless time units, bounded limit cycle; all gates pass, with the real-time rescaling disclosed as out of scope",
      hidden="existence and period of the limit cycle", regime="dimensionless Selkov normalisation"),
 dict(sub="energy", mod="mitochondrial_oxphos", type="EMPIRICAL", risk="HIGH", status="PROVEN",
      q="Do the proton-pumping and c-ring stoichiometries compose into the measured P/O ratios and membrane potential?",
      mech="stoichiometric gear-ratio model of complexes I-IV and ATP synthase, composed into a P/O grid",
      anchors="PMID 14977419, 20847295, 22392981; DOI 10.1016/j.bbabio.2004.09.004",
      decisive="mammalian c8/3 gear ratio gives H+/ATP = 3.667 and P/O(NADH) = 2.45-2.73; OVERALL PASS including the composite P/O-plus-membrane-potential falsifier",
      hidden="H+/ATP stoichiometry of the rotary machine", regime="mammalian mitochondria"),
 dict(sub="growth", mod="cytoskeleton_critical_concentration", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does two-end polymer kinetics give a lower barbed-end than pointed-end critical concentration, and does the asymmetry collapse without hydrolysis?",
      mech="two-end kon/koff polymer kinetics with a measured ADP-actin anchor and a noise-perturbation robustness test",
      anchors="Pollard 1986 J Cell Biol 103:2747 as quoted in code",
      decisive="Cc barbed 0.121 uM vs pointed 0.615 uM, ratio 5.10; hydrolysis-blocked ratio 1.0 at the measured 1.8 uM; all gates pass",
      hidden="end-specific critical concentrations", regime="Mg-ATP and Mg-ADP actin in vitro"),
 dict(sub="growth", mod="dna_replication_kinetic_proofreading_drake", type="EMPIRICAL", risk="LOW", status="PROVEN",
      q="Does Drake's cross-species mutations-per-genome invariant reproduce the stated E. coli per-base-pair fidelity?",
      mech="invariant divided by measured genome sizes, with a log-uniform redraw void floor",
      anchors="PMID 1831267",
      decisive="E. coli per-bp fidelity 7.112e-10 against the stated 7.11e-10; all gates pass",
      hidden="per-base-pair replication error rate", regime="microbial genomes, per replication"),
 dict(sub="growth", mod="genetic_code_error_minimization", type="THEORY", risk="MED", status="REFUTED",
      q="Does the standard codon table itself reproduce the claimed wobble redundancy, second-position clustering and optimisation ratio?",
      mech="deterministic enumeration of the codon table plus scrambled-relabelling and 2-opt void floors",
      anchors="the standard genetic code table as encoded in the module",
      decisive="wobble synonymous fraction 0.667 against the claimed 0.689 and eta2 = 0.1818 against the claimed 0.578-0.756: two headline gates fail, all_gates_pass false",
      hidden="error-minimisation margin of the standard code", regime="standard code, no experimental input"),
 dict(sub="growth", mod="somite_hes7_lewis_zeiser_dde_rebuild", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does a from-scratch delay-differential rebuild of the HES7 oscillator reproduce the segmentation-clock period?",
      mech="delay differential equations with Hill repression, numeric period plus an analytic cross-check and two void floors",
      anchors="PMID 16504083, 16432209",
      decisive="numeric period 120.120 min vs analytic 120.076 min vs anchor 120.0 min (0.04% cross-check); VERDICT CONFIRMED",
      hidden="intrinsic period of the delayed negative-feedback loop", regime="mouse somitogenesis"),
 dict(sub="growth", mod="spindle_sac_biomd186_rebuild", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does a rebuild of the spindle-assembly-checkpoint model clear the checkpoint inside the measured mitotic window?",
      mech="mass-action ODE rebuild of a curated checkpoint model with a never-attach void floor",
      anchors="PMID 18253502, 7642709",
      decisive="t90 brackets 1.650 to 48.286 min around the 23.0 min anchor, mass-conservation drift 1.05e-14; VERDICT CONFIRMED",
      hidden="checkpoint silencing time", regime="mammalian mitosis"),
 dict(sub="immune", mod="complement_cascade", type="EMPIRICAL", risk="HIGH", status="PROVEN",
      q="Does a reduced compartmental model of the three complement pathways discriminate activator from host surface?",
      mech="8-state compartmental ODE with autocatalytic C3b amplification and regulator-driven surface discrimination",
      anchors="PMID 14730550, 15544620, 16715088; DOI 10.1002/jcla.10095",
      decisive="unregulated vs regulated host membrane-attack-complex ratio 579.76 (gate > 2.0), forced adversary ratio 1.0; OVERALL PASS",
      hidden="surface-dependent amplification rate", regime="human plasma concentrations"),
 dict(sub="immune", mod="ferroptosis_cert", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Is ferroptosis pharmacologically orthogonal to apoptosis and necroptosis, with a lipid-peroxidation causal signature?",
      mech="machine-checked extraction gates over primary-source datapoints, no simulation",
      anchors="PMID 22632970, 24439385, 25402683, 27842066, 31634900",
      decisive="both required falsifiers pass, including the iron-specific potentiation of 3 sources against 0 of 4 other metals",
      hidden="death-mode identity of the lipid-peroxidation pathway", regime="cell-line and in-vivo studies as cited"),
 dict(sub="immune", mod="necroptosis_cert", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Is necroptosis a distinct RIP1/RIP3/MLKL-dependent death mode rather than generic toxicity?",
      mech="machine-checked extraction gates over primary-source datapoints, no simulation",
      anchors="PMID 16408008, 19498109, 19524512, 19524513",
      decisive="both required falsifiers pass across pharmacological and genetic dissection, including the kinase-activity requirement",
      hidden="death-mode identity of the necrosome pathway", regime="cell-line and in-vivo studies as cited"),
 dict(sub="immune", mod="tcell_activation_exhaustion", type="EMPIRICAL", risk="HIGH", status="ASSUMED",
      q="Do two-signal activation, clonal-expansion kinetics and graded exhaustion hold in one reduced model?",
      mech="kinetic-proofreading discrimination gate plus piecewise growth/death ODE cross-checked against LSODA",
      anchors="PMID 11602708, 11877489, 12663797, 16382236, 16917489",
      decisive="parts 1-4 all pass but the robustness sweep fails, so the run reports overall_pass False",
      hidden="functional exhaustion depth under chronic antigen", regime="mouse and human chronic-antigen studies"),
 dict(sub="musculoskeletal", mod="bone_calcium_flux_reconciliation", type="EMPIRICAL", risk="HIGH", status="ASSUMED",
      q="Does osteoclast activity converted bottom-up to resorbed mineral match the tracer-kinetic whole-skeleton calcium flux?",
      mech="two decorrelated conversion legs (histomorphometric and geometric) against a re-sourced tracer-kinetic anchor",
      anchors="PMID 12740946, 18988698; ICRP Publication 70 reference skeletal calcium",
      decisive="leg A 180.8-202.7 mg/day vs the 300-630 mg/day tracer anchor, gap ratio 1.48x, same sign as the rapid-exchange mechanism",
      hidden="daily resorbed mineral mass", regime="adult human skeleton"),
 dict(sub="musculoskeletal", mod="bone_rankl_opg_lemaire2004", type="EMPIRICAL", risk="MED", status="ASSUMED",
      q="Does an independent reimplementation of the RANKL/RANK/OPG ODE reproduce its own claimed baseline and perturbation response?",
      mech="from-scratch ODE reimplementation from the disclosed equations, free-run drift and Jacobian check",
      anchors="the disclosed Lemaire-2004 equations and parameters as quoted in the module",
      decisive="osteoclast peak +304.2% against the recorded +304%, 2000-day free-run drift below 0.02%; no overall-pass token is printed",
      hidden="osteoclast pool response to a RANKL stimulus", regime="dimensionless pM signalling units"),
 dict(sub="musculoskeletal", mod="meniscus_hoop_stress_model", type="EMPIRICAL", risk="HIGH", status="ASSUMED",
      q="Does the hoop-stress mechanism, not retained tissue mass, explain meniscal load transmission?",
      mech="closed-form hoop-tension model plus three forced-adversary tests against published measurements",
      anchors="PMID 11415635, 18762653, 37986646, 657636",
      decisive="3 of 3 adversary tests fall (root tear +25.0% peak pressure, hoop strain to ~0 when fibres are cut), but the held-out 50%-width test mismatches: predicted 1.0 vs measured 1.75",
      hidden="fraction of joint load carried as circumferential tension", regime="human and porcine knee, axial load"),
 dict(sub="musculoskeletal", mod="muscle_pcsa_crossbridge_specific_tension_rebuild", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does a molecular cross-bridge force budget predict whole-muscle specific tension?",
      mech="cross-bridge force times duty ratio times filament density, composed to a PCSA-normalised stress",
      anchors="PMID 39169839, 6511546, 8139653",
      decisive="specific tension 25.35 N/cm2 (sweep 14.32-41.57) against anchors 26.8 and 22.5 N/cm2, void floor max 8.31; VERDICT CONFIRMED",
      hidden="fraction of attached cross-bridges in isometric tetanus", regime="adult human skeletal muscle"),
 dict(sub="musculoskeletal", mod="reproduce_motor_unit_force_factorization", type="EMPIRICAL", risk="LOW", status="PROVEN",
      q="Does motor-unit count times mean twitch torque reproduce the recorded undershoot against maximal voluntary contraction?",
      mech="closed-form factorisation of count and twitch torque with a twitch-to-tetanus fusion mechanism check",
      anchors="PMID 15685623, 9415831",
      decisive="worst relative difference 0.2621% against a 1% tolerance; twitch:tetanus 0.25 predicts the stated ~75% undershoot; REPRODUCES",
      hidden="fusion factor between twitch and tetanic force", regime="tibialis anterior, young adults"),
 dict(sub="nervous", mod="corticospinal_motor_command", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does a population-vector decode of M1 activity beat the strongest fair labeled-line adversary?",
      mech="population-vector decode with Poisson trial noise on a dense held-out direction grid, plus a conduction-velocity gate",
      anchors="PMID 7143039, 8996498, 24872533",
      decisive="all gates pass: population vector under 15 degrees error at published cell counts and beating the labeled-line adversary at every N; overall_pass True",
      hidden="movement direction encoded in the population", regime="primate M1 reaching"),
 dict(sub="nervous", mod="hodgkin_huxley_squid_axon_ap_rebuild", type="EMPIRICAL", risk="LOW", status="PROVEN",
      q="Does a from-scratch Hodgkin-Huxley rebuild show all-or-none threshold, overshoot and a refractory period?",
      mech="four-variable conductance ODE integrated with bisected threshold search and a zero-gNa void floor",
      anchors="Hodgkin and Huxley 1952 J Physiol as quoted in the module",
      decisive="threshold 6.921 uA/cm2, peak 40.51 mV (residual +0.51 mV vs anchor), relative refractory window ~10-12 ms; VERDICT CONFIRMED",
      hidden="gating-variable kinetics behind the spike", regime="squid giant axon, 6.3 C"),
 dict(sub="nervous", mod="na_k_atpase", type="EMPIRICAL", risk="MED", status="ASSUMED",
      q="Do the Goldman-Hodgkin-Katz equation and measured pump kinetics jointly reproduce the resting membrane potential and the pump's share of the energy budget?",
      mech="GHK voltage equation over a concentration sweep plus Q10-corrected partial-reaction turnover arithmetic",
      anchors="PMID 18075585, 2410761, 12381678, 16439665; DOI 10.1172/JCI106703",
      decisive="6/8 pre-registered gates PASS, including the kidney-over-average ordinal gate at 40.0% vs 20.2%; the failing gate is diagnosed and retested",
      hidden="intracellular ion activities and per-site pump turnover", regime="adult mammalian neuron, 37 C"),
 dict(sub="nervous", mod="rhodopsin_thermal_dark_noise_arrhenius_gate", type="EMPIRICAL", risk="MED", status="ASSUMED",
      q="Does an Arrhenius barrier, rather than a small attempt frequency, explain the thermal dark-noise rate of rhodopsin?",
      mech="Arrhenius arithmetic relating the measured thermal rate, the photon energy and the implied prefactor, with a barrier-less void floor",
      anchors="PMID 26061742",
      decisive="thermal rate 1e-11/s with an activation energy of 22.0 kcal/mol implies a prefactor of 2.52e5/s; the naive-prior gate fails, so the verdict is CONFIRMED_WITH_DOCUMENTED_ANOMALY",
      hidden="activation barrier for spontaneous isomerisation", regime="vertebrate rod photoreceptor, body temperature"),
 dict(sub="nervous", mod="sleep_two_process_recovery_rebound", type="EMPIRICAL", risk="LOW", status="PROVEN",
      q="Does the two-process model predict the measured slow-wave rebound after sleep deprivation?",
      mech="single-exponential homeostatic process on a fixed sleep-wake schedule, derived independently of the recorded figure",
      anchors="PMID 6696142, 39458438",
      decisive="rebound +89.0 min, inside the anchor band and reached without the node's previously recorded +107.14 min figure; overall_pass true",
      hidden="homeostatic sleep pressure", regime="adult human, fixed deprivation schedule"),
 dict(sub="nervous", mod="synaptic_tsodyks_markram_ppr", type="EMPIRICAL", risk="LOW", status="PROVEN",
      q="Does the depression-only resource model reproduce the measured paired-pulse ratio?",
      mech="single depleting resource pool with recovery time constant and fixed release fraction, plus a no-dynamics void floor",
      anchors="PMID 30256194; Tsodyks and Markram 1997 PNAS as quoted",
      decisive="rebuild consistency gates all pass; the no-dynamics void floor gives a paired-pulse ratio of 1.0 in the wrong direction as required",
      hidden="release-ready resource fraction", regime="neocortical depressing synapse"),
 dict(sub="nervous", mod="synaptic_vesicle_release", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does a Ca-triggered SNARE release model reproduce the measured release probability and sub-millisecond latency?",
      mech="cooperative Ca-binding release model over a microdomain transient, composed to an excitatory postsynaptic current",
      anchors="PMID 10972290, 16794037, 16990140, 18046404",
      decisive="all five parts pass, including sub-millisecond action-potential-to-current latency; overall_pass True",
      hidden="local Ca concentration at the sensor", regime="central presynaptic terminal"),
 dict(sub="organ_systems", mod="countercurrent_multiplier", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does the single effect multiplied by loop geometry reproduce the corticomedullary osmotic gradient?",
      mech="iterated single-effect multiplier over a hairpin loop with washout, plus a geometry-free adversary",
      anchors="PMID 22237592, 22914749, 8760217",
      decisive="papilla 1254.7 mOsm from a 300 mOsm cortex, coupled-over-adversary gradient ratio 9.84; OVERALL_PASS true",
      hidden="axial interstitial osmolality profile", regime="mammalian kidney, antidiuresis"),
 dict(sub="organ_systems", mod="nephron_transport_cell", type="EMPIRICAL", risk="HIGH", status="ASSUMED",
      q="Does one nephron model hold filtration flat under a pressure sweep, span the urine-osmolality range and match measured clearance?",
      mech="Starling-divider filtration with a constant-elasticity autoregulator, swept over perfusion pressure",
      anchors="PMID 12791588, 15415454, 17728380",
      decisive="gate I fails at 102.0% worst-case deviation over 80-180 mmHg while holding the band over 90-120 mmHg; inulin clearance 122.8 mL/min; OVERALL_PASS false",
      hidden="pressure-dependent autoregulatory gain", regime="adult human kidney"),
 dict(sub="respiratory", mod="reproduce_alveolar_surface_factorization", type="EMPIRICAL", risk="LOW", status="PROVEN",
      q="Does the measured per-alveolus surface area exceed the equal-volume sphere prediction by the recorded factor?",
      mech="closed-form equal-volume-sphere comparison against morphometric counts and areas",
      anchors="PMID 14512270, 644146",
      decisive="ratio 2.3665 against the recorded 2.37, worst relative difference 0.1466% at a 1% tolerance; REPRODUCES",
      hidden="alveolar surface-to-volume amplification", regime="adult human lung morphometry"),
 dict(sub="signalling", mod="circadian_q10_compensation_goodwin", type="THEORY", risk="MED", status="REFUTED",
      q="Does a uniformly Q10-scaled Goodwin oscillator stay temperature compensated?",
      mech="3-variable Hill-repression oscillator with every rate scaled by one Q10, plus a broken-loop void floor",
      anchors="Goodwin 1965 formulation as encoded in the module",
      decisive="measured period Q10 = 6.303, outside the claimed compensation band, so the naive uniform-scaling null is not compensated; all_gates_pass false",
      hidden="temperature dependence of the loop period", regime="dimensionless Goodwin oscillator"),
 dict(sub="signalling", mod="mapk_erk_cascade", type="EMPIRICAL", risk="HIGH", status="PROVEN",
      q="Does the three-tier cascade convert a graded input into an ultrasensitive and, with feedback, bistable output?",
      mech="tiered covalent-modification ODE with Goldbeter-Koshland tiers and a feedback bistability sweep",
      anchors="PMID 6947258, 8816754, 9166761, 9228083, 10823939, 21768338",
      decisive="tier slopes multiply to 1.880 against a measured overall slope of 1.880 (0.000% relative error); overall_pass True",
      hidden="effective Hill coefficient of the cascade", regime="Xenopus extract and mammalian parameter sets as cited"),
 dict(sub="signalling", mod="mtorc1_feedback_rebuild", type="EMPIRICAL", risk="HIGH", status="PROVEN",
      q="Does a literal rebuild from the recorded topology reproduce the direction of the mTORC1 feedback responses?",
      mech="mass-action rebuild of the named feedback and drug edges over conserved pools; sign-level comparison only",
      anchors="the recorded constants and topology quoted in the module",
      decisive="SELFTEST PASS: all three sign checks agree (rapamycin +AKT, kinase-inhibitor S473 collapse, T308 over-recovery) while the magnitudes differ from the recorded values",
      hidden="strength of the S6K-IRS negative feedback", regime="sign-level only; the stoichiometry is under-determined by the recorded spec"),
 dict(sub="signalling", mod="smad_nucleocytoplasmic_shuttling_biomd173", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does a transcribed curated SMAD shuttling model reproduce the nuclear accumulation time course and dose response?",
      mech="compartmental ODE transcribed from a curated model file, with a gated half-maximum band",
      anchors="Schmierer, Tischer and Bird 2008 PNAS; the curated model file fetched as described in the module",
      decisive="overall_pass true, with the half-maximum outside-band check gated as designed",
      hidden="nuclear-cytoplasmic SMAD distribution", regime="HaCaT-like parameter set"),
 dict(sub="signalling", mod="wnt_betacatenin_rebuild", type="EMPIRICAL", risk="MED", status="PROVEN",
      q="Does a rebuild from the recorded destruction-complex equations reproduce the recorded beta-catenin steady state and the tankyrase null result?",
      mech="15-species mass-action destruction-complex ODE built only from the stated fluxes and constants",
      anchors="the recorded Lee-2003 lineage equations and constants quoted in the module",
      decisive="SELFTEST PASS: steady-state beta-catenin 1645.9144 unchanged to 2.3e-13 under a 1000-fold tankyrase change, matching the recorded 0.0% rescue",
      hidden="free beta-catenin steady state", regime="dimensionless model units"),
]

SUBSYSTEMS = ["aging", "cardiovascular", "energy", "growth", "immune",
              "musculoskeletal", "nervous", "organ_systems", "respiratory", "signalling"]


def module_sources():
    src = {}
    for sub in SUBSYSTEMS:
        for p in sorted((CELLS / sub).glob("*.py")):
            if p.name == "__init__.py":
                continue
            src[p.stem] = p.read_text(encoding="utf-8")
    return src


def import_edges(src):
    """R1: an edge for every intra-set module import."""
    edges = set()
    for mod, text in src.items():
        for name in re.findall(r"^\s*(?:from|import)\s+([\w\.]+)", text, re.M):
            leaf = name.split(".")[-1]
            if leaf in src and leaf != mod:
                edges.add((mod, leaf))
    return edges


def quantity_edges(src):
    """R2 + R3: shared-quantity edges, cycle-broken by producer precedence."""
    edges = set()
    census = {}
    for qname, pattern, producer in QUANTITIES:
        hits = sorted(m for m, text in src.items() if re.search(pattern, text))
        census[qname] = hits
        if producer is None:
            continue
        for consumer in hits:
            if consumer == producer:
                continue
            if consumer in PRECEDENCE and PRECEDENCE.index(producer) >= PRECEDENCE.index(consumer):
                continue          # R3: never point back up the supply order
            edges.add((consumer, producer))
    return edges, census


def build_catalog():
    src = module_sources()
    missing = {r["mod"] for r in CELL_ROWS} - set(src)
    if missing:
        raise SystemExit(f"cell table names modules that are not on disk: {sorted(missing)}")
    edges = import_edges(src)
    qedges, census = quantity_edges(src)
    edges |= qedges

    deps = {}
    for a, b in sorted(edges):
        deps.setdefault(a, []).append(b)

    clusters = {}
    for i, row in enumerate(CELL_ROWS, start=1):
        cell = {
            "node": {
                "id": row["mod"],
                "claim": f"{row['q']} Model: {row['mech']}. Decisive number from the cell's own run: {row['decisive']}.",
                "type": row["type"],
                "status": row["status"],
                "evidence": [] if row["status"] in ("OPEN", "DEFERRED") else [f"CELL-{i:03d}"],
                "depends_on": sorted(deps.get(row["mod"], [])),
                "risk": row["risk"],
                "regime_note": row["regime"],
            },
            "cert": {
                "hidden": row["hidden"],
                "occluded": row["hidden"] != "none (library)",
                "legs": [["model computation", [row["mod"]]], ["literature anchors", [row["anchors"]]]],
                "anchor": row["anchors"],
                "mode": "agreement",
                "_verify": {"verdict": row["status"], "note": row["decisive"]},
            },
        }
        clusters.setdefault(row["sub"], {"cells": []})["cells"].append(cell)

    for sub in SUBSYSTEMS:
        members = sorted(r["mod"] for r in CELL_ROWS if r["sub"] == sub)
        clusters[sub]["cells"].append({
            "node": {
                "id": f"{sub}_subsystem",
                "claim": f"Aggregate node for the {sub} subsystem package: the {len(members)} cells "
                         f"under it are the evidence this subsystem is modelled at all.",
                "type": "GOAL",
                "status": "OPEN",
                "evidence": [],
                "depends_on": members,
                "risk": "MED",
                "regime_note": f"stands over src/bodytwin/cells/{sub}/",
            },
            "cert": {"hidden": f"whole-subsystem behaviour of {sub}", "occluded": True,
                     "legs": [["member cells", members]], "anchor": "the member cells' own anchors",
                     "mode": "aggregate",
                     "_verify": {"verdict": "aggregate", "note": "no gate of its own"}},
        })

    catalog = {"foundational": [], "clusters": clusters}
    (HERE / "designs").mkdir(exist_ok=True)
    (HERE / "designs" / "cell_designs.json").write_text(
        json.dumps(catalog, indent=2) + "\n", encoding="utf-8")
    return census, edges


def write_ledger():
    lines = []
    src_index = {r["mod"]: r for r in CELL_ROWS}
    consumers = {}
    graph = json.loads((HERE / "ANCHOR_GRAPH.json").read_text(encoding="utf-8"))
    for n in graph["nodes"]:
        for dep in n["depends_on"]:
            consumers.setdefault(dep, []).append(n["id"])
    for i, row in enumerate(CELL_ROWS, start=1):
        path = f"src/bodytwin/cells/{row['sub']}/{row['mod']}.py"
        lines.append(json.dumps({
            "id": f"CELL-{i:03d}",
            "node": row["mod"],
            "claim": row["q"],
            "verdict": row["status"],
            "decisive_number": row["decisive"],
            "gate_cmd": f"BODYTWIN_OUT=./outputs .venv-bodytwin/bin/python {path}",
            "source": path,
            "consumers": sorted(consumers.get(row["mod"], [])),
            "supersedes": None,
        }, ensure_ascii=False))
    (HERE / "FOLD_LEDGER.jsonl").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    census, edges = build_catalog()
    env = dict(os.environ)
    env.update({
        "BODYTWIN_ROOT": str(HERE),
        "BODYTWIN_DESIGNS": str(HERE / "designs" / "cell_designs.json"),
        "BODYTWIN_GRAPH": str(HERE / "ANCHOR_GRAPH.json"),
        "BODYTWIN_INVENTORY": str(Path(tempfile.gettempdir()) / "bodytwin_cell_inventory.md"),
    })
    builder = REPO / "src" / "bodytwin" / "framework" / "bodytwin_build_anchor_graph.py"
    graph_path = HERE / "ANCHOR_GRAPH.json"
    if graph_path.exists():
        graph_path.unlink()          # the builder refuses to bootstrap over a populated graph
    r = subprocess.run([sys.executable, str(builder)], env=env, capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    sys.stderr.write(r.stderr)
    if r.returncode != 0:
        raise SystemExit(r.returncode)
    write_ledger()
    for q, hits in census.items():
        print(f"quantity {q}: {len(hits)} modules match")
    print(f"cell-to-cell edges: {len(edges)}")


if __name__ == "__main__":
    main()
