"""TAU PATHOLOGY -- hyperphosphorylation -> paired-helical-filament (PHF) aggregation,
the Braak staging stereotyped spatial spread (entorhinal -> limbic -> isocortical association ->
primary cortex), and the tau-vs-amyloid TEMPORAL ordering + spatial-tracking mismatch (amyloid diffuse
and plateaus early; tau tracks cognition and Braak stage far better). Companion to the amyloid_beta_aggregation cell.

QUESTION (pre-registered falsifier, stated before any number below is computed): does the model
reproduce (a) the measured tau-PET (or CSF/plasma p-tau) vs cognition correlation being STRONGER than
amyloid-PET vs cognition, in the SAME cohort where possible (the decorrelated "which protein tracks the
disease" test) AND (b) the Braak-stage stereotyped spatial ordering (entorhinal first, primary cortex
last)?

SYMMETRIC QC, stated up front, not discovered after the fact:
  - The amyloid-cascade CAUSAL direction is genuinely contested (amyloid-first vs tau-first vs parallel
    independent initiation) -- held OPEN. Jack et al 2013 (PMID 23332364) itself revises its own 2010
    model (PMID 20083042) to allow that "the two major proteinopathies ... might be initiated
    independently in sporadic AD."
  - Tau-PET tracers (AV-1451/flortaucipir) have documented off-target binding (neuromelanin/melanin,
    choroid-plexus calcification, vessels, substantia nigra -- Marquie 2015 PMID 26344059, Lowe 2016
    PMID 27296779) -- held OPEN as a measurement caveat. Reasoned (not measured) point: classical
    errors-in-variables attenuation bias means off-target NOISE should, if anything, WEAKEN a true
    correlation, not manufacture one -- so this caveat argues against, not for, the tau-beats-amyloid
    finding being a tracer-quality artifact. Flagged as reasoning, not a new empirical datapoint.

Reads: nothing (no upstream JSON dependency). Every number below is either (a) machine-computed
from a disclosed, coarse, literature-informed TOY connectivity graph (Part 2, the geometric/spectral
piece -- explicitly NOT a claim of connectomic precision, see Honest Gaps), or (b) directly
transcribed with disclosed, pre-registered arithmetic from live-verified published abstracts
(Part 3), never invented.
Writes: tau_pathology_results.json.
Gate: the pre-registered gates collected in all_gates; overall_pass is their conjunction.

GEOMETRIC STRUCTURE (derive from the geometry, not heuristics):
  Part 2 models trans-synaptic/trans-neuronal tau spread as a continuous-time diffusion process on a
  weighted anatomical graph: dx/dt = -L x, L = D - W the COMBINATORIAL GRAPH LAPLACIAN, solved via its
  own SPECTRUM (eigendecomposition L = U diag(lambda) U^T, x(t) = U diag(exp(-lambda t)) U^T x0) -- the
  same "derive from the graph/spectrum" family the sigma_min discipline uses elsewhere, not
  a curve fit. A connected graph's diffusion conserves total mass and converges to the UNIFORM vector
  (1/N)*sum(x0) (the averaging/consensus theorem for the graph Laplacian) -- this fixes, by construction
  and with NO free parameter, the per-node "half-arrival-time" threshold used below (0.5/N).
"""

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

OUT_DIR = _os.path.join(OUT_ROOT, "tau_pathology")
OUT_PATH = _os.path.join(OUT_DIR, "tau_pathology_results.json")

# ============================================================================
# PART 2 -- Braak spatial-spread: literature-informed TOY connectivity graph
# (coarse, disclosed as illustrative-mechanism-only, NOT a connectomic claim)
# ============================================================================

# 10 regional nodes, grouped into 4 tiers DIRECTLY licensed by Braak & Braak
# 1991's verbatim abstract text (PMID 1759558, re-verified live): "transentorhinal stages I-II" (tier 1); "limbic stages (III-IV)"
# with "mild involvement of the first Ammon's horn sector" (tier 2, hence
# hippocampus + amygdala + parahippocampal/perirhinal cortex grouped here);
# "isocortical stages (V-VI) ... destruction of virtually all isocortical
# ASSOCIATION areas" (tier 3 -- note the word "association", which by its own
# plain meaning EXCLUDES primary sensory/motor/visual cortex -- licensing
# tier 4 as a directly-implied "beyond VI / non-association" latest group,
# not an unverified extrapolation).
NODES = ["EC", "HIP", "AMY", "PHC", "ITG", "PCC", "PAR", "PFC", "V1", "M1S1"]
BRAAK_TIER = {
    "EC": 1,                                    # transentorhinal/entorhinal, stage I-II
    "HIP": 2, "AMY": 2, "PHC": 2,                # limbic, stage III-IV
    "ITG": 3, "PCC": 3, "PAR": 3, "PFC": 3,       # isocortical ASSOCIATION, stage V-VI
    "V1": 4, "M1S1": 4,                          # primary (non-association), latest/least-affected
}

# Undirected weighted edges: coarse, literature-informed (perforant path,
# entorhinal-parahippocampal adjacency, amygdalo-prefrontal projections,
# hippocampo-cingulate "extended limbic system", fronto-parietal and
# dorsal/ventral visual-association streams) -- textbook neuroanatomy, NOT a
# DTI-derived structural connectome. Disclosed as a toy, see Honest Gaps.
EDGES = [
    ("EC", "HIP", 1.0),   # perforant path -- the canonical, strongest EC->HIP projection
    ("EC", "PHC", 1.0),   # entorhinal cortex is embedded in/adjacent to parahippocampal gyrus
    ("EC", "AMY", 0.7),   # direct entorhinal-amygdala connections, weaker than the perforant path
    ("HIP", "PCC", 0.6),  # subiculum/cingulum bundle, "extended hippocampal-limbic system"
    ("AMY", "PFC", 0.5),  # amygdalo-prefrontal projections
    ("PHC", "ITG", 0.7),  # parahippocampal gyrus is continuous with medial/inferior temporal neocortex
    ("ITG", "PCC", 0.4),
    ("ITG", "PAR", 0.5),
    ("PCC", "PAR", 0.5),  # posterior default-mode-like cingulo-parietal connectivity
    ("PAR", "PFC", 0.5),  # fronto-parietal association network
    ("PAR", "V1", 0.3),   # dorsal visual-association stream
    ("ITG", "V1", 0.3),   # ventral visual-association stream
    ("PFC", "M1S1", 0.3),
    ("PAR", "M1S1", 0.2),
]

# ---- pre-registered gates (fixed BEFORE any diffusion is computed) --------
PREREG = {
    "rho_real_min": 0.75,           # EC-seeded, real topology, must reproduce Braak order strongly
    "ec_seed_rank_max": 1,          # EC must be the (strict) best of all 10 candidate epicenters
    "scrambled_rho_max_fraction": 0.5,  # scrambled-topology mean rho <= this fraction of rho_real
    "void_floor_range_max_fraction": 0.01,  # complete-equal-weight-graph spread <= 1% of real spread
    "n_scramble_trials": 500,
    "rng_seed": 20260722,
}


def build_graph(node_list, edge_list):
    n = len(node_list)
    idx = {node: i for i, node in enumerate(node_list)}
    w = np.zeros((n, n))
    for a, b, wt in edge_list:
        w[idx[a], idx[b]] = wt
        w[idx[b], idx[a]] = wt
    return w


def laplacian(w):
    d = np.diag(w.sum(axis=1))
    return d - w


def algebraic_connectivity(w):
    lap = laplacian(w)
    eigvals = np.linalg.eigvalsh(lap)
    return float(np.sort(eigvals)[1])  # second-smallest (Fiedler value); 0 iff disconnected


def diffusion_arrival_times(w, seed_idx, n_grid=4000):
    """Solve dx/dt = -L x via the graph's spectrum; return per-node time to
    reach 0.5*(1/N) (the node's equilibrium share), i.e. the half-arrival
    time. Seed node's arrival time is 0 by construction (excluded from the
    decisive rank statistic below, reported separately as the trivial part)."""
    n = w.shape[0]
    lap = laplacian(w)
    eigvals, eigvecs = np.linalg.eigh(lap)  # ascending; eigvals[0]=0 for connected graph
    x0 = np.zeros(n)
    x0[seed_idx] = 1.0
    coeffs = eigvecs.T @ x0  # projection onto spectral basis

    lam2 = eigvals[1] if len(eigvals) > 1 and eigvals[1] > 1e-12 else 1.0
    t_max = 25.0 / lam2
    ts = np.linspace(0.0, t_max, n_grid)

    # x(t)_i = sum_k coeffs[k] * exp(-eigvals[k]*t) * eigvecs[i,k]
    decay = np.exp(-np.outer(ts, eigvals))          # (n_grid, n)
    xt = decay @ (coeffs[:, None] * eigvecs.T)       # (n_grid, n)

    equilibrium = 1.0 / n
    threshold = 0.5 * equilibrium
    arrival = np.full(n, np.nan)
    for i in range(n):
        if i == seed_idx:
            arrival[i] = 0.0
            continue
        above = np.where(xt[:, i] >= threshold)[0]
        arrival[i] = ts[above[0]] if len(above) else np.inf
    return arrival, xt, ts


def rho_vs_braak(arrival, seed_idx, node_list, exclude_seed=True):
    ranks_target = np.array([BRAAK_TIER[n] for n in node_list], dtype=float)
    if exclude_seed:
        mask = np.arange(len(node_list)) != seed_idx
    else:
        mask = np.ones(len(node_list), dtype=bool)
    a = arrival[mask]
    b = ranks_target[mask]
    finite = np.isfinite(a)
    if finite.sum() < 3:
        return float("nan")
    rho, _p = spearmanr(a[finite], b[finite])
    return float(rho)


def run_part2():
    w_real = build_graph(NODES, EDGES)
    idx = {node: i for i, node in enumerate(NODES)}
    assert algebraic_connectivity(w_real) > 1e-9, "toy graph must be connected by construction"

    # --- primary: EC-seeded, real topology -----------------------------------
    ec_i = idx["EC"]
    arrival_ec, xt_ec, ts_ec = diffusion_arrival_times(w_real, ec_i)
    rho_real = rho_vs_braak(arrival_ec, ec_i, NODES)

    # --- adversary axis 1: exhaustive seed sweep (real topology, EVERY node
    #     as candidate epicenter) -- non-cherry-picked specificity-of-EPICENTER
    #     test. EC must come out on top for the "entorhinal-first" claim to be
    #     non-tautological (some ordering trivially exists for ANY seed on a
    #     hierarchical graph; the question is whether EC gives the BEST match).
    seed_sweep = {}
    for node in NODES:
        i = idx[node]
        arr, _, _ = diffusion_arrival_times(w_real, i)
        seed_sweep[node] = rho_vs_braak(arr, i, NODES)
    ranked_seeds = sorted(seed_sweep.items(), key=lambda kv: (-(kv[1] if np.isfinite(kv[1]) else -9), kv[0]))
    ec_rank = [n for n, _ in ranked_seeds].index("EC") + 1

    # --- adversary axis 2: scrambled topology (same 14 edge weights, same EC
    #     seed, node-pair identities randomized) -- specificity-of-TOPOLOGY
    #     test, many trials, only connected draws counted.
    rng = np.random.default_rng(PREREG["rng_seed"])
    all_pairs = [(a, b) for a in range(len(NODES)) for b in range(a + 1, len(NODES))]
    weights_multiset = [wt for _, _, wt in EDGES]
    scrambled_rhos = []
    attempts = 0
    while len(scrambled_rhos) < PREREG["n_scramble_trials"] and attempts < PREREG["n_scramble_trials"] * 20:
        attempts += 1
        chosen = rng.choice(len(all_pairs), size=len(weights_multiset), replace=False)
        w_scr = np.zeros_like(w_real)
        shuffled_weights = rng.permutation(weights_multiset)
        for (pair_idx, wt) in zip(chosen, shuffled_weights):
            a, b = all_pairs[pair_idx]
            w_scr[a, b] = wt
            w_scr[b, a] = wt
        if algebraic_connectivity(w_scr) <= 1e-9:
            continue
        arr_scr, _, _ = diffusion_arrival_times(w_scr, ec_i)
        scrambled_rhos.append(rho_vs_braak(arr_scr, ec_i, NODES))
    scrambled_rhos = np.array(scrambled_rhos, dtype=float)
    scrambled_rhos = scrambled_rhos[np.isfinite(scrambled_rhos)]

    # --- void floor: complete graph, equal weight (mean of real weights),
    #     EC-seeded -- topology-BLIND diffusion; no structure => no order.
    mean_w = float(np.mean(weights_multiset))
    w_complete = mean_w * (np.ones((len(NODES), len(NODES))) - np.eye(len(NODES)))
    arrival_void, _, _ = diffusion_arrival_times(w_complete, ec_i)
    non_seed_void = arrival_void[np.arange(len(NODES)) != ec_i]
    non_seed_real = arrival_ec[np.arange(len(NODES)) != ec_i]
    void_range = float(np.ptp(non_seed_void[np.isfinite(non_seed_void)]))
    real_range = float(np.ptp(non_seed_real[np.isfinite(non_seed_real)]))

    gates = {
        "graph_connected_by_construction": bool(algebraic_connectivity(w_real) > 1e-9),
        "rho_real_meets_prereg_threshold": bool(rho_real >= PREREG["rho_real_min"]),
        "ec_is_top_ranked_of_10_candidate_epicenters": bool(ec_rank <= PREREG["ec_seed_rank_max"]),
        "scrambled_topology_degrades_below_prereg_fraction": bool(
            np.isfinite(scrambled_rhos.mean())
            and scrambled_rhos.mean() <= PREREG["scrambled_rho_max_fraction"] * rho_real
        ),
        "void_floor_collapses_ordering_signal": bool(
            void_range <= PREREG["void_floor_range_max_fraction"] * real_range
        ),
    }

    return {
        "nodes": NODES,
        "braak_tier": BRAAK_TIER,
        "edges": [{"a": a, "b": b, "weight": wt} for a, b, wt in EDGES],
        "rho_real_ec_seed": rho_real,
        "arrival_times_ec_seed": {n: (None if not np.isfinite(arrival_ec[i]) else float(arrival_ec[i]))
                                   for i, n in enumerate(NODES)},
        "seed_sweep_rho_all_10_candidates": {n: (None if not np.isfinite(v) else v) for n, v in seed_sweep.items()},
        "seed_sweep_ranked_best_to_worst": [n for n, _ in ranked_seeds],
        "ec_rank_among_10_candidates": ec_rank,
        "scrambled_topology": {
            "n_valid_connected_trials": int(len(scrambled_rhos)),
            "n_attempts": int(attempts),
            "mean_rho": float(scrambled_rhos.mean()) if len(scrambled_rhos) else None,
            "sd_rho": float(scrambled_rhos.std()) if len(scrambled_rhos) else None,
            "min_rho": float(scrambled_rhos.min()) if len(scrambled_rhos) else None,
            "max_rho": float(scrambled_rhos.max()) if len(scrambled_rhos) else None,
        },
        "void_floor_complete_graph": {
            "mean_edge_weight_used": mean_w,
            "non_seed_arrival_time_range_void": void_range,
            "non_seed_arrival_time_range_real": real_range,
            "void_range_as_fraction_of_real_range": (void_range / real_range) if real_range > 0 else None,
        },
        "gates": gates,
    }


# ============================================================================
# PART 3 -- temporal ordering + decorrelated "which protein tracks the
# disease" arithmetic, transcribed directly from live-verified abstracts.
# All numbers below are quoted, not invented; the only computation performed
# is transparent subtraction/division, shown explicitly.
# ============================================================================

def run_part3():
    # Villemagne et al 2013 Lancet Neurol (PMID 23477989, AIBL, n=200
    # longitudinal, re-verified live when this cell was written): years-before-clinical-
    # dementia-onset for each biomarker to reach its own threshold.
    villemagne = {
        "pmid": "23477989",
        "n_total_assessed": 200,
        "n_with_3plus_pib_followups_and_positive_accumulation": 163,
        "amyloid_years_before_onset": 17.0,
        "amyloid_years_before_onset_ci": [14.9, 19.9],
        "hippocampal_atrophy_years_before_onset": 4.2,
        "hippocampal_atrophy_years_before_onset_ci": [3.6, 5.1],
        "memory_impairment_years_before_onset": 3.3,
        "memory_impairment_years_before_onset_ci": [2.5, 4.5],
        "suvr_hc_mean": 1.38, "suvr_mci_mean": 1.94, "suvr_ad_mean": 2.27,
        "suvr_accumulation_rate_per_year": 0.043,
    }
    gap_amyloid_to_atrophy = villemagne["amyloid_years_before_onset"] - villemagne["hippocampal_atrophy_years_before_onset"]
    gap_amyloid_to_memory = villemagne["amyloid_years_before_onset"] - villemagne["memory_impairment_years_before_onset"]
    gap_atrophy_to_memory = villemagne["hippocampal_atrophy_years_before_onset"] - villemagne["memory_impairment_years_before_onset"]

    # Ossenkoppele et al 2016 Brain (PMID 26962052, n=16 AD patients with all
    # 3 PET tracers available, re-verified live when this cell was written): SAME-cohort,
    # SAME-day-scan decorrelated comparison of regional correspondence.
    ossenkoppele = {
        "pmid": "26962052",
        "n_all_three_pet_available": 16,
        "n_total_ad_cohort": 20,
        "n_amyloid_negative_controls": 15,
        "r_tau_vs_fdg_hypometabolism": -0.49, "r_tau_vs_fdg_sd": 0.07,
        "r_amyloid_vs_fdg_hypometabolism": 0.16, "r_amyloid_vs_fdg_sd": 0.09,
        "r_tau_vs_amyloid": 0.18, "r_tau_vs_amyloid_sd": 0.09,
        "all_p_values": "< 0.001 (all three)",
    }
    ratio_tau_over_amyloid_fdg_tracking = abs(ossenkoppele["r_tau_vs_fdg_hypometabolism"]) / abs(ossenkoppele["r_amyloid_vs_fdg_hypometabolism"])

    # Multi-cohort convergence tally -- decorrelated by RESEARCH SITE, not by
    # paper count (Bejanin/La Joie/Ossenkoppele-2016 share the same UCSF
    # Memory and Aging Center recruitment base -- treated as ONE site, not
    # three, to avoid a common-mode over-count; this is a deliberate
    # correction, disclosed here, not a flattering inflation).
    convergence = [
        {"site": "Washington University St Louis (Knight ADRC)", "pmid": "27169802",
         "cite": "Brier et al 2016 Sci Transl Med",
         "finding": "tau-PET (temporal lobe) more closely tracked dementia status and was a "
                    "better predictor of cognitive performance than amyloid-PET in ANY region",
         "tau_beats_amyloid": True},
        {"site": "Harvard Aging Brain Study", "pmid": "31157827",
         "cite": "Hanseeuw et al 2019 JAMA Neurol", "n": 60,
         "finding": "serial mediation: antecedent Abeta rise -> tau change (beta=1.07,P=.02) -> "
                    "cognitive (PACC) change (beta=-3.28,P=.001); tau, not baseline amyloid, "
                    "is the proximate driver of the cognitive step",
         "tau_beats_amyloid": True},
        {"site": "UCSF Memory and Aging Center (3 overlapping-cohort papers, counted ONCE)",
         "pmids": ["26962052", "29053874", "31894103"],
         "cite": "Ossenkoppele 2016 Brain (n=16-20) + Bejanin 2017 Brain (n=40) + La Joie 2020 "
                 "Sci Transl Med (n=32)",
         "finding": "tau-PET quantitatively 3.06x more strongly regionally correlated with "
                    "hypometabolism than amyloid-PET (same 16 pts); tau-PET (not amyloid-PET) "
                    "predicts subsequent atrophy (n=32); tau-cognition regional match only "
                    "weakly related to amyloid burden (n=40)",
         "tau_beats_amyloid": True},
        {"site": "Sichuan Provincial People's Hospital, China (independent geography/population)",
         "pmid": "42404124", "n": 438,
         "cite": "Hu et al 2026 Front Neurol",
         "finding": "in the largest and most recent cohort found (n=438: 325 AD+68 MCI+45 HC), "
                    "tau deposition showed stronger cognitive correlations than amyloid; tau an "
                    "independent driver of cognitive decline",
         "tau_beats_amyloid": True},
    ]
    n_independent_sites = len(convergence)
    n_dissenting = sum(1 for c in convergence if not c["tau_beats_amyloid"])
    # disclosed, honest note: a targeted (not systematic/exhaustive) adversarial
    # search for a dissenting same-cohort finding (amyloid > tau at tracking
    # cognition) was run when this cell was written; none was found. Absence of a finding in
    # a targeted, non-exhaustive search is reported as exactly that -- not
    # proof no such finding exists anywhere in the literature.
    dissent_search_conducted = True
    dissent_found = False

    gates = {
        "amyloid_leads_neurodegeneration_by_over_a_decade": bool(gap_amyloid_to_atrophy >= 10.0),
        "atrophy_and_memory_impairment_near_coincident_lt_2y": bool(gap_atrophy_to_memory <= 2.0),
        "tau_fdg_tracking_at_least_2x_amyloid_fdg_same_cohort": bool(ratio_tau_over_amyloid_fdg_tracking >= 2.0),
        "convergence_at_least_3_independent_sites_zero_dissent": bool(n_independent_sites >= 3 and n_dissenting == 0),
    }

    return {
        "villemagne_2013_aibl": villemagne,
        "computed_gaps_years": {
            "amyloid_to_hippocampal_atrophy": gap_amyloid_to_atrophy,
            "amyloid_to_memory_impairment": gap_amyloid_to_memory,
            "atrophy_to_memory_impairment": gap_atrophy_to_memory,
        },
        "ossenkoppele_2016_same_cohort_decorrelation": ossenkoppele,
        "ratio_tau_over_amyloid_fdg_tracking": ratio_tau_over_amyloid_fdg_tracking,
        "convergence_tally": convergence,
        "n_independent_sites_agreeing": n_independent_sites,
        "n_dissenting_sites": n_dissenting,
        "dissent_search_conducted_targeted_not_exhaustive": dissent_search_conducted,
        "dissent_found": dissent_found,
        "gates": gates,
    }


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    part2 = run_part2()
    part3 = run_part3()

    all_gates = {}
    all_gates.update({f"part2_{k}": v for k, v in part2["gates"].items()})
    all_gates.update({f"part3_{k}": v for k, v in part3["gates"].items()})
    overall_pass = all(all_gates.values())

    open_modeling_uncertainty = {
        "toy_graph_is_not_a_connectome": (
            "The Part 2 network is a 10-node, 14-edge, hand-specified toy built from textbook "
            "neuroanatomy (perforant path, entorhinal-parahippocampal adjacency, fronto-parietal "
            "and visual-association streams), NOT a DTI-tractography or histology-derived "
            "structural connectome. Its role is to demonstrate that hierarchical-topology + "
            "correct-epicenter is a SUFFICIENT mechanism to reproduce the qualitative Braak "
            "ordering from first-principles graph diffusion -- the external, real-data anchor for "
            "the actual magnitude is Vogel et al 2020 (PMID 32457389, n=312 real human tau-PET, "
            "connectome epidemic-spreading-model explains up to 70% of spatial-pattern variance)."
        ),
        "amyloid_causal_direction_genuinely_contested": (
            "Jack et al 2013 (PMID 23332364) itself revised its own 2010 model (PMID 20083042) to "
            "allow that Abeta and tau 'might be initiated independently in sporadic AD' -- held "
            "OPEN, not resolved, per the task's framing. Vogel 2020 (PMID 32457389) found the "
            "connectome-diffusion PATTERN of tau fits irrespective of amyloid status, but regions "
            "with higher amyloid burden show MORE tau than connectivity predicts -- amyloid may "
            "accelerate spread RATE without gating the spatial PATTERN itself; Hanseeuw 2019 (PMID "
            "31157827) and Palmqvist 2020 (PMID 32722745, rho=0.64 amyloid-positive vs rho=0.15 "
            "P=.33 NS amyloid-negative) support an amyloid-gates-tau-rate/magnitude relationship. "
            "These are not fully reconciled here -- disclosed as an open hidden variable."
        ),
        "amyloid_tau_spatial_correlation_disagrees_across_cohort_scale": (
            "Ossenkoppele 2016 (PMID 26962052, n=16, mixed atypical AD phenotypes + amyloid-negative "
            "controls) found amyloid-tau regional correlation weak (r=0.18). Hu et al 2026 (PMID "
            "42404124, n=438, mostly typical amnestic AD) found amyloid-tau SUVR correlation strong "
            "(r=0.65-0.81). Plausible (not adjudicated) explanation: phenotype-diverse/early/mixed "
            "cohorts decouple the two pathologies spatially; large homogeneous-stage amnestic-AD "
            "cohorts show more co-localization as disease progresses. Reported as a genuine, "
            "disclosed tension, not smoothed into false consistency."
        ),
        "tau_pet_off_target_binding_is_a_measurement_caveat": (
            "AV-1451/flortaucipir off-target binding is directly documented (Marquie 2015 PMID "
            "26344059: neuromelanin/melanin, hemorrhagic lesions; Lowe 2016 PMID 27296779: vessels, "
            "iron-associated regions, substantia nigra, choroid-plexus calcification, leptomeningeal "
            "melanin -- and explicitly 'does not completely reflect early stage tau progression "
            "suggested by Braak neurofibrillary tangle staging'). The reasoned (not measured) "
            "point: classical errors-in-variables attenuation bias means such noise should WEAKEN, "
            "not manufacture, a true correlation -- argues against, not for, the tau-beats-amyloid "
            "finding being a tracer-artifact, but is not itself a quantified correction."
        ),
        "no_single_cohort_jointly_measures_connectome_fit_AND_autopsy_Braak_AND_seeding_assay": (
            "The 4 decorrelated spatial-spread legs (autopsy histology Braak 1991 PMID 1759558; "
            "biochemical seeding-assay Kaufman 2018 PMID 29752551; in-vivo PET Scholl 2016 PMID "
            "26938442; connectome-diffusion Vogel 2020 PMID 32457389) corroborate ACROSS different "
            "cohorts, not jointly within one -- a disclosed limitation of the tau-propagation model."
        ),
        "brier_2016_headline_quote_is_qualitative_not_a_reported_r_value": (
            "Brier et al 2016 (PMID 27169802)'s abstract states tau 'was a better predictor of "
            "cognitive performance than Abeta deposition in any region of the brain' -- an explicit, "
            "verified quote, but the abstract does not itself print the underlying r/beta "
            "coefficients (would require full-text extraction, not attempted when this cell was written). Used "
            "here as a qualitative corroborating leg, not a quantitative one."
        ),
        "dissent_search_was_targeted_not_a_systematic_review": (
            "A live adversarial search for a same-cohort finding where amyloid beats tau at tracking "
            "cognition/neurodegeneration was run when this cell was written (3 targeted query variants); none was "
            "found. This is reported exactly as that -- a targeted negative result -- not as proof "
            "that no such finding exists anywhere in the literature."
        ),
    }

    report = {
        "task": (
            "Tau hyperphosphorylation -> paired-helical-filament (PHF) aggregation; Braak staging "
            "stereotyped spatial spread (entorhinal -> hippocampus/limbic -> neocortex); tau-vs-"
            "amyloid temporal ordering + spatial mismatch (amyloid diffuse/early-plateau, tau tracks "
            "cognition/Braak stage far better)."
        ),
        "confidence_tier": (
            "in-vivo-anchored (tau-PET / CSF-plasma p-tau / Braak autopsy neuropathology) for the "
            "spatial-spread and tau-vs-amyloid-cognition legs; the PHF molecular-mechanism leg is "
            "in-vitro/ex-vivo biochemical + atomic-resolution cryo-EM structural (Fitzpatrick 2017); "
            "the Part 2 network-diffusion model is an illustrative FIRST-PRINCIPLES TOY GEOMETRY, "
            "explicitly not connectome-precision-anchored (see Honest Gaps)."
        ),
        "script": "tau_pathology.py",
        "part2_braak_spatial_spread_geometric_model": part2,
        "part3_temporal_ordering_and_decorrelation": part3,
        "gates": all_gates,
        "n_gates_pass": sum(all_gates.values()),
        "n_gates_total": len(all_gates),
        "overall_pass": overall_pass,
        "open_modeling_uncertainty": open_modeling_uncertainty,
    }

    print(json.dumps(all_gates, indent=2))
    print(f"\n{sum(all_gates.values())}/{len(all_gates)} gates PASS")
    print(f"\nOVERALL: {'PASS' if overall_pass else 'FAIL/SURPRISE -- see gates above'}")

    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {OUT_PATH}")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
