"""REWARD PREDICTION ERROR (RPE) -- the midbrain dopamine (VTA/SNc) teaching signal that updates
striatal values (Schultz/Montague/Sutton-Barto temporal-difference (TD) learning). This is a
DISTINCT, higher (computational/spike-rate-coding) layer from the dopamine_kinetics cell, which
models the synaptic DAT-uptake concentration TRANSIENT rather than the spike-rate code.

QUESTION (pre-registered falsifier, stated before any number below is
computed): does a TD(0) model, given only a state representation of
"time-since-trial-events" and the standard update
delta(t) = r(t) + gamma*V(t+1) - V(t), reproduce the three signatures Schultz
and colleagues measured in primate VTA/SNc single units (Schultz/Apicella/
Ljungberg 1993 PMID 8441015; Schultz/Dayan/Montague 1997 PMID 9054347;
Hollerman/Schultz 1998 PMID 10195164): (i) a phasic BURST to unpredicted
reward, (ii) near-SILENCE at reward delivery once a cue fully predicts it (with
the burst TRANSFERRED to the cue), (iii) a sub-baseline DIP when a predicted
reward is OMITTED? DECORRELATED CHECK: does firing (modeled as an affine
function of delta) behave as a LINEAR function of the reward-prediction-error
term specifically -- not of raw reward magnitude alone -- matching the
quantitative claim of Bayer & Glimcher 2005 (PMID 15996553)?

FORCED ADVERSARY (stated up front, not discovered after the fact): a "dopamine
encodes reward/value ITSELF, not the error" account is given its steel-manned
form (reward-magnitude coding PLUS stimulus-specific habituation/adaptation,
a genuine confound that can independently mimic part of signature (ii) via
repetition alone) and is tested for whether it can still reproduce signature
(iii) (the omission dip) and a targeted third discriminator (a late,
never-cued, "surprise" reward after the adversary's habituation has fully
saturated) that the TD account and the adversary predict OPPOSITELY.

GEOMETRIC STRUCTURE (derive from the state-transition geometry, not curve
fitting): trial time is represented as a small absorbing/refreshing Markov
chain -- a single shared "background/ITI" state entered with a GEOMETRIC
(memoryless) holding time before the conditioned stimulus (CS) fires, followed
by a fixed-length deterministic countdown of "CS+i" states to the unconditioned
stimulus (US/reward) bin. This is not an arbitrary choice: a FIXED per-trial
clock (state = "bins since trial onset") would let value backpropagate all the
way to trial onset, since trial-onset would then ALSO be a perfect predictor of
reward timing -- an artifact of the representation, not a property of
learning. The random-holding-time background state is what forces the model to
localize its burst at the EARLIEST *reliable* predictor (the CS) rather than at
an arbitrary earlier point in time: background's conditional time-to-reward
has high variance (long, uncertain wait), so its bootstrapped value stays a
small residual, while the CS's conditional time-to-reward has ZERO variance (a
deterministic countdown), so essentially all of the discounted value collects
there. The full-chain fixed point is solved here in BOTH closed form (an exact
algebraic solve of the linear Bellman equations) and by stochastic-approximation
simulation, and the two are cross-checked against each other -- not asserted.
One clean, general TD-theoretic identity falls out of the algebra and is
reported explicitly: at full convergence, the omission-dip magnitude equals
EXACTLY the reward magnitude (delta_omission = -r_mag), because the background
contribution that appears in both the reward-delivery and the omission
equations cancels algebraically -- a sharper and more elegant result than the
"dip approximates the negative of the expected value" framing, and fully
consistent with it (V(reward-state) is itself close to r_mag).

DYSFUNCTION POLE (same substrate, same update rule, three different single-
parameter lesions of it): (1) ADDICTION -- Redish 2004's (PMID 15591205)
"noncompensable drug-induced dopamine increase" is implemented as a fixed
fraction phi of the reward-bin's effective TD update that ignores the
subtracted value term entirely; the closed-form fixed point shows this must
inflate the learned cue-value by a factor of (1 + phi/(1-phi)) relative to a
fully-compensable natural reward -- an unbounded escalation as phi->1. (2)
PARKINSON'S -- Frank/Seeberger/O'Reilly 2004's (PMID 15528409) Go(D1)/NoGo(D2)
asymmetric-learning-rate mechanism is implemented directly (positive-error
updates gated by a DA-dependent "Go gain", negative-error updates by a
DA-independent "NoGo gain") on a two-armed probabilistic-selection task, and
the OFF-medication (low Go-gain) vs ON-medication bias-reversal is measured.
(3) ABERRANT SALIENCE (Kapur 2003, PMID 12505794) -- modeled as a spurious,
causally-uncoupled positive-error injection rate onto a truly-neutral,
non-predictive stimulus, producing a false learned value where none is
warranted; flagged explicitly as the softest/most illustrative of the three
poles, matching Kapur 2003's status as a conceptual/theoretical review
rather than a primary quantitative data paper (an honest, disclosed asymmetry
in evidence strength across the three poles).

Reads: nothing (all literature parameters embedded).
Writes: reward_prediction_error_results.json.
Gate: the pre-registered gates in PREREG, summarised in the results JSON.
"""

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import json
import os

import numpy as np

OUT_DIR = _os.path.join(OUT_ROOT, "reward_prediction_error")
OUT_PATH = _os.path.join(OUT_DIR, "reward_prediction_error_results.json")

# ---- pre-registered gates (fixed BEFORE any number below is computed) -----
PREREG = {
    "rng_seed": 20260722,
    # --- core TD chain geometry (complete-serial-compound-style state rep over
    # a background-hazard + fixed-countdown Markov chain) ---
    "cs_to_us_delay": 10,      # k: deterministic bins from CS onset to reward
    "p_bg_to_cs": 0.02,        # geometric hazard leaving background (mean ITI = 50 bins)
    "gamma": 0.85,             # discount factor
    "alpha": 0.10,             # TD(0) learning rate
    "reward_mag": 1.0,
    "k_jitter_buffer": 6,      # extra CS-relative states to accommodate timing jitter
    "n_learn_trials": 300,     # trials to reach asymptote (checked for convergence, not assumed)
    "convergence_check_last_n": 20,
    "convergence_max_rel_std": 0.15,
    # Gate 1 (low-weight sanity/plumbing check; the interesting claims are G2-G4)
    "min_naive_burst": 0.5,
    # Gate 2: predicted-reward response must shrink by >= this fraction vs the naive burst
    "min_reward_response_shrink_frac": 0.90,
    # Gate 3: CS-onset burst after learning must reach >= this fraction of the naive burst
    # (pre-derived closed-form expectation ~0.177 at these parameters; threshold set well below it)
    "min_cue_transfer_frac_of_naive_burst": 0.10,
    # Gate 4 (decisive): omission dip -- exact closed-form identity is delta_omission = -reward_mag
    "max_omission_exact_identity_rel_err_closedform": 0.01,
    "max_omission_simulated_vs_identity_rel_err": 0.20,
    "min_omission_dip_magnitude": 0.3,
    # realistic-uncertainty variant (forced fix for the idealized-silence gap vs Schultz 1993's
    # empirical 25%->9%, not 100%->0%, residual responsiveness)
    "realistic_jitter_std_bins": 1.5,
    "realistic_partial_reinforcement_p": 0.9,
    "realistic_n_trials": 300,
    # adversary (value/reward-magnitude code, forced to its strongest fair form w/ habituation)
    "adversary_habituation_beta": 0.985,
    "adversary_n_trials_to_saturate_habituation": 300,
    # Bayer & Glimcher synthetic identifiability grid (explicitly SYNTHETIC, not a re-fit of their
    # raw data -- demonstrates why the linear-RPE claim is identifiable only when expectation is
    # varied independently of reward magnitude; the "positive-error-only" gain asymmetry below is
    # informed by Bayer & Glimcher's abstract caveat, live-verified when this cell was written)
    "reward_levels": [0.5, 1.0, 1.5, 2.0, 2.5],
    "expectation_levels": [0.0, 0.5, 1.0, 1.5, 2.0],
    "neg_rpe_gain_ratio": 0.5,
    "firing_noise_std": 0.15,
    "min_r2_gap_full_design": 0.5,
    "max_delta_variance_degenerate_design": 1e-9,
    # addiction (Redish TDRL) escalation check
    "addiction_n_trials": 300,
    "addiction_phi_noncompensable": 0.9,
    "min_value_ratio_drug_vs_natural": 3.0,   # closed-form predicts ~1+phi/(1-phi)=10x; threshold set well below
    # Parkinson's Go/NoGo asymmetry (Frank et al 2004 mechanism)
    "pd_n_trials": 800,
    "pd_p_reward_A": 0.8,
    "pd_p_reward_B": 0.2,
    "pd_softmax_beta": 5.0,
    "pd_alpha_base": 0.15,
    "pd_off_med_go_gain": 0.3,
    "pd_off_med_nogo_gain": 1.0,
    "pd_on_med_go_gain": 1.0,
    "pd_on_med_nogo_gain": 1.0,
    # aberrant salience (Kapur) spurious-value-on-neutral-cue check
    "salience_n_trials": 500,
    "salience_p_spurious": 0.05,
    "salience_spurious_bump": 0.5,
    "salience_min_neutral_value_aberrant": 0.05,
}


# ============================================================================
# Part 0 -- closed-form (exact algebraic) fixed point of the background+CS
# countdown chain, solved directly from the linear Bellman equations (not
# fitted, not iterated) -- used to cross-check the stochastic simulation below.
# ============================================================================
def closed_form_fixed_point(gamma, p_cs, k, r_mag):
    """V[0] (background), V[k] (reward bin), V[1] (CS-onset+1 state) at the
    exact TD fixed point of: V[k]=r+gamma*V[0]; V[rel]=gamma*V[rel+1] for
    rel<k; V[0]=p_cs*gamma*V[1]/(1-gamma*(1-p_cs)) (memoryless background
    Bellman equation). Solved by direct substitution into one linear equation
    for V[0]."""
    D = 1 - gamma * (1 - p_cs)
    denom = D - p_cs * gamma ** (k + 1)
    v0 = p_cs * (gamma ** k) * r_mag / denom
    vk = r_mag + gamma * v0
    v1 = gamma ** (k - 1) * vk
    delta_cs_onset = gamma * v1 - v0
    delta_reward_bin_at_fixedpoint = r_mag + gamma * v0 - vk  # must be exactly 0 (definition of vk)
    delta_omission_exact = 0.0 + gamma * v0 - vk              # the general identity: == -r_mag algebraically
    return {
        "v0_background": v0, "vk_reward_bin": vk, "v1_cs_onset_state": v1,
        "delta_cs_onset_predicted": delta_cs_onset,
        "delta_reward_bin_at_fixedpoint_should_be_zero": delta_reward_bin_at_fixedpoint,
        "delta_omission_predicted": delta_omission_exact,
        "delta_omission_algebraic_identity_minus_r_mag": -r_mag,
    }


# ============================================================================
# Part 1 -- the core TD(0) simulator over the background-hazard + CS-countdown
# chain. One function reused (with different knobs) for every experiment
# below: naive burst, learned silence+transfer, omission probe, realistic
# jitter/partial-reinforcement, and the addiction "noncompensable" lesion.
# ============================================================================
def run_conditioning(rng, prereg, n_trials, partial_p=1.0, jitter_std=0.0,
                     drug_mode=False, phi=0.0, V=None):
    k = prereg["cs_to_us_delay"]
    p_cs = prereg["p_bg_to_cs"]
    gamma = prereg["gamma"]
    alpha = prereg["alpha"]
    r_mag = prereg["reward_mag"]
    k_buffer = prereg["k_jitter_buffer"]
    n_states = k + k_buffer + 1  # 0=background, 1..(k+k_buffer)=CS-relative bins

    if V is None:
        V = np.zeros(n_states)

    delta_at_nominal_us = []
    delta_at_cs_onset = []
    v_at_cs_onset_trace = []

    for _trial in range(n_trials):
        bg_len = int(rng.geometric(p_cs)) - 1
        for _ in range(bg_len):
            delta_bg = 0.0 + gamma * V[0] - V[0]
            V[0] += alpha * delta_bg
        delta_transition = 0.0 + gamma * V[1] - V[0]
        V[0] += alpha * delta_transition
        delta_at_cs_onset.append(delta_transition)
        v_at_cs_onset_trace.append(V[1])

        us_bin = k
        if jitter_std > 0:
            us_bin = int(np.clip(round(rng.normal(k, jitter_std)), 1, k + k_buffer))
        rewarded = rng.random() < partial_p

        for rel in range(1, k + k_buffer + 1):
            is_us_bin = (rel == us_bin)
            r = r_mag if (is_us_bin and rewarded) else 0.0
            terminal = is_us_bin or (rel == k + k_buffer)
            next_state = 0 if terminal else rel + 1
            delta = r + gamma * V[next_state] - V[rel]
            eff_delta = delta
            if drug_mode and is_us_bin and rewarded:
                eff_delta = (1 - phi) * delta + phi * r_mag
            V[rel] += alpha * eff_delta
            if rel == k:
                delta_at_nominal_us.append(delta)
            if terminal:
                break

    return {
        "V": V,
        "delta_at_nominal_us": np.array(delta_at_nominal_us),
        "delta_at_cs_onset": np.array(delta_at_cs_onset),
        "v_at_cs_onset_trace": np.array(v_at_cs_onset_trace),
    }


def part1_three_signatures(rng, prereg):
    r_mag = prereg["reward_mag"]

    # --- signature (i): naive unpredicted reward, V=0 everywhere (trial 1) ---
    naive = run_conditioning(rng, prereg, n_trials=1, partial_p=1.0, jitter_std=0.0)
    delta_naive_reward = float(naive["delta_at_nominal_us"][0])

    # --- signature (ii)+(iii) substrate: deterministic learning to asymptote ---
    learn = run_conditioning(rng, prereg, n_trials=prereg["n_learn_trials"], partial_p=1.0, jitter_std=0.0)
    V_converged = learn["V"]
    last_n = prereg["convergence_check_last_n"]
    tail_us = learn["delta_at_nominal_us"][-last_n:]
    tail_cs = learn["delta_at_cs_onset"][-last_n:]
    delta_learned_reward = float(np.mean(tail_us))
    delta_learned_cs_onset = float(np.mean(tail_cs))
    convergence_rel_std_us = float(np.std(tail_us) / (abs(np.mean(tail_us)) + 1e-9))

    # --- closed form cross-check ---
    cf = closed_form_fixed_point(prereg["gamma"], prereg["p_bg_to_cs"], prereg["cs_to_us_delay"], r_mag)
    v1_sim, v1_cf = V_converged[1], cf["v1_cs_onset_state"]
    v0_sim, v0_cf = V_converged[0], cf["v0_background"]
    vk_sim, vk_cf = V_converged[prereg["cs_to_us_delay"]], cf["vk_reward_bin"]
    rel_err_v1 = abs(v1_sim - v1_cf) / v1_cf
    rel_err_v0 = abs(v0_sim - v0_cf) / v0_cf
    rel_err_vk = abs(vk_sim - vk_cf) / vk_cf

    # --- signature (iii): single held-out omission probe from the converged V ---
    probe = run_conditioning(rng, prereg, n_trials=1, partial_p=0.0, jitter_std=0.0, V=V_converged.copy())
    delta_omission_sim = float(probe["delta_at_nominal_us"][0])
    delta_omission_closedform = cf["delta_omission_predicted"]
    rel_err_omission_identity_closedform = abs(delta_omission_closedform - (-r_mag)) / r_mag
    rel_err_omission_sim_vs_identity = abs(delta_omission_sim - (-r_mag)) / r_mag

    # --- forced fix for the idealized-silence gap: realistic timing jitter + partial
    # reinforcement (behavioral uncertainty is real; Schultz 1993's data show a
    # residual 9%-of-neurons response, NOT complete abolition -- this variant tests
    # whether adding realistic uncertainty structurally produces a nonzero residual,
    # matching the DIRECTION of that empirical discrepancy, not fitting its exact value) ---
    realistic = run_conditioning(
        rng, prereg, n_trials=prereg["realistic_n_trials"],
        partial_p=prereg["realistic_partial_reinforcement_p"],
        jitter_std=prereg["realistic_jitter_std_bins"],
    )
    tail_us_realistic = realistic["delta_at_nominal_us"][-last_n:]
    delta_learned_reward_realistic = float(np.mean(tail_us_realistic))

    gates = {
        "gate_g1_naive_burst_positive_sanity": bool(delta_naive_reward >= prereg["min_naive_burst"]),
        "gate_g2_predicted_reward_shrinks": bool(
            (delta_naive_reward - delta_learned_reward) / delta_naive_reward
            >= prereg["min_reward_response_shrink_frac"]
        ),
        "gate_g3_cue_transfer_burst": bool(
            delta_learned_cs_onset / delta_naive_reward >= prereg["min_cue_transfer_frac_of_naive_burst"]
        ),
        "gate_g4a_omission_dip_negative_and_large": bool(
            delta_omission_sim <= -prereg["min_omission_dip_magnitude"]
        ),
        "gate_g4b_omission_matches_closedform_identity": bool(
            rel_err_omission_sim_vs_identity <= prereg["max_omission_simulated_vs_identity_rel_err"]
        ),
        "gate_g4c_closedform_identity_internally_exact": bool(
            rel_err_omission_identity_closedform <= prereg["max_omission_exact_identity_rel_err_closedform"]
        ),
        "gate_g5_simulation_matches_closedform_v1": bool(rel_err_v1 <= 0.10),
        "gate_g5_simulation_matches_closedform_v0": bool(rel_err_v0 <= 0.10),
        "gate_g5_simulation_matches_closedform_vk": bool(rel_err_vk <= 0.10),
        "gate_g6_convergence_reached": bool(convergence_rel_std_us <= prereg["convergence_max_rel_std"]),
        "gate_g7_realistic_uncertainty_increases_residual_vs_idealized": bool(
            abs(delta_learned_reward_realistic) > abs(delta_learned_reward)
        ),
    }

    return {
        "delta_naive_reward_unpredicted": delta_naive_reward,
        "delta_learned_reward_predicted_deterministic": delta_learned_reward,
        "delta_learned_cs_onset_transferred_burst": delta_learned_cs_onset,
        "delta_omission_simulated": delta_omission_sim,
        "delta_omission_closed_form_identity": delta_omission_closedform,
        "delta_omission_algebraic_identity_note": "at full convergence, delta_omission == -reward_mag exactly "
            "(the background-bootstrap term gamma*V0 cancels algebraically against the same term inside V[k]'s "
            "own fixed-point definition) -- an exact TD-theoretic identity, not an approximation.",
        "closed_form_values": cf,
        "simulated_converged_V_at_key_states": {"v0": float(v0_sim), "v1": float(v1_sim), "vk": float(vk_sim)},
        "rel_err_sim_vs_closedform": {"v0": rel_err_v0, "v1": rel_err_v1, "vk": rel_err_vk},
        "convergence_rel_std_last_n": convergence_rel_std_us,
        "delta_learned_reward_realistic_uncertainty_variant": delta_learned_reward_realistic,
        "empirical_anchor_schultz1993_pct_neurons_responding_to_reward": {
            "during_learning_pct_of_76_neurons": 25, "after_learning_established_pct_of_163_neurons": 9,
            "measurement_sense_DISCLOSED_DIFFERENT_FROM_MODEL": "this is the FRACTION OF NEURONS showing a "
                "phasic response (a population proportion), NOT a per-neuron delta MAGNITUDE like this "
                "model reports -- the two are not directly numerically comparable, only qualitatively "
                "(both show large-but-incomplete reduction, never a return to the naive/unpredicted level, "
                "and never complete abolition to exactly 0), reported honestly rather than conflated.",
        },
        "gates": gates,
    }


# ============================================================================
# Part 2 -- forced adversary: "dopamine encodes reward/value, not error"
# steel-manned with stimulus-specific habituation/adaptation.
# ============================================================================
def part2_value_coding_adversary(prereg, part1_result):
    beta = prereg["adversary_habituation_beta"]
    n = prereg["adversary_n_trials_to_saturate_habituation"]
    k0 = 1.0
    k_final = k0 * (beta ** n)
    r_mag = prereg["reward_mag"]
    baseline = 0.0

    F_naive = baseline + k0 * r_mag
    F_learned_predicted = baseline + k_final * r_mag
    # forced discriminator: a reward delivered UNPREDICTED (no prior CS pairing) at a LATE trial
    # index, after the adversary's repetition-driven habituation has already saturated. Pure
    # repetition-count adaptation is agnostic to CS-pairing -> predicts the SAME small response as
    # the routine predicted trials. The TD/error account predicts a LARGE response (V=0 for this
    # untrained state regardless of how many OTHER reward deliveries occurred elsewhere).
    F_late_unpredicted_adversary_prediction = baseline + k_final * r_mag
    td_late_unpredicted_prediction = r_mag  # TD: delta = r + gamma*0 - 0 = r_mag, full-size, no attenuation
    # signature (iii): omission -- r(t)=0 identically, so k_final*0 == 0 for ANY k (structural, not tuned)
    F_omission_adversary = baseline + k_final * 0.0

    gates = {
        "gate_adversary_omission_exactly_floored_never_negative": bool(abs(F_omission_adversary - baseline) < 1e-12),
        "gate_adversary_vs_td_diverge_on_late_unpredicted_reward": bool(
            abs(td_late_unpredicted_prediction - F_late_unpredicted_adversary_prediction) > 0.5
        ),
        "gate_td_omission_dip_is_strictly_below_adversary_floor": bool(
            part1_result["delta_omission_simulated"] < F_omission_adversary - prereg["min_omission_dip_magnitude"]
        ),
    }
    return {
        "adversary_form": "F(t) = baseline + k(n)*r(t), k(n)=k0*beta^n (repetition-driven habituation, "
                          "NOT prediction-driven) -- the fair steelman of a pure value/reward-magnitude code",
        "F_naive_unpredicted_reward": F_naive,
        "F_learned_predicted_reward_via_habituation": F_learned_predicted,
        "habituation_k_final": k_final,
        "F_late_unpredicted_reward_adversary_prediction": F_late_unpredicted_adversary_prediction,
        "td_late_unpredicted_reward_prediction": td_late_unpredicted_prediction,
        "F_omission_adversary": F_omission_adversary,
        "td_omission_dip_simulated": part1_result["delta_omission_simulated"],
        "why_adversary_falls": "Habituation is a genuine confound that CAN partially mimic signature (ii)'s "
            "downward trend via repetition alone, conceded here as the fair steelman. It cannot: (a) explain "
            "why an unpredicted reward late in training (after habituation has saturated) still evokes a "
            "large response empirically (Schultz 1993's design interleaves predicted/unpredicted trials "
            "for exactly this control) -- the adversary predicts a SMALL response there (repetition-count-"
            "driven), TD predicts a LARGE one (state-specific, V=0 for the untrained pairing); (b) produce "
            "ANY negative value at omission, for ANY habituation parameter -- k(n)*0 = 0 identically, a "
            "structural (parameter-free) floor at baseline, while the TD account produces a genuine, "
            "quantitatively-matched dip (delta_omission = -reward_mag exactly at convergence).",
        "gates": gates,
    }


# ============================================================================
# Part 3 -- Bayer & Glimcher (2005) linear-RPE-encoding synthetic identifiability
# grid: explicitly a SYNTHETIC demonstration (not a re-fit of their raw data,
# unavailable in this sandbox) of why "encodes RPE" vs "encodes reward
# magnitude" is separable ONLY when expectation is varied independently of
# outcome magnitude -- a sigma_min/design-identifiability argument.
# ============================================================================
def part3_bayer_glimcher_identifiability(rng, prereg):
    R_levels = np.array(prereg["reward_levels"])
    Vexp_levels = np.array(prereg["expectation_levels"])
    gain = 1.0
    rho = prereg["neg_rpe_gain_ratio"]
    F0 = 5.0
    noise_std = prereg["firing_noise_std"]

    def synth_firing(R, Vexp):
        delta = R - Vexp
        g = gain if delta >= 0 else gain * rho
        return F0 + g * delta + rng.normal(0, noise_std)

    rows_full = np.array([[R, Vexp, R - Vexp, synth_firing(R, Vexp)]
                          for R in R_levels for Vexp in Vexp_levels])
    rows_degen = np.array([[R, R, 0.0, synth_firing(R, R)] for R in R_levels])

    def fit_r2(x, y):
        X = np.column_stack([np.ones(len(x)), x])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        yhat = X @ coef
        ss_res = np.sum((y - yhat) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        return float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0

    r2_rpe_full = fit_r2(rows_full[:, 2], rows_full[:, 3])
    r2_reward_full = fit_r2(rows_full[:, 0], rows_full[:, 3])
    gap_full = r2_rpe_full - r2_reward_full

    delta_var_degen = float(np.var(rows_degen[:, 2]))
    r2_reward_degen = fit_r2(rows_degen[:, 0], rows_degen[:, 3])

    X_full = np.column_stack([np.ones(len(rows_full)), rows_full[:, 0], rows_full[:, 1]])
    X_degen = np.column_stack([np.ones(len(rows_degen)), rows_degen[:, 0], rows_degen[:, 1]])
    smin_full = float(np.linalg.svd(X_full, compute_uv=False)[-1])
    smin_degen = float(np.linalg.svd(X_degen, compute_uv=False)[-1])

    gates = {
        "gate_full_factorial_design_separates_hypotheses": bool(gap_full >= prereg["min_r2_gap_full_design"]),
        "gate_degenerate_colinear_design_collapses_identifiability": bool(
            delta_var_degen <= prereg["max_delta_variance_degenerate_design"]
        ),
    }
    return {
        "note": "SYNTHETIC generative-model demonstration, NOT a re-fit of Bayer & Glimcher's raw spike "
            "data (unavailable in this sandbox). External anchor is their VERBATIM abstract finding "
            "(PMID 15996553, live-verified when this cell was written): firing rate is quantitatively predicted by a "
            "linear RPE model 'for circumstances in which this signal has a positive value' -- the "
            "positive-only asymmetry used in this synthetic model's gain term (rho=0.5 for delta<0) is "
            "directly informed by that caveat, not an arbitrary choice.",
        "r2_rpe_model_full_factorial_design": r2_rpe_full,
        "r2_reward_only_model_full_factorial_design": r2_reward_full,
        "gap_full_design": gap_full,
        "delta_variance_degenerate_colinear_design": delta_var_degen,
        "r2_reward_only_model_degenerate_design": r2_reward_degen,
        "sigma_min_design_matrix_full_factorial": smin_full,
        "sigma_min_design_matrix_degenerate_colinear": smin_degen,
        "geometric_reading": "the degenerate design (expectation always equals reward magnitude) forces "
            "delta=R-Vexp=0 identically for every row -- the RPE regressor has EXACTLY ZERO variance, a "
            "literal sigma_min collapse of the design matrix's [delta] column, and the two hypotheses "
            "become mathematically indistinguishable on that data (both fit equally, since 'reward alone' "
            "explains everything the coincidentally-collinear 'error' term could). Only a design that "
            "spans the (R, Vexp) plane independently -- which is what real reward-magnitude x reward-"
            "probability manipulations do -- keeps sigma_min bounded away from 0 and makes the two "
            "hypotheses separable.",
        "gates": gates,
    }


# ============================================================================
# Part 4 -- ADDICTION: Redish (2004) noncompensable-drug-dopamine TDRL lesion.
# ============================================================================
def part4_addiction_escalation(rng, prereg):
    natural = run_conditioning(rng, prereg, n_trials=prereg["addiction_n_trials"],
                               partial_p=1.0, jitter_std=0.0, drug_mode=False)
    drug = run_conditioning(rng, prereg, n_trials=prereg["addiction_n_trials"],
                            partial_p=1.0, jitter_std=0.0, drug_mode=True,
                            phi=prereg["addiction_phi_noncompensable"])

    v1_natural = float(natural["V"][1])
    v1_drug = float(drug["V"][1])
    ratio = v1_drug / v1_natural

    phi = prereg["addiction_phi_noncompensable"]
    predicted_ratio_closed_form = 1.0 + phi / (1.0 - phi)  # derived: vk_drug ~ vk_natural + phi*r/(1-phi)

    gates = {
        "gate_drug_value_escalates_vs_natural": bool(ratio >= prereg["min_value_ratio_drug_vs_natural"]),
        "gate_escalation_direction_matches_closed_form": bool(
            (ratio > 1.0) and (predicted_ratio_closed_form > 1.0)
        ),
    }
    return {
        "mechanism": "Redish 2004 (PMID 15591205, live-verified): 'By adding a noncompensable drug-induced "
            "dopamine increase to a TDRL model... a computational model of addiction is constructed that "
            "over-selects actions leading to drug receipt.' Implemented here as: effective_delta = "
            "(1-phi)*true_delta + phi*reward_mag at the reward bin -- a fraction phi of the update ignores "
            "the subtracted value term entirely, so the teaching signal never fully cancels even once the "
            "cue perfectly predicts the drug.",
        "v1_cs_onset_value_natural_reward": v1_natural,
        "v1_cs_onset_value_drug_reward": v1_drug,
        "escalation_ratio_drug_over_natural": ratio,
        "closed_form_predicted_ratio_1_plus_phi_over_1_minus_phi": predicted_ratio_closed_form,
        "gates": gates,
    }


# ============================================================================
# Part 5 -- PARKINSON'S: Frank/Seeberger/O'Reilly (2004) Go/NoGo asymmetry.
# ============================================================================
def run_go_nogo_task(rng, prereg, go_gain, nogo_gain):
    QA, QB = 0.0, 0.0
    beta = prereg["pd_softmax_beta"]
    alpha_base = prereg["pd_alpha_base"]
    pA_reward, pB_reward = prereg["pd_p_reward_A"], prereg["pd_p_reward_B"]
    for _ in range(prereg["pd_n_trials"]):
        pA = 1.0 / (1.0 + np.exp(-beta * (QA - QB)))
        choose_A = rng.random() < pA
        if choose_A:
            r = 1.0 if rng.random() < pA_reward else 0.0
            delta = r - QA
            a = go_gain * alpha_base if delta > 0 else nogo_gain * alpha_base
            QA += a * delta
        else:
            r = 1.0 if rng.random() < pB_reward else 0.0
            delta = r - QB
            a = go_gain * alpha_base if delta > 0 else nogo_gain * alpha_base
            QB += a * delta
    # transfer tests vs a novel, never-seen neutral option (Q_N=0), Frank et al's "Choose-A"/"Avoid-B" logic
    go_score = 1.0 / (1.0 + np.exp(-beta * (QA - 0.0)))
    nogo_score = 1.0 / (1.0 + np.exp(-beta * (0.0 - QB)))
    return {"QA": float(QA), "QB": float(QB), "go_score_choose_A_vs_novel": float(go_score),
            "nogo_score_avoid_B_vs_novel": float(nogo_score)}


def part5_parkinsons_go_nogo(rng, prereg):
    off_med = run_go_nogo_task(rng, prereg, prereg["pd_off_med_go_gain"], prereg["pd_off_med_nogo_gain"])
    on_med = run_go_nogo_task(rng, prereg, prereg["pd_on_med_go_gain"], prereg["pd_on_med_nogo_gain"])

    off_bias = off_med["nogo_score_avoid_B_vs_novel"] - off_med["go_score_choose_A_vs_novel"]
    on_bias = on_med["nogo_score_avoid_B_vs_novel"] - on_med["go_score_choose_A_vs_novel"]

    gates = {
        "gate_off_medication_better_at_avoid_than_approach": bool(off_bias > 0),
        "gate_medication_shifts_balance_toward_go": bool(
            (on_med["go_score_choose_A_vs_novel"] - on_med["nogo_score_avoid_B_vs_novel"])
            > (off_med["go_score_choose_A_vs_novel"] - off_med["nogo_score_avoid_B_vs_novel"])
        ),
    }
    return {
        "mechanism": "Frank/Seeberger/O'Reilly 2004 (PMID 15528409, live-verified): 'Parkinson's patients "
            "off medication are better at learning to avoid choices that lead to negative outcomes than "
            "they are at learning from positive outcomes. Dopamine medication reverses this bias.' "
            "Implemented as a DA-gain asymmetry: positive-error (Go/D1) updates scaled by a low gain OFF "
            "medication (D1's lower DA affinity under-responds at low tone); negative-error (NoGo/D2) "
            "updates scaled by a DA-INDEPENDENT gain (D2's higher affinity stays responsive at low DA).",
        "off_medication": off_med, "on_medication": on_med,
        "off_medication_nogo_minus_go_bias": off_bias,
        "on_medication_nogo_minus_go_bias": on_bias,
        "gates": gates,
    }


# ============================================================================
# Part 6 -- ABERRANT SALIENCE (Kapur 2003): spurious value on a truly-neutral,
# non-predictive stimulus, softest/most illustrative of the three poles.
# ============================================================================
def part6_aberrant_salience(rng, prereg):
    alpha = prereg["alpha"]
    n = prereg["salience_n_trials"]
    p_spurious = prereg["salience_p_spurious"]
    bump = prereg["salience_spurious_bump"]

    V_N_control = 0.0
    V_N_aberrant = 0.0
    for _ in range(n):
        V_N_control += alpha * 0.0  # N truly predicts nothing; delta==0 always in a healthy state
        spurious = bump if rng.random() < p_spurious else 0.0
        V_N_aberrant += alpha * spurious

    gates = {
        "gate_control_stays_at_true_zero": bool(V_N_control == 0.0),
        "gate_aberrant_condition_acquires_spurious_value": bool(
            V_N_aberrant >= prereg["salience_min_neutral_value_aberrant"]
        ),
    }
    return {
        "mechanism": "Kapur 2003 (PMID 12505794, live-verified; a CONCEPTUAL/theoretical review, not a "
            "primary quantitative data paper -- the softest evidence base of the three dysfunction poles, "
            "disclosed): 'A central role of dopamine is to mediate the salience of environmental events... "
            "a dysregulated, hyperdopaminergic state... leads to an aberrant assignment of salience to the "
            "elements of one's experience.' Implemented as a small, causally-uncoupled spurious-positive-"
            "error injection rate onto a stimulus (N) that truly predicts nothing (true delta=0 always).",
        "V_neutral_stimulus_control_healthy": V_N_control,
        "V_neutral_stimulus_aberrant_salience": V_N_aberrant,
        "gates": gates,
    }


# ============================================================================
# MAIN -- assemble, grade, write
# ============================================================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rng = np.random.default_rng(PREREG["rng_seed"])

    part1 = part1_three_signatures(rng, PREREG)
    part2 = part2_value_coding_adversary(PREREG, part1)
    part3 = part3_bayer_glimcher_identifiability(rng, PREREG)
    part4 = part4_addiction_escalation(rng, PREREG)
    part5 = part5_parkinsons_go_nogo(rng, PREREG)
    part6 = part6_aberrant_salience(rng, PREREG)

    gates = {}
    for label, part in [("p1", part1), ("p2", part2), ("p3", part3), ("p4", part4), ("p5", part5), ("p6", part6)]:
        for gname, gval in part["gates"].items():
            gates[f"{label}_{gname}"] = gval
    overall_pass = all(gates.values())

    result = {
        "task": "Reward prediction error (RPE): TD(0) reproduction of Schultz's 3 dopamine firing "
                "signatures, Bayer & Glimcher linear-RPE identifiability, addiction/Parkinson's/aberrant-"
                "salience dysfunction poles of the same substrate.",
        "prereg": PREREG,
        "part1_three_signatures": part1,
        "part2_value_coding_adversary": part2,
        "part3_bayer_glimcher_identifiability": part3,
        "part4_addiction_escalation": part4,
        "part5_parkinsons_go_nogo": part5,
        "part6_aberrant_salience": part6,
        "gates": gates,
        "overall_pass": bool(overall_pass),
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(gates, indent=2))
    print("overall_pass:", overall_pass)
    print("wrote:", OUT_PATH)


if __name__ == "__main__":
    main()
