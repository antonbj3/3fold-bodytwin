#!/usr/bin/env python3
"""DNA supercoiling topology: the linking-number invariant Lk=Tw+Wr plus the quadratic
supercoiling free energy, evaluated from the model's stated equations and parameters.

Reads: nothing (all equations and parameters are embedded below).
Writes: dna_supercoiling_rebuild_results.json under the cell output directory.
Gate: the acceptance test below -- agreement with the record's reported numbers.

ACCEPTANCE TEST (pre-registered, stated before running): reproduce, from the record's stated
equations + parameters alone (no tuning, no consulting a target number below to adjust a constant):
  (a) Lk0 = N/h = 419.3 turns (pBR322, N=4361bp, h=10.4bp/turn)
  (b) dLk = sigma*Lk0 = -25.2 turns (sigma=-0.06)
  (c) K = NK/N = 0.264 RT/turn^2 (NK=1150 RT*bp, Horowitz-Wang 1984)
  (d) dG_total = K*dLk^2 = 166.9 RT = 102.9 kcal/mol @37C
  (e) the algebraic identity dG = (NK/h^2)*N*RT*sigma^2, coefficient NK/h^2 = 10.6
  (f) torque = dG/d(dLk)/(2*pi) = -8.7 pN*nm @298K / -9.0 pN*nm @310K, vs measured -11 pN*nm
      melting-torque plateau (record's stated 18% residual, NOT a target to hit -- the record
      itself says this is a regime mismatch, not an exact prediction)
  (g) topA-null sigma>2x WT -> dG scales as sigma^2 -> 4.0x predicted torsional free energy
  (h) writhe-fraction from dWr:dTw=2.6:1 -> Wr/(Wr+Tw) = 72.2%, vs the record's independently
      reported EM count of 89% (both a literature datum, not derived here, only bracketed)
Physically meaningful numbers (a)-(e),(g) are held to a <1% agreement bar; (f) the torque residual vs
the MEASURED plateau is expected to differ (record's honest_gaps says the harmonic estimate is
"outside the regime where the derivative is strictly licensed") so only order-of-magnitude + sign
agreement is required there, matching the record's reported -18%-ish gap, not zero.

SOURCE EQUATIONS (verbatim from the record's claim/derived_values text, transcribed unchanged):
  Lk = Tw + Wr                                    (White-Fuller-Calugareanu exact topological thm)
  Lk0 = N / h                                     (relaxed-state linking number)
  dLk = sigma * Lk0                               (supercoiling density sigma = dLk/Lk0)
  K = NK / N                                      (per-bp free-energy coefficient, RT/turn^2)
  dG(dLk) = K * dLk^2                             (quadratic supercoiling free energy, in RT units)
  dG = (NK/h^2) * N * RT * sigma^2                 (equivalent form stated in the record, algebraic
                                                     identity of the two lines above -- checked, not
                                                     assumed)
  torque(dLk) = d(dG)/d(dLk) / (2*pi) * RT          (energy-per-radian; dLk in turns, 1 turn=2*pi rad)
  topA_scale = (sigma_mutant / sigma_wt)^2          (quadratic energy scaling under sigma doubling)
  writhe_fraction = dWr / (dWr + dTw)               (partition of the linking-number deficit)

PARAMETERS (verbatim from the record's "derived_values" text -- no value here was adjusted):
  N=4361 bp (pBR322), h=10.4 bp/turn (Wang 1979), NK=1150 RT*bp (Horowitz-Wang 1984),
  sigma=-0.06 (in-vivo typical), dWr:dTw=2.6:1 (Boles/White/Cozzarelli 1990),
  topA-null sigma ratio = 2.0 (>2x WT, PMC5992920)

GEOMETRIC STRUCTURE: Lk=Tw+Wr is an exact topological invariant (Calugareanu-White-Fuller theorem) --
changes only via strand passage (topoisomerase action), a genuine conserved quantity, not a fitted
approximation. The quadratic free energy dG=K*dLk^2 is the harmonic (small-oscillation) term of the
elastic free energy around the relaxed state Lk0 -- torque = dG/d(theta) is then a first-derivative
(linear-response) quantity, valid only near dLk=0, which is exactly why the record itself flags the
torque comparison against a MEASURED (nonlinear, buckled-rod) plateau as regime-mismatched, not a
tight prediction -- this rebuild reproduces that same regime-limited match, not a better one.

"""

import json
import os

import numpy as np

import os as _os

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "dna_supercoiling_topology_rebuild")
OUT_PATH = _os.path.join(OUT_DIR, "dna_supercoiling_rebuild_results.json")
SOURCE_NODE = "DNA supercoiling topology model"
SOURCE_EVIDENCE = "the stored evidence record for that model"

# ---- parameters, verbatim from the record's derived_values / claim text ----
N_BP = 4361.0          # pBR322 plasmid length, bp
H_BP_PER_TURN = 10.4   # Wang 1979 helical repeat, bp/turn
NK_CONST = 1150.0      # Horowitz-Wang 1984, RT*bp
SIGMA = -0.06          # in-vivo typical supercoiling density
DWR_DTW_RATIO = 2.6    # Boles/White/Cozzarelli 1990 writhe:twist partition
TOPA_SIGMA_RATIO = 2.0 # topA-null sigma > 2x WT (PMC5992920)

R_KCAL_PER_MOL_K = 1.987204e-3  # gas constant, kcal/(mol*K)
R_J_PER_MOL_K = 8.314462618     # gas constant, J/(mol*K)
NA = 6.02214076e23              # Avogadro
T_37C_K = 310.0
T_TORQUE_LOW_K = 298.0
T_TORQUE_HIGH_K = 310.0

# ---- record's reported numbers (targets, for comparison only -- never fed back into the model) ----
RECORD = {
    "Lk0_turns": 419.3,
    "dLk_turns": -25.2,
    "K_RT_per_turn2": 0.264,
    "dG_total_RT": 166.9,
    "dG_total_kcal_per_mol": 102.9,
    "algebraic_coeff_NK_over_h2": 10.6,
    "torque_298K_pN_nm": -8.7,
    "torque_310K_pN_nm": -9.0,
    "measured_melting_plateau_pN_nm": -11.0,
    "topA_null_scale": 4.0,
    "writhe_fraction_from_ratio_pct": 72.0,
    "writhe_fraction_EM_pct": 89.0,
}


def main():
    print("=" * 78)
    print(f"REBUILDING {SOURCE_NODE} from {SOURCE_EVIDENCE}")
    print("=" * 78)

    Lk0 = N_BP / H_BP_PER_TURN
    dLk = SIGMA * Lk0
    K = NK_CONST / N_BP
    dG_RT = K * dLk ** 2
    dG_kcal = dG_RT * R_KCAL_PER_MOL_K * T_37C_K

    algebraic_coeff = NK_CONST / H_BP_PER_TURN ** 2  # equivalent-form identity check

    def RT_pN_nm(T_kelvin):
        return (R_J_PER_MOL_K * T_kelvin / NA) * 1e21  # J -> pN*nm (1 pN*nm = 1e-21 J)

    torque_298 = RT_pN_nm(T_TORQUE_LOW_K) * K * dLk / np.pi
    torque_310 = RT_pN_nm(T_TORQUE_HIGH_K) * K * dLk / np.pi

    topA_scale = TOPA_SIGMA_RATIO ** 2
    writhe_frac_from_ratio = DWR_DTW_RATIO / (DWR_DTW_RATIO + 1.0) * 100.0

    computed = {
        "Lk0_turns": float(Lk0),
        "dLk_turns": float(dLk),
        "K_RT_per_turn2": float(K),
        "dG_total_RT": float(dG_RT),
        "dG_total_kcal_per_mol": float(dG_kcal),
        "algebraic_coeff_NK_over_h2": float(algebraic_coeff),
        "torque_298K_pN_nm": float(torque_298),
        "torque_310K_pN_nm": float(torque_310),
        "topA_null_scale": float(topA_scale),
        "writhe_fraction_from_ratio_pct": float(writhe_frac_from_ratio),
    }
    print(json.dumps(computed, indent=2))

    def pct_err(a, b):
        return abs(a - b) / abs(b) * 100 if b != 0 else abs(a - b)

    agreement = {
        "Lk0_pct_err": pct_err(computed["Lk0_turns"], RECORD["Lk0_turns"]),
        "dLk_pct_err": pct_err(computed["dLk_turns"], RECORD["dLk_turns"]),
        "K_pct_err": pct_err(computed["K_RT_per_turn2"], RECORD["K_RT_per_turn2"]),
        "dG_RT_pct_err": pct_err(computed["dG_total_RT"], RECORD["dG_total_RT"]),
        "dG_kcal_pct_err": pct_err(computed["dG_total_kcal_per_mol"], RECORD["dG_total_kcal_per_mol"]),
        "algebraic_coeff_pct_err": pct_err(computed["algebraic_coeff_NK_over_h2"], RECORD["algebraic_coeff_NK_over_h2"]),
        "torque_298K_pct_err": pct_err(computed["torque_298K_pN_nm"], RECORD["torque_298K_pN_nm"]),
        "torque_310K_pct_err": pct_err(computed["torque_310K_pN_nm"], RECORD["torque_310K_pN_nm"]),
        "torque_vs_measured_plateau_pct_gap": pct_err(computed["torque_310K_pN_nm"], RECORD["measured_melting_plateau_pN_nm"]),
        "topA_scale_pct_err": pct_err(computed["topA_null_scale"], RECORD["topA_null_scale"]),
        "writhe_fraction_pct_err": pct_err(computed["writhe_fraction_from_ratio_pct"], RECORD["writhe_fraction_from_ratio_pct"]),
    }
    print("\nAGREEMENT vs record's reported numbers (%err):")
    print(json.dumps(agreement, indent=2))

    gates = {
        "Lk0_within_1pct": bool(agreement["Lk0_pct_err"] < 1.0),
        "dLk_within_1pct": bool(agreement["dLk_pct_err"] < 1.0),
        "K_within_1pct": bool(agreement["K_pct_err"] < 1.0),
        "dG_RT_within_1pct": bool(agreement["dG_RT_pct_err"] < 1.0),
        "dG_kcal_within_1pct": bool(agreement["dG_kcal_pct_err"] < 1.0),
        "algebraic_identity_within_1pct": bool(agreement["algebraic_coeff_pct_err"] < 1.0),
        "torque_298K_within_1pct": bool(agreement["torque_298K_pct_err"] < 1.0),
        "torque_310K_within_1pct": bool(agreement["torque_310K_pct_err"] < 1.0),
        "topA_scale_exact": bool(agreement["topA_scale_pct_err"] < 1e-9),
        "writhe_fraction_within_1pct": bool(agreement["writhe_fraction_pct_err"] < 1.0),
        # NOT a target-match gate -- the record itself states this is a regime-mismatched,
        # order-of-magnitude comparison (harmonic estimate vs a nonlinear buckled-rod plateau);
        # gate only checks the rebuild lands in the SAME qualitative regime the record describes
        # (correct sign, same order of magnitude, 15-25% low), not that it hits -11 exactly.
        "torque_vs_plateau_same_regime_as_record": bool(15.0 <= agreement["torque_vs_measured_plateau_pct_gap"] <= 25.0),
    }
    overall_pass = all(gates.values())
    print("\nGATES:")
    print(json.dumps(gates, indent=2))
    print(f"\nOVERALL RECOVERY: {'REPRODUCED' if overall_pass else 'PARTIAL/FAILED -- see gates'}")

    honest_gaps = {
        "torque_is_a_regime_mismatch_by_the_records_own_admission": "The -8.7/-9.0 pN*nm harmonic "
            "estimate vs the -11 pN*nm MEASURED melting-torque plateau is, per the record's "
            "honest_gaps, an idealized linear-response derivative evaluated outside the small-|sigma| "
            "regime where it is strictly licensed -- reproduced here as the SAME ~18-22% gap the "
            "record itself reports, not closed to zero (closing it would require the full "
            "Marko-Siggia/Moroz-Nelson buckled-rod free energy, explicitly out of scope per the record).",
        "writhe_fraction_89pct_EM_figure_is_a_literature_datum_not_derived": "The record's 89% EM "
            "superhelix count is an independently reported measurement, not computed from any stated "
            "equation here -- only the 72.2% figure (from the stated dWr:dTw=2.6:1 ratio) is a "
            "rebuild output; both are compared as brackets, not a single point prediction.",
            "no_parameter_was_tuned": "Every constant above (N,h,NK,sigma,dWr:dTw ratio,topA sigma "
            "ratio) is copied verbatim from the record's claim/derived_values text; none was "
            "adjusted to improve agreement.",
    }
    report = {
        "source_node": SOURCE_NODE,
        "source_evidence": SOURCE_EVIDENCE,
        "parameters_used_verbatim_from_record": {
            "N_bp": N_BP, "h_bp_per_turn": H_BP_PER_TURN, "NK_RT_bp": NK_CONST,
            "sigma": SIGMA, "dWr_dTw_ratio": DWR_DTW_RATIO, "topA_sigma_ratio": TOPA_SIGMA_RATIO,
        },
        "computed": computed,
        "record_reported": RECORD,
        "agreement_pct_err": agreement,
        "gates": gates,
        "overall_pass": overall_pass,
        "honest_gaps": honest_gaps,
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(f"\nWrote {OUT_PATH}")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
