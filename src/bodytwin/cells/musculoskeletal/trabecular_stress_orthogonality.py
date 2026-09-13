"""Trabecular alignment vs single-load-case principal-stress-trajectory theory: does the naive
theory (trabeculae align with the principal stress directions of ONE load instant) survive the
measured femoral-neck trabecular crossing angle, given that ANY single load case's principal
directions are forced to be exactly 90 degrees apart (a spectral-theorem fact about symmetric 2x2
stress tensors, independent of beam geometry)?

Scope: this cell tests only the two facts that do not depend on a full curved-beam/FE machinery
(both model-independent); a reduced axial+shear-only dominant-angle estimate is reported as
INFORMATIONAL ONLY and is not gated, because it omits the bending-moment term.

Model-independent facts tested:
  1. ORTHOGONALITY: for any symmetric 2D stress tensor [[sigma, tau],[tau, 0]] (axial + in-plane
     shear, zero transverse normal stress -- the standard beam cross-section state), the two
     principal directions are exactly 90 degrees apart, regardless of magnitudes or load case.
  2. REFUTATION ARITHMETIC: Skedros & Baucom 2007 measured a femoral-neck crossing angle of ~70
     degrees (non-orthogonal); |90 - 70| = 20 deg > the 15 deg falsifier threshold.

Reads: nothing. Writes: trabecular_stress_orthogonality.json.

Gates (fixed before computing):
  G1 orthogonality: for 200 random (sigma, tau) draws (sigma in [-200,200] MPa, tau in
     [-100,100] MPa, tau != 0), the eigenvector angle difference is 90.0 +/- 1e-6 degrees in
     every draw.
  G2 refutation arithmetic: |90 - 70| > 15 -> REFUTED.
  G3 void floor: break the SYMMETRY of the stress tensor (tau_12 != tau_21, an
     equilibrium-violating matrix) for the same 200 draws -> the eigenvector angle difference must
     deviate from 90 degrees in at least 90% of draws.
  G4 (informational, not gated) reduced axial+shear-only dominant angle for both disclosed sign
     branches, reported alongside Hadad et al 2025's measured 31.8-32.3 deg, explicitly
     scope-limited (omits bending).
Verdict WEAKENED for the orthogonality sub-claim iff G1-G3 pass.
"""
import json
import os

import numpy as np

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT = _os.path.join(OUT_ROOT, "trabecular_stress_orthogonality", "trabecular_stress_orthogonality.json")

MEASURED_CROSSING_DEG = 70.0  # Skedros & Baucom 2007, femoral neck
ORTHOGONAL_DEG = 90.0
FALSIFIER_THRESHOLD_DEG = 15.0
HADAD_DOMINANT_DEG_RANGE = (31.8, 32.3)

NSA_DEG = 125.8  # Hadad2025-measured neck-shaft angle
THETA_R_DEG = 18.4  # Pauwels central hip-joint-reaction angle from vertical

RNG_SEED = 20260728


def eigvec_angle_diff_deg(sigma, tau, tau21=None):
    """Angle between the two eigenvectors of [[sigma, tau], [tau21_or_tau, 0]], in degrees."""
    if tau21 is None:
        tau21 = tau
    M = np.array([[sigma, tau], [tau21, 0.0]])
    if not np.allclose(M, M.T):
        # non-symmetric: eigenvectors of a general 2x2 real matrix may be complex/non-orthogonal;
        # use np.linalg.eig and take the real part's angle if real eigenvectors exist.
        w, v = np.linalg.eig(M)
        if np.iscomplexobj(w) and np.any(np.abs(w.imag) > 1e-9):
            return None  # complex eigenvalues: no real principal directions at all
        v = v.real
    else:
        w, v = np.linalg.eigh(M)
    v1 = v[:, 0] / np.linalg.norm(v[:, 0])
    v2 = v[:, 1] / np.linalg.norm(v[:, 1])
    cos_ang = np.clip(abs(np.dot(v1, v2)), -1.0, 1.0)
    ang = np.degrees(np.arccos(cos_ang))
    # angle BETWEEN the two directions (undirected lines), take the acute-complement-safe value
    return ang


def reduced_dominant_angle_deg(phi_deg):
    """theta_p = 0.5*arctan(2*tan(phi)) -- beam-area-independent principal-angle formula
    for a pure axial(N)+shear(V)-only stress state with V/N = tan(phi)."""
    phi = np.radians(phi_deg)
    return 0.5 * np.degrees(np.arctan(2.0 * np.tan(phi)))


def main():
    rng = np.random.default_rng(RNG_SEED)
    n_draws = 200
    sigmas = rng.uniform(-200.0, 200.0, n_draws)
    taus = rng.uniform(-100.0, 100.0, n_draws)
    taus = np.where(np.abs(taus) < 1.0, 1.0, taus)  # avoid tau==0 degenerate draws

    # G1: symmetric tensor -> exactly 90deg every time
    sym_diffs = [eigvec_angle_diff_deg(s, t) for s, t in zip(sigmas, taus)]
    g1 = all(abs(d - ORTHOGONAL_DEG) <= 1e-6 for d in sym_diffs if d is not None)
    gates = {"G1_orthogonality_exact_90deg_all_draws": bool(g1)}

    # G2: refutation arithmetic
    gap = abs(ORTHOGONAL_DEG - MEASURED_CROSSING_DEG)
    gates["G2_refutation_gap_exceeds_15deg"] = bool(gap > FALSIFIER_THRESHOLD_DEG)

    # G3: void floor -- asymmetric tensor (broken equilibrium: tau21 = -tau, sign-flipped)
    asym_diffs = []
    n_deviant = 0
    for s, t in zip(sigmas, taus):
        d = eigvec_angle_diff_deg(s, t, tau21=-t)
        if d is None or abs(d - ORTHOGONAL_DEG) > 1e-3:
            n_deviant += 1
    frac_deviant = n_deviant / n_draws
    gates["G3_void_floor_asymmetric_deviates_ge_90pct"] = bool(frac_deviant >= 0.90)

    # G4: informational reduced-model dominant angle, both sign branches (NOT gated)
    branch_A_phi = abs((180.0 - NSA_DEG) - THETA_R_DEG)  # neck-tilt and R-tilt same side
    branch_B_phi = (180.0 - NSA_DEG) + THETA_R_DEG        # opposite sides
    theta_p_A = reduced_dominant_angle_deg(branch_A_phi)
    theta_p_B = reduced_dominant_angle_deg(branch_B_phi)

    gates = {k: bool(v) for k, v in gates.items()}
    verdict = "WEAKENED" if all(gates.values()) else "DISAGREE"

    result = {
        "node_id": "MODEL-TRABECULAR-WOLFF-STRESS-ALIGNMENT",
        "orthogonality_gap_deg": gap,
        "fraction_deviant_under_asymmetry_void_floor": frac_deviant,
        "reduced_model_informational_only": {
            "branch_A_phi_deg": branch_A_phi,
            "branch_B_phi_deg": branch_B_phi,
            "theta_p_A_deg": theta_p_A,
            "theta_p_B_deg": theta_p_B,
            "hadad_measured_range_deg": list(HADAD_DOMINANT_DEG_RANGE),
            "note": (
                "NOT gated: this reduced axial+shear-only model omits the bending-moment term "
                "the source's curved-beam+FE model includes (per instruction to check before "
                "concluding a disagreement). Reported for transparency only."
            ),
        },
        "gates": gates,
        "claimed_by_narrative": {
            "verdict": "WEAKENED",
            "orthogonality_forced_90deg": True,
            "measured_crossing_deg": 70.0,
            "gap_exceeds_15deg_falsifier": True,
        },
        "verdict": verdict,
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict,
                       "orthogonality_gap_deg": gap,
                       "fraction_deviant_void_floor": frac_deviant}, indent=2))
    return result


if __name__ == "__main__":
    main()
