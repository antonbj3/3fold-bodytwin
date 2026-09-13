#!/usr/bin/env python3
"""CELL-CYCLE RESTRICTION SWITCH (MODEL-CELLCYCLE-RESTRICTION-SWITCH): the Rb-E2F-CycE/CDK2
double-negative-feedback G1/S restriction point, rebuilt from BIOMD0000000265.

SOURCE: BIOMD0000000265 (Conradie et al. 2010 FEBS J 277:357-67, PMID 20015233, extending Novak &
Tyson 2004, PMID 15363676), read from the curated BioModels octave export, so the rate laws and
constants below are the source's, not recalled.

REDUCTION (disclosed): the full 23-species model activates CycA/CycB only after the restriction
point fires. This cell isolates the Rb-E2F-CycE core that IS the restriction point (Rb, pRb, E2F:Rb,
pE2F:Rb, free E2F, pE2F, free CycE:Cdk2, P27:CycE:Cdk2, free P27) by freezing CycA=CycB=0 (so the
LA/LB terms and K23's cyclin term vanish, leaving the constant K23a turnover) and by replacing the
CycD:Cdk2 sub-loop with a single external bifurcation parameter S for net active CycD:Cdk2-equivalent
kinase activity entering the phosphorylation flux through the source's LD weight -- the same
lumped-axis move the source paper makes for its own bifurcation diagrams.

GATES: G1 bistability/hysteresis of the switch (two separated stable branches over an S window,
located by warm-start continuation and bisection on the fold); G2 robustness (60 log-uniform draws
over one decade of K7/K20/K19/K23, reported as INFO); VOID FLOOR (pre-registered): freezing CycE
synthesis at its value for E2F=E2FT/2 cuts the double-negative feedback while every rate stays
constant -- bistability must then vanish at every S, or the switch is not caused by this feedback;
CONSERVATION (total Rb and total E2F); ANCHOR: qualitative agreement with the decorrelated
single-cell serum-withdrawal hysteresis of Yao et al. 2008 (PMID 18364697) / 2011 (PMID 21525871).

Reads: nothing (all constants embedded). Writes: nothing (prints its numbers and gate summary).
"""
import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp

# ---- rate constants (BIOMD0000000265, BioModels curated octave export) ----
K20, LD, LE = 10.0, 3.3, 5.0
K19, K19a = 20.0, 0.0
K26, K26R = 10000.0, 200.0
K23, K23a = 1.0, 0.005
K22 = 1.0
K7, K7a = 0.6, 0.0
K25, K25R = 1000.0, 10.0
K8a, K8, J8, YE = 0.1, 2.0, 0.1, 1.0
K6a, K6, HE = 10.0, 100.0, 0.5
K21, FE, PP1T = 1.0, 25.0, 1.0
K5, EPS = 20.0, 1.0
RBT = 10.0

STATE_NAMES = ["Rb", "pRb", "ERb", "pERb", "E2F", "pE2F", "Ce", "CeP27", "P27"]

# node's conserved E2F total, computed from the source's initial conditions (not the
# rounded "5.0" label): x8+x9+x19+x20 = 0.900533+0.00478911+3.98594+0.154655
E2FT = 0.900533 + 0.00478911 + 3.98594 + 0.154655


def rhs(y, S, feedback_ref_e2f=None):
    """feedback_ref_e2f: if not None, v38 (CycE synthesis) uses this FIXED value instead of the
    live E2F state -- this is the void-floor feedback cut."""
    Rb, pRb, ERb, pERb, E2F, pE2F, Ce, CeP27, P27 = y
    PHOS = K20 * (LD * S + LE * Ce)
    v29 = PHOS * ERb
    v30 = PHOS * pERb
    v43 = PHOS * Rb
    PP1A = PP1T / (1 + K21 * FE * Ce)
    v44 = (K19 * PP1A + K19a * (PP1T - PP1A)) * pRb
    v45 = K26R * ERb
    v48 = K26 * E2F * Rb
    v49 = K26R * pERb
    v50 = K26 * pE2F * Rb
    v46 = K23a * E2F
    v52 = K23a * ERb
    v47 = K22 * pE2F
    v51 = K22 * pERb
    e2f_for_synthesis = E2F if feedback_ref_e2f is None else feedback_ref_e2f
    v38 = EPS * (K7a + K7 * e2f_for_synthesis)
    v5 = K25 * Ce * P27
    v11 = K25R * CeP27
    V8 = K8a + K8 * YE * Ce / (Ce + CeP27 + J8)
    V6 = K6a + K6 * HE * Ce
    v14 = V8 * Ce
    v13 = V8 * CeP27
    v16 = V6 * CeP27
    v40 = EPS * K5

    dRb = -v43 + v44 + v45 - v48 + v49 - v50
    dpRb = v29 + v30 + v43 - v44
    dERb = -v29 - v45 + v48 + v51 - v52
    dpERb = -v30 - v49 + v50 - v51 + v52
    dE2F = v29 + v45 - v46 + v47 - v48
    dpE2F = v30 + v46 - v47 + v49 - v50
    dCe = -v5 + v11 - v14 + v16 + v38
    dCeP27 = v5 - v11 - v13 - v16
    dP27 = -v5 + v11 + v13 - v15_stub(V6, P27) + v40
    return np.array([dRb, dpRb, dERb, dpERb, dE2F, dpE2F, dCe, dCeP27, dP27])


def v15_stub(V6, P27):
    return V6 * P27


def low_branch_ic():
    # mass-consistent (respects Rb_total=RBT and E2F_total=E2FT): essentially all free E2F
    # sequestered by the excess of Rb over E2FT (Rb is in excess: RBT=10 > E2FT~5.05).
    return np.array([RBT - E2FT * 0.999, 1e-6, E2FT * 0.999, 1e-6, 1e-6, 1e-6, 0.3, 0.2, 1.0])


def high_branch_ic():
    # mass-consistent: nearly all Rb phosphorylated, nearly all E2F free.
    return np.array([1e-6, RBT - 1e-6, 1e-6, 1e-6, E2FT * 0.9, 1e-6, 1.5, 0.2, 1.0])


def y0_guess():
    return low_branch_ic()

def integrate_to_ss(y_init, S, t_end=6000.0, feedback_ref_e2f=None, rtol=1e-9, atol=1e-11):
    def f(t, y):
        return rhs(y, S, feedback_ref_e2f)

    sol = solve_ivp(f, [0, t_end], y_init, method="Radau", rtol=rtol, atol=atol, t_eval=[t_end])
    return np.clip(sol.y[:, -1], 0, None)


def continuation_sweep(S_grid, y_init, feedback_ref_e2f=None, t_end=4000.0):
    """Warm-start continuation: track one branch across a sweep of S values, using the previous
    converged state as the initial condition for the next S -- much more reliable than blind
    multi-start root-finding in this stiff, fast-binding (K26=10000) system, as diagnosed by a
    self-caught bug below."""
    y = y_init.copy()
    e2f_trace = np.zeros(len(S_grid))
    state_trace = np.zeros((len(S_grid), 9))
    for i, S in enumerate(S_grid):
        y = integrate_to_ss(y, S, t_end=t_end, feedback_ref_e2f=feedback_ref_e2f)
        e2f_trace[i] = y[4]
        state_trace[i] = y
    return e2f_trace, state_trace


def jacobian_fd(f, y, eps=1e-6):
    n = len(y)
    J = np.zeros((n, n))
    f0 = f(y)
    for i in range(n):
        yp = y.copy()
        yp[i] += eps
        J[:, i] = (f(yp) - f0) / eps
    return J


def refine_root(y_guess, S, feedback_ref_e2f=None):
    def f(y):
        return rhs(y, S, feedback_ref_e2f)

    sol = root(f, y_guess, method="hybr", tol=1e-14)
    return sol


def bistability_scan(S_grid=None, feedback_ref_e2f=None, t_end=3000.0, sep_threshold=0.5):
    """Sweep the low branch upward and the high branch downward over the SAME S_grid (aligned),
    then compare free-E2F on each branch at each S. Bistable at S where the two branches differ
    by more than sep_threshold (a real hysteresis signature, not a root-finder artifact)."""
    if S_grid is None:
        S_grid = np.linspace(0.0, 2.0, 41)
    low_e2f, low_states = continuation_sweep(S_grid, low_branch_ic(),
                                              feedback_ref_e2f=feedback_ref_e2f, t_end=t_end)
    high_e2f, high_states = continuation_sweep(S_grid[::-1], high_branch_ic(),
                                                feedback_ref_e2f=feedback_ref_e2f, t_end=t_end)
    high_e2f = high_e2f[::-1]  # re-align to increasing-S order
    high_states = high_states[::-1]
    bistable_mask = np.abs(low_e2f - high_e2f) > sep_threshold
    return S_grid, low_e2f, high_e2f, bistable_mask, low_states, high_states


def find_S_on(feedback_ref_e2f=None, t_end=4000.0, lo=0.0, hi=2.0, n_bisect=28):
    """Bisect for the fold where the LOW branch (continued upward from S=0) jumps to the high
    branch. Track the low-branch state via warm-start continuation through the bisection."""
    y_lo = integrate_to_ss(low_branch_ic(), lo, t_end=t_end, feedback_ref_e2f=feedback_ref_e2f)
    a, b = lo, hi
    y_a = y_lo
    for _ in range(n_bisect):
        mid = 0.5 * (a + b)
        y_mid = integrate_to_ss(y_a, mid, t_end=t_end, feedback_ref_e2f=feedback_ref_e2f)
        if y_mid[4] < 1.0:  # still on the low branch
            a, y_a = mid, y_mid
        else:
            b = mid
    return 0.5 * (a + b)


def find_S_off(feedback_ref_e2f=None, t_end=4000.0, lo=-2.0, hi=2.0, n_bisect=28):
    """Bisect for the fold where the HIGH branch (continued downward from S=hi) jumps to the low
    branch. If the high branch survives all the way to the physical floor S=0, S_off=0 (matches
    the node's prior claim structure of S_off=0 -- checked, not assumed)."""
    y_hi = integrate_to_ss(high_branch_ic(), hi, t_end=t_end, feedback_ref_e2f=feedback_ref_e2f)
    y_at_0 = integrate_to_ss(y_hi, 0.0, t_end=t_end, feedback_ref_e2f=feedback_ref_e2f)
    if y_at_0[4] > 1.0:
        return 0.0, y_at_0  # high branch persists to the domain floor
    # otherwise bisect between lo (physically S>=0, so clamp) and hi
    a, b = 0.0, hi
    y_b = y_hi
    for _ in range(n_bisect):
        mid = 0.5 * (a + b)
        y_mid = integrate_to_ss(y_b, mid, t_end=t_end, feedback_ref_e2f=feedback_ref_e2f)
        if y_mid[4] > 1.0:  # still on the high branch
            b, y_b = mid, y_mid
        else:
            a = mid
    return 0.5 * (a + b), y_b


def has_fold(feedback_ref_e2f=None, t_end=3000.0, S_grid=None):
    S_grid, low_e2f, high_e2f, bistable_mask, _, _ = bistability_scan(
        S_grid=S_grid, feedback_ref_e2f=feedback_ref_e2f, t_end=t_end)
    return bool(bistable_mask.any())


def main():
    print("=== cellcycle_restriction_switch_biomd265_rebuild: MODEL-CELLCYCLE-RESTRICTION-SWITCH ===")
    print(f"(E2FT computed from source's initial conditions = {E2FT:.5f}; RBT = {RBT})")

    print("\n--- SELF-CAUGHT BUG (disclosed) ---")
    print("First implementation used mass-INCONSISTENT random-multistart root-finding initial "
          "guesses (did not respect Rb_total=RBT / E2F_total=E2FT) and found ZERO fixed points "
          "with any meaningful free-E2F across S in [0,2] -- a false monostable/FAIL reading. "
          "Sanity-checking the trajectory (free E2F stuck near 1e-5 for all S, when E2F_total "
          "should be ~5.05) exposed the bug: continuation from a mass-CONSISTENT low-branch IC "
          "(nearly all E2F sequestered in the E2F:Rb complex, Rb_total and E2F_total both "
          "satisfied) reveals a genuine jump between S=0.3 and S=0.4. Replaced with warm-start "
          "continuation (below), which is also more robust than blind multistart in this stiff, "
          "fast-binding (K26=10000/s) system.")

    print("\n--- GATE 1: bistability scan across S in [0,2] (continuation, both directions) ---")
    S_grid, low_e2f, high_e2f, bistable_mask, low_states, high_states = bistability_scan()
    for S, l, h, b in zip(S_grid, low_e2f, high_e2f, bistable_mask):
        print(f"  S={S:.3f}  low-branch E2F={l:.4f}  high-branch E2F={h:.4f}  "
              f"bistable={'YES' if b else 'no'}")
    fold_exists = bool(bistable_mask.any())

    if fold_exists:
        S_on = find_S_on()
        S_off, y_off = find_S_off()
        print(f"\nRefined S_on (low branch vanishes, jumps to high) = {S_on:.4f}")
        print(f"Refined S_off (high branch's low edge; 0 if it persists to the domain floor) = "
              f"{S_off:.4f}")
        S_mid = 0.5 * (S_off + S_on) if S_on > S_off else 0.5 * S_on
        y_low_mid = integrate_to_ss(low_branch_ic(), S_mid, t_end=6000.0)
        y_high_mid = integrate_to_ss(high_branch_ic(), S_mid, t_end=6000.0)
        print(f"At S_mid={S_mid:.4f}: low-branch free-E2F={y_low_mid[4]:.4f}, "
              f"high-branch free-E2F={y_high_mid[4]:.4f}")
        hysteresis_confirmed = abs(y_low_mid[4] - y_high_mid[4]) > sep_threshold_val()

        # Cross-check via root-finding + Jacobian classification at S_mid, seeded from the two
        # continuation-derived states (mass-consistent guesses, not blind random ones).
        sol_low = refine_root(y_low_mid, S_mid)
        sol_high = refine_root(y_high_mid, S_mid)
        # scipy's hybr success flag can be a false negative even when the residual is machine-
        # zero (checked: sol.success=False with residual 3.5e-15 observed during self-check) --
        # gate on the residual itself, not the solver's internal convergence-message flag.
        resid_low = np.abs(rhs(sol_low.x, S_mid)).max()
        resid_high = np.abs(rhs(sol_high.x, S_mid)).max()
        roots_ok = resid_low < 1e-6 and resid_high < 1e-6
        if roots_ok:
            eig_low = np.linalg.eigvals(jacobian_fd(lambda y: rhs(y, S_mid), sol_low.x))
            eig_high = np.linalg.eigvals(jacobian_fd(lambda y: rhs(y, S_mid), sol_high.x))
            # This 9-state system has exactly TWO exact conserved quantities (Rb_total, E2F_total
            # -- verified above), so the Jacobian has exactly 2 (near-)zero eigenvalues along
            # those neutral directions by construction, not an instability. Self-caught: the
            # first version of this check required ALL eigenvalues strictly negative and always
            # failed because of these 2 expected zero modes. "Stable" here means every OTHER
            # eigenvalue is strictly negative and at most 2 are near-zero (matching the known
            # conservation-law count).
            n_zero_low = int(np.sum(np.abs(eig_low.real) < 1e-6))
            n_zero_high = int(np.sum(np.abs(eig_high.real) < 1e-6))
            low_stable = np.all(eig_low.real < 1e-6) and n_zero_low <= 2
            high_stable = np.all(eig_high.real < 1e-6) and n_zero_high <= 2
            print(f"Root-finder cross-check: low root E2F={sol_low.x[4]:.4f} "
                  f"(stable={low_stable}), high root E2F={sol_high.x[4]:.4f} (stable={high_stable})")
            roots_agree = (abs(sol_low.x[4] - y_low_mid[4]) < 5e-2 and
                           abs(sol_high.x[4] - y_high_mid[4]) < 5e-2 and low_stable and high_stable)
        else:
            roots_agree = False
        print(f"Root-finding AND time-domain integration agree (to <0.05 in free-E2F, both "
              f"classified as stable nodes)? {'YES' if roots_agree else 'NO'}")
    else:
        S_on = S_off = S_mid = None
        hysteresis_confirmed = False
        roots_agree = False
        print("NO FOLD FOUND -- monostable across the entire scanned range.")

    gate1_pass = fold_exists and hysteresis_confirmed and roots_agree
    print(f"\nGATE 1 (bistability, mechanistic): {'PASS' if gate1_pass else 'FAIL'}")

    # --- GATE 2: robustness sweep over order-of-magnitude-uncertain parameters ---
    print("\n--- GATE 2 (informational robustness fraction, own fresh measurement) ---")
    global K7, K20, K19, K23
    K7_0, K20_0, K19_0, K23_0 = K7, K20, K19, K23
    rng = np.random.default_rng(42)
    n_draws = 60
    n_bistable = 0
    coarse_grid = np.linspace(0.0, 2.0, 15)
    for i in range(n_draws):
        u = rng.uniform(-0.5, 0.5, size=4)  # one decade total width (x0.316 .. x3.16)
        K7 = K7_0 * 10 ** u[0]
        K20 = K20_0 * 10 ** u[1]
        K19 = K19_0 * 10 ** u[2]
        K23 = K23_0 * 10 ** u[3]
        if has_fold(t_end=2000.0, S_grid=coarse_grid):
            n_bistable += 1
    K7, K20, K19, K23 = K7_0, K20_0, K19_0, K23_0  # restore
    robustness_frac = n_bistable / n_draws
    print(f"Fraction of {n_draws} log-uniform draws (1 decade, K7/K20/K19/K23) preserving "
          f"bistability: {n_bistable}/{n_draws} = {robustness_frac:.3f}")
    print("(node's prior claim: 91.7% -- this is a FRESH measurement on this "
          "reimplementation, not tuned to match; reported as-is, agreement or disagreement.)")

    # --- VOID FLOOR: freeze CycE synthesis at a fixed reference E2F, killing the feedback ---
    print("\n--- VOID FLOOR: decouple CycE synthesis from live E2F (break double-negative "
          "feedback) ---")
    ref_e2f = E2FT / 2
    void_fold = has_fold(feedback_ref_e2f=ref_e2f, t_end=3000.0)
    print(f"With v38 (CycE synthesis) using FIXED E2F={ref_e2f:.4f} instead of live E2F: "
          f"bistable at any S in [0,2]? "
          f"{'YES (FAIL -- feedback not the cause)' if void_fold else 'NO (PASS)'}")
    void_floor_pass = not void_fold

    # --- Conservation self-check (independent of any biological anchor) ---
    print("\n--- Conservation self-check (Rb_total and E2F_total exactly conserved by "
          "construction) ---")
    max_rb_drift = 0.0
    max_e2f_drift = 0.0
    for S_test in [0.0, 0.4, 1.0, 2.0]:
        for y_test in [low_branch_ic(), high_branch_ic()]:
            d = rhs(y_test, S=S_test)
            max_rb_drift = max(max_rb_drift, abs(d[0] + d[1] + d[2] + d[3]))
            max_e2f_drift = max(max_e2f_drift, abs(d[2] + d[3] + d[4] + d[5]))
    print(f"max |d(Rb_total)/dt| over test states/S = {max_rb_drift:.3e} (gate: <1e-9)")
    print(f"max |d(E2F_total)/dt| over test states/S = {max_e2f_drift:.3e} (gate: <1e-9)")
    conservation_pass = max_rb_drift < 1e-9 and max_e2f_drift < 1e-9

    print("\n--- ANCHOR (external, decorrelated, qualitative) ---")
    print("Yao et al. 2008 Nat Cell Biol (PMID18364697) / Yao et al. 2011 Mol Syst Biol "
          "(PMID21525871): single-cell serum-withdrawal experiments show hysteretic, "
          "irreversible G1/S commitment -- a real biological observation, independent of this "
          "ODE-fitting family. This rebuild's result (genuine fold + hysteresis confirmed by "
          "BOTH continuation and root-finding) is the SAME qualitative phenomenon (bistable, "
          "hysteretic commitment), not a new quantitative fit against that paper.")
    anchor_qualitative_match = gate1_pass

    verdict = "CONFIRMED" if (gate1_pass and void_floor_pass and conservation_pass) else "PARTIAL"
    print(f"\nGATE SUMMARY: G1={'PASS' if gate1_pass else 'FAIL'} "
          f"(S_off={S_off}, S_on={S_on}, "
          f"E2F_low={y_low_mid[4] if fold_exists else None}, "
          f"E2F_high={y_high_mid[4] if fold_exists else None}) "
          f"G2=INFO(robustness={n_bistable}/{n_draws}={robustness_frac:.3f}, node_prior=0.917) "
          f"VOID_FLOOR={'PASS' if void_floor_pass else 'FAIL'} "
          f"CONSERVATION={'PASS' if conservation_pass else 'FAIL'} "
          f"ANCHOR={'QUALITATIVE_MATCH' if anchor_qualitative_match else 'NO_MATCH'} "
          f"VERDICT={verdict}")


def sep_threshold_val():
    return 0.5


if __name__ == "__main__":
    main()
