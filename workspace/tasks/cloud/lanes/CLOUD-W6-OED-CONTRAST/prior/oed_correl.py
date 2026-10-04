#!/usr/bin/env python3
"""CLOUD-W3-OED-CORREL: correlated measurement error in a minimal measurement protocol.

Fully SYNTHETIC linear-Gaussian study (Python standard library only, no external data).

Parameters theta (mm): [pelvis width, thigh length, shank length] of one subject.
Candidate measurements (8): 5 landmark-palpation measures (L*) and 3 imaging measures (I*).
Measurement error of candidate j in modality g for subject s:
    e_j = mu_j + b_{g,s} + eps_{j,s}
where b_{g,s} is a per-subject SYSTEMATIC bias shared by all measurements of modality g
(soft-tissue offset for palpation, calibration offset for imaging) -> within-modality
correlation rho_g = var(b_g) / (var(b_g) + var(eps)).

Pipeline per scenario:
 1. Calibration cohort (N_CAL synthetic subjects with known truth, all 8 measured):
    estimate mean bias mu_j and error covariance
      - covariance-aware: structured estimate (shared + independent per modality)
      - naive: diagonal of the same sample covariance (independence assumption)
 2. Optimal design: choose M of 8 candidates minimising posterior variance of the target
    c'theta (c-optimality), exhaustively, under each assumed covariance.
 3. Held-out cohort (N_TEST new subjects, independent seed): posterior mean and 90%
    interval of the target on each model's own selected design; empirical coverage.
 4. Exact (analytic) actual coverage of each interval given the TRUE covariance.
 5. Replicate calibrations to quantify variability of held-out coverage.

Usage:  python3 oed_correl.py            (writes results.json)
"""
import itertools
import json
import math
import platform
import random
import sys

Z90 = 1.6448536269514722  # standard normal 0.95 quantile
N_CAL = 100
N_TEST = 5000
N_REP = 200          # replicate calibrations
N_REP_TEST = 1000    # held-out subjects per replicate
M = 3                # measurement budget (minimal protocol)
SEED = 20260924

# ---------------------------------------------------------------- linear algebra
def matmul(A, B):
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(r, c)) for c in Bt] for r in A]

def T(A):
    return [list(r) for r in zip(*A)]

def inv(A):
    n = len(A)
    M_ = [list(map(float, r)) + [1.0 if i == j else 0.0 for j in range(n)] for i, r in enumerate(A)]
    for i in range(n):
        p = max(range(i, n), key=lambda k: abs(M_[k][i]))
        M_[i], M_[p] = M_[p], M_[i]
        piv = M_[i][i]
        if abs(piv) < 1e-14:
            raise ValueError("singular")
        M_[i] = [v / piv for v in M_[i]]
        for k in range(n):
            if k != i:
                f = M_[k][i]
                if f:
                    M_[k] = [a - f * b for a, b in zip(M_[k], M_[i])]
    return [r[n:] for r in M_]

def add(A, B):
    return [[a + b for a, b in zip(r, s)] for r, s in zip(A, B)]

def sub(A, B):
    return [[a - b for a, b in zip(r, s)] for r, s in zip(A, B)]

def chol(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = A[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = math.sqrt(s) if i == j else s / L[j][j]
    return L

def mv(A, x):
    return [sum(a * b for a, b in zip(r, x)) for r in A]

def quad(c, A):
    return sum(ci * v for ci, v in zip(c, mv(A, c)))

def Phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def wilson(k, n, z=1.959963984540054):
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [centre - half, centre + half]

# ---------------------------------------------------------------- synthetic world
NAMES = ["L1_pelvis", "L2_thigh", "L3_shank", "L2b_thigh_repeat", "L3b_shank_repeat",
         "I1_pelvis", "I2_thigh", "I3_shank"]
H_ALL = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 0], [0, 0, 1],
         [1, 0, 0], [0, 1, 0], [0, 0, 1]]
MOD = ["L", "L", "L", "L", "L", "I", "I", "I"]
MU0 = [240.0, 430.0, 410.0]
SD0 = [18.0, 24.0, 22.0]
CORR0 = [[1, .3, .3], [.3, 1, .6], [.3, .6, 1]]
SIG0 = [[CORR0[i][j] * SD0[i] * SD0[j] for j in range(3)] for i in range(3)]
L0 = chol(SIG0)
SD_TOT = {"L": 7.0, "I": 9.0}      # total error SD per measurement (mm)
MEAN_BIAS = {"L": 8.0, "I": -3.0}  # systematic mean offset (mm), removed via calibration

def true_R(rho, mult_imaging_scale=None):
    """True error covariance for additive shared bias with within-modality correlation rho[g].
    If mult_imaging_scale is set, imaging bias is multiplicative (scale error s ~ N(0, sd^2))
    and the covariance returned is the population (marginal over theta) covariance."""
    R = [[0.0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            gi, gj = MOD[i], MOD[j]
            if gi != gj:
                continue
            tot = SD_TOT[gi] ** 2
            if gi == "I" and mult_imaging_scale is not None:
                # e = s * (h'theta) + eps ; indep part keeps (1-rho) share of SD_TOT
                ind = (1 - rho[gi]) * tot
                hi, hj = H_ALL[i], H_ALL[j]
                Eij = quad_pair(hi, hj)
                R[i][j] = mult_imaging_scale ** 2 * Eij + (ind if i == j else 0.0)
            else:
                R[i][j] = tot * (1.0 if i == j else rho[gi])
    return R

def quad_pair(hi, hj):
    # E[(hi'theta)(hj'theta)] under the prior
    m = sum(a * b for a, b in zip(hi, MU0)) * sum(a * b for a, b in zip(hj, MU0))
    return m + sum(hi[a] * SIG0[a][b] * hj[b] for a in range(3) for b in range(3))

def draw_subject(rng, rho, mult_imaging_scale=None):
    z = [rng.gauss(0, 1) for _ in range(3)]
    th = [MU0[i] + sum(L0[i][k] * z[k] for k in range(3)) for i in range(3)]
    b = {g: rng.gauss(0, math.sqrt(rho[g]) * SD_TOT[g]) for g in "LI"}
    s = rng.gauss(0, mult_imaging_scale) if mult_imaging_scale is not None else None
    y = []
    for j in range(8):
        g = MOD[j]
        ht = sum(a * b_ for a, b_ in zip(H_ALL[j], th))
        eps = rng.gauss(0, math.sqrt(1 - rho[g]) * SD_TOT[g])
        bias = (s * ht) if (g == "I" and s is not None) else b[g]
        y.append(ht + MEAN_BIAS[g] + bias + eps)
    return th, y

# ---------------------------------------------------------------- calibration
def calibrate(rng, rho, mult):
    res = []
    for _ in range(N_CAL):
        th, y = draw_subject(rng, rho, mult)
        res.append([y[j] - sum(a * b for a, b in zip(H_ALL[j], th)) for j in range(8)])
    mu = [sum(r[j] for r in res) / N_CAL for j in range(8)]
    S = [[sum((r[i] - mu[i]) * (r[j] - mu[j]) for r in res) / (N_CAL - 1) for j in range(8)]
         for i in range(8)]
    R_aware = [[0.0] * 8 for _ in range(8)]
    for g in "LI":
        idx = [j for j in range(8) if MOD[j] == g]
        off = [S[i][j] for i in idx for j in idx if i != j]
        shared = max(0.0, sum(off) / len(off))
        for i in idx:
            for j in idx:
                R_aware[i][j] = shared if i != j else max(S[i][i] - shared, 1e-6) + shared
    R_naive = [[S[i][j] if i == j else 0.0 for j in range(8)] for i in range(8)]
    return mu, R_aware, R_naive

# ---------------------------------------------------------------- design + inference
SIG0_INV = inv(SIG0)

def sub_R(R, idx):
    return [[R[i][j] for j in idx] for i in idx]

def posterior(idx, R):
    H = [H_ALL[j] for j in idx]
    Ri = inv(sub_R(R, idx))
    P = inv(add(SIG0_INV, matmul(matmul(T(H), Ri), H)))
    K = matmul(matmul(P, T(H)), Ri)
    return P, K

def select(R, c):
    best = None
    for idx in itertools.combinations(range(8), M):
        v = quad(c, posterior(list(idx), R)[0])
        if best is None or v < best[1] - 1e-12:
            best = (list(idx), v)
    return best

def actual_var(idx, K, R_true, c):
    """Exact error variance of c'(posterior mean) - c'theta given TRUE error covariance
    (additive model; for multiplicative scenario the population covariance)."""
    H = [H_ALL[j] for j in idx]
    I3 = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    A = sub(I3, matmul(K, H))
    V = add(matmul(matmul(A, SIG0), T(A)), matmul(matmul(K, sub_R(R_true, idx)), T(K)))
    return quad(c, V)

def run_test(rng, n, rho, mult, designs, c):
    """designs: dict name -> (idx, mu_hat, R_assumed). Returns hits, widths per design."""
    prepared = {}
    for name, (idx, mu, R) in designs.items():
        P, K = posterior(idx, R)
        # posterior mean m = P Sig0^-1 mu0 + K (y - mu)
        a = mv(matmul(P, SIG0_INV), MU0)
        prepared[name] = (idx, mu, K, a, math.sqrt(quad(c, P)))
    hits = {k: 0 for k in designs}
    sqerr = {k: 0.0 for k in designs}
    for _ in range(n):
        th, y = draw_subject(rng, rho, mult)
        tgt = sum(ci * t for ci, t in zip(c, th))
        for name, (idx, mu, K, a, sd) in prepared.items():
            yy = [y[j] - mu[j] for j in idx]
            m = [a[i] + sum(K[i][k] * yy[k] for k in range(len(idx))) for i in range(3)]
            est = sum(ci * v for ci, v in zip(c, m))
            hits[name] += abs(est - tgt) <= Z90 * sd
            sqerr[name] += (est - tgt) ** 2
    return {k: {"coverage": hits[k] / n, "coverage_wilson95": wilson(hits[k], n),
                "claimed_halfwidth_mm": Z90 * prepared[k][4],
                "rmse_mm": math.sqrt(sqerr[k] / n)} for k in designs}

def scenario(name, rho, c, mult=None, seed_off=0):
    rng_cal = random.Random(SEED + 1000 * seed_off + 1)
    rng_test = random.Random(SEED + 1000 * seed_off + 2)
    Rt = true_R(rho, mult)
    mu, Ra, Rn = calibrate(rng_cal, rho, mult)
    sel = {"aware": select(Ra, c), "naive": select(Rn, c), "oracle": select(Rt, c)}
    designs = {
        "aware": (sel["aware"][0], mu, Ra),
        "naive": (sel["naive"][0], mu, Rn),
        # decomposition: covariance-aware inference on the naive-selected design
        "aware_on_naive_design": (sel["naive"][0], mu, Ra),
    }
    test = run_test(rng_test, N_TEST, rho, mult, designs, c)
    out = {"rho": rho, "target_c": c, "imaging_multiplicative_scale_sd": mult,
           "selected": {k: [NAMES[j] for j in v[0]] for k, v in sel.items()},
           "claimed_target_sd_mm": {k: math.sqrt(v[1]) for k, v in sel.items()},
           "heldout": test, "analytic": {}}
    for k, (idx, _, R) in designs.items():
        P, K = posterior(idx, R)
        va, vc = actual_var(idx, K, Rt, c), quad(c, P)
        out["analytic"][k] = {"actual_sd_mm": math.sqrt(va), "claimed_sd_mm": math.sqrt(vc),
                              "coverage_exact_given_calibrated_R":
                                  2 * Phi(Z90 * math.sqrt(vc / va)) - 1}
    # replicate calibrations: variability of held-out coverage
    rng_rep = random.Random(SEED + 1000 * seed_off + 3)
    covs = {"aware": [], "naive": []}
    same = 0
    for _ in range(N_REP):
        mu_r, Ra_r, Rn_r = calibrate(rng_rep, rho, mult)
        da, dn = select(Ra_r, c)[0], select(Rn_r, c)[0]
        same += da == dn
        r = run_test(rng_rep, N_REP_TEST, rho, mult,
                     {"aware": (da, mu_r, Ra_r), "naive": (dn, mu_r, Rn_r)}, c)
        for k in covs:
            covs[k].append(r[k]["coverage"])
    def summ(v):
        v = sorted(v)
        return {"median": v[len(v) // 2], "p2.5": v[int(0.025 * len(v))],
                "p97.5": v[int(0.975 * len(v)) - 1],
                "frac_in_85_95": sum(0.85 <= x <= 0.95 for x in v) / len(v)}
    out["replicates"] = {"n_rep": N_REP, "n_test_each": N_REP_TEST,
                         "frac_same_design": same / N_REP,
                         **{k: summ(v) for k, v in covs.items()}}
    print(f"[{name}] sel aware={out['selected']['aware']} naive={out['selected']['naive']} "
          f"cov aware={test['aware']['coverage']:.3f} naive={test['naive']['coverage']:.3f}",
          file=sys.stderr)
    return out

def main():
    SUM, DIFF = [0, 1, 1], [0, 1, -1]
    res = {"meta": {"python": platform.python_version(), "seed": SEED, "N_CAL": N_CAL,
                    "N_TEST": N_TEST, "N_REP": N_REP, "N_REP_TEST": N_REP_TEST, "budget_M": M,
                    "candidates": NAMES, "H": H_ALL, "prior_mean": MU0, "prior_cov": SIG0,
                    "sd_total": SD_TOT, "mean_bias": MEAN_BIAS, "data": "synthetic only"},
           "scenarios": {}}
    sc = res["scenarios"]
    sc["S1_high_corr_sum"] = scenario("S1", {"L": 0.8, "I": 0.8}, SUM, seed_off=1)
    sc["S2_moderate_corr_sum"] = scenario("S2", {"L": 0.5, "I": 0.5}, SUM, seed_off=2)
    sc["S3_no_corr_sum"] = scenario("S3", {"L": 0.0, "I": 0.0}, SUM, seed_off=3)
    sc["S4_high_corr_contrast"] = scenario("S4", {"L": 0.8, "I": 0.8}, DIFF, seed_off=4)
    sc["S5_high_corr_sum_mult_imaging"] = scenario("S5", {"L": 0.8, "I": 0.8}, SUM,
                                                   mult=0.02, seed_off=5)
    s1 = sc["S1_high_corr_sum"]["heldout"]
    ca, cn = s1["aware"]["coverage"], s1["naive"]["coverage"]
    res["prereg"] = {
        "criterion": "PASS if covariance-aware 90% intervals cover 85-95% and naive intervals "
                     "fail (coverage outside 85-95%) in high-correlation case across >=1000 trials",
        "case": "S1_high_corr_sum", "n_trials": N_TEST,
        "aware_coverage": ca, "naive_coverage": cn,
        "aware_in_85_95": 0.85 <= ca <= 0.95, "naive_outside_85_95": not (0.85 <= cn <= 0.95),
        "verdict": "PASS (synthetic)" if (0.85 <= ca <= 0.95 and not 0.85 <= cn <= 0.95)
                   else "FAIL"}
    with open("results.json", "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res["prereg"], indent=1))

if __name__ == "__main__":
    main()
