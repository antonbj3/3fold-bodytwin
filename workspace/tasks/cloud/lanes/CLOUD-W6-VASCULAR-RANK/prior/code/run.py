"""Design optimisation + held-out Monte Carlo coverage. Usage: python3 run.py [n_sims]"""
import json, sys, time
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm
import model as m

SEED = 20260924
N_SIMS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
LT0 = np.log(m.THETA_NOM)
Z90 = norm.ppf(0.95)
E = np.array([])


def logdet_reg(F):
    return np.linalg.slogdet(F + 1e-9 * np.trace(F) * np.eye(len(F)))[1]


def greedy(tp, ti, n_add, pool):
    """Greedy D-optimal addition (with replacement) of n_add points from pool ('p' or 'i')."""
    tp, ti = list(tp), list(ti)
    for _ in range(n_add):
        best = None
        for t in m.CAND:
            a, b = (tp + [t], ti) if pool == "p" else (tp, ti + [t])
            v = logdet_reg(m.fisher(LT0, np.array(a), np.array(b)))
            if best is None or v > best[0]:
                best = (v, t)
        (tp if pool == "p" else ti).append(best[1])
    return np.array(sorted(tp)), np.array(sorted(ti))


def summarise(tp, ti, lt=LT0):
    F = m.fisher(lt, tp, ti)
    r, s = m.fim_rank(F)
    return {"tp": tp.tolist(), "ti": ti.tolist(), "rank": r,
            "sv_ratio_min": float(s[-1] / s[0]), "logdet_reg": float(logdet_reg(F))}


def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(c - h), float(c + h)]


def coverage(tp, ti, n, rng, cv=0.25):
    cov = np.zeros(4); ok = 0; fails = 0; widths = []
    for _ in range(n):
        lt_true = LT0 + rng.normal(0, np.log(1 + cv), 4)      # held-out truth
        y = m.observe(lt_true, tp, ti)
        y = y + np.concatenate([rng.normal(0, m.SD_P, len(tp)), rng.normal(0, m.SD_I, len(ti))])
        w = np.concatenate([np.full(len(tp), 1 / m.SD_P), np.full(len(ti), 1 / m.SD_I)])
        try:
            res = least_squares(lambda lt: (m.observe(lt, tp, ti) - y) * w, LT0,
                                bounds=(LT0 - 4, LT0 + 4), x_scale=1.0, xtol=1e-10, ftol=1e-10)
            cv_ = np.linalg.inv(m.fisher(res.x, tp, ti))
            se = np.sqrt(np.diag(cv_))
            if not np.all(np.isfinite(se)):
                raise np.linalg.LinAlgError
        except (np.linalg.LinAlgError, ValueError):
            fails += 1; continue
        ok += 1
        cov += np.abs(res.x - lt_true) <= Z90 * se
        widths.append(2 * Z90 * se)
    return {"n_requested": n, "n_fitted": ok, "n_failed": fails,
            "coverage": dict(zip(m.NAMES, (cov / max(ok, 1)).round(4).tolist())),
            "coverage_wilson95": {k: wilson(c, ok) for k, c in zip(m.NAMES, cov)},
            "median_ci_width_log": dict(zip(m.NAMES, np.median(widths, 0).round(4).tolist())) if widths else None}


def main():
    t0 = time.time(); out = {"seed": SEED, "n_sims": N_SIMS, "theta_nominal": dict(zip(m.NAMES, m.THETA_NOM.tolist())),
                             "noise_sd": {"plasma": m.SD_P, "interstitial": m.SD_I}, "rank_tol": m.RANK_TOL}
    # conservation check
    Cp, Ci, y = m.simulate(m.THETA_NOM, m.CAND)
    out["conservation_max_drift"] = {"volume_L": float(np.ptp(y[0] + y[1])), "mass_g": float(np.ptp(y[2] + y[3]))}
    # designs
    P6, _ = greedy([], [], 6, "p")
    P8, _ = greedy(P6, [], 2, "p")
    _, I2 = greedy(P6, [], 2, "i")
    _, I1 = greedy(P6, [], 1, "i")
    d = {"P6_plasma_only": summarise(P6, E), "P8_plasma_only_equal_budget": summarise(P8, E),
         "P48_plasma_all_candidates": summarise(m.CAND, E),
         "P6_plus_I1": summarise(P6, I1), "P6_plus_I2": summarise(P6, I2)}
    # minimal plasma count with added interstitial points reaching rank 4
    minimal = None
    for n in range(1, 7):
        tp, _ = greedy([], [], n, "p"); _, ti = greedy(tp, [], 1, "i")
        s = summarise(tp, ti)
        if s["rank"] == 4:
            minimal = {"n_plasma": n, "n_interstitial": 1, **s}; break
    d["minimal_full_rank_design"] = minimal
    # counterexamples
    d["counterexample_I2_early_only"] = summarise(P6, np.array([5.0, 10.0]))
    lt_kf0 = LT0.copy(); lt_kf0[0] = np.log(1e-6)
    d["counterexample_nearly_no_fluid_shift_kf1e-6_P6_plus_I2"] = summarise(P6, I2, lt_kf0)
    out["designs"] = d
    # null direction of plasma-only FIM
    F = m.fisher(LT0, m.CAND, E); w, V = np.linalg.eigh(F)
    out["plasma_only_weakest_direction"] = dict(zip(m.NAMES, V[:, 0].round(4).tolist()))
    # Monte Carlo coverage (held-out truths)
    rng = np.random.default_rng(SEED)
    out["mc_P6_plus_I2"] = coverage(P6, I2, N_SIMS, rng)
    out["mc_P8_plasma_only"] = coverage(P8, E, min(N_SIMS, 300), rng)
    out["mc_counterexample_high_noise_x5_P6_plus_I2"] = None
    sp, si = m.SD_P, m.SD_I; m.SD_P, m.SD_I = 5 * sp, 5 * si
    out["mc_counterexample_high_noise_x5_P6_plus_I2"] = coverage(P6, I2, min(N_SIMS, 300), rng)
    m.SD_P, m.SD_I = sp, si
    cov = out["mc_P6_plus_I2"]["coverage"].values()
    rank_up = d["P6_plus_I2"]["rank"] > d["P6_plasma_only"]["rank"] and d["P6_plus_I2"]["rank"] > d["P48_plasma_all_candidates"]["rank"]
    cov_ok = all(0.85 <= c <= 0.95 for c in cov) and out["mc_P6_plus_I2"]["n_fitted"] >= 1000
    out["criterion"] = {"rank_increase": bool(rank_up), "coverage_85_95_all_params_ge1000_fits": bool(cov_ok),
                        "synthetic_PASS": bool(rank_up and cov_ok), "real_world_optimal_schedule": "UNKNOWN"}
    out["runtime_s"] = round(time.time() - t0, 1)
    json.dump(out, open("../results.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
