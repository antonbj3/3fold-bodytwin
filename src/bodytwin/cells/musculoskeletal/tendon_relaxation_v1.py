"""Data-anchored tendon stress-relaxation analysis (normalized).

Scope and honesty boundary
--------------------------
This module reuses a parsed copy of the UniPD human foot-tendon series, which contains
force in gram-force (gf) plus either time (stress-relaxation sheets) or actuator
position (failure sheets). The stress-relaxation data have NO displacement, NO
strain protocol and NO specimen geometry (area, gauge length). Therefore only the
NORMALIZED relaxation shape

    r(x) = F(x) / F_peak ,   x = t - t_peak

is identifiable. Absolute material parameters (elastic modulus, viscosity,
relaxation time as a physical constant of the tissue) are NOT identifiable from
these data and the module returns ABSTAIN whenever the geometry/protocol
prerequisites are absent instead of reporting a number.

Design
------
* Units are explicit: gf <-> N via the exact standard-gravity definition; any
  absolute force reporting converts to N. A unit-swap guard refuses a declared
  frame whose magnitude is inconsistent with gf.
* Models: constant (M0), single exponential (M1), stretched exponential (M2),
  power law (M3). Per-hold offset and amplitude (c, A) are profiled out, so only
  the shared shape parameter(s) are fitted across a calibration set.
* Validation is GROUPED: calibration and validation never share a specimen or a
  donor. The grouping is asserted, not assumed.
* Identifiability is reported as parameter interval width + bound-hit fraction +
  held-out margin. The tau/beta correlation is recorded as a diagnostic but is
  explicitly NOT the identifiability gate (correlation alone is not a universal
  gate).
* Unknown load protocols and unknown geometry stay declared as unmeasured
  constraints in every result envelope.

It performs no writes and no implicit I/O; the caller
loads data explicitly. Run as a script it prints a hermetic synthetic check (no
data files are read).

Attribution: the hold detector and the shared-shape fitter are derived from the UniPD
relaxation and held-out protocol analysis scripts of the 2026-09 BodyTwin review (not part
of this repository).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

try:  # reuse the installed framework contract when the repo src is importable
    from bodytwin.framework.coupling_evidence_v1 import CouplingEvidence
except Exception:  # pragma: no cover - standalone fallback, never silently accepts
    CouplingEvidence = None

__all__ = [
    "GF_TO_N", "gf_to_newton", "newton_to_gf", "assert_declared_force_frame",
    "TAU_REF_S", "TAU_REL_TOL", "LIN2_REL_TOL",
    "m0_const", "m1_exp", "m2_stretch", "m3_power",
    "basis", "lin2", "fit_shared_shape", "eval_heldout",
    "detect_holds_in_series", "load_unipd_holds",
    "make_folds", "assert_no_leakage", "grouped_validation",
    "shape_param_bootstrap_ci", "shape_param_correlation", "hit_bounds",
    "shared_at_bound",
    "practically_identifiable", "assess_geometry_prerequisites", "ABSTAIN", "PASS",
    "FAIL",
]

# --------------------------------------------------------------------- units ----
# 1 gram-force is, by definition, the weight of 1 gram under standard gravity.
GF_TO_N = 9.80665e-3

# Fixed numerical-reference tolerances, declared BEFORE any run. Unknown
# biological tolerance stays unknown and is never inferred from these.
TAU_REF_S = 13.7          # analytic single-exponential reference
TAU_REL_TOL = 1e-3        # relative tolerance on recovered reference tau
LIN2_REL_TOL = 1e-9       # relative tolerance of closed-form 2-param LSQ


def _as_float_array(values: Any, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        raise ValueError(f"{name} must be non-empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def gf_to_newton(force_gf: Any) -> np.ndarray:
    """Convert gram-force to newton (exact factor 9.80665e-3 N/gf)."""
    return _as_float_array(force_gf, "force_gf") * GF_TO_N


def newton_to_gf(force_n: Any) -> np.ndarray:
    """Convert newton to gram-force."""
    return _as_float_array(force_n, "force_n") / GF_TO_N


def assert_declared_force_frame(force: Any, declared_unit: str) -> str:
    """Guard against a gf/N frame swap for absolute force reporting.

    The UniPD values are recorded in gf with a median non-zero step of about
    1.26 gf. If a caller declares ``gf`` but the numeric magnitude is a
    plausible newton-scale absolute load, or vice versa, refuse instead of
    silently rescaling. This is a failure-mode guard for reported absolute
    force, not a physics claim about the dataset.
    """
    force = _as_float_array(force, "force")
    unit = str(declared_unit).strip()
    peak = float(np.max(np.abs(force)))
    if unit == "gf":
        # a relaxation peak of many kN recorded in gf would be implausible; a
        # newton-scale load mislabelled as gf shows up as a tiny peak (< 10 gf)
        if peak < 10.0:
            raise ValueError(
                "declared frame 'gf' but peak magnitude < 10 gf; likely a "
                "newton-valued array mislabelled as gram-force")
    elif unit == "N":
        if peak > 2000.0:
            raise ValueError(
                "declared frame 'N' but peak magnitude > 2000 N; likely a "
                "gram-force-valued array mislabelled as newton")
    else:
        raise ValueError("declared_unit must be 'gf' or 'N'")
    return unit


# -------------------------------------------------------------------- models ----
def m0_const(x: np.ndarray, c: float) -> np.ndarray:
    return np.full_like(np.asarray(x, float), c)


def m1_exp(x: np.ndarray, c: float, A: float, tau: float) -> np.ndarray:
    return c + A * np.exp(-np.asarray(x, float) / tau)


def m2_stretch(x: np.ndarray, c: float, A: float, tau: float, beta: float) -> np.ndarray:
    return c + A * np.exp(-np.power(np.asarray(x, float) / tau, beta))


def m3_power(x: np.ndarray, c: float, A: float, tau: float, beta: float) -> np.ndarray:
    return c + A * np.power(1.0 + np.asarray(x, float) / tau, -beta)


MODEL_FORMS = ("M0", "M1", "M2", "M3")
K_SHAPE = {"M0": 0, "M1": 1, "M2": 2, "M3": 2}
BASIS_MODELS = ("M1", "M2", "M3")


def basis(name: str, x: np.ndarray, tau: float, beta: float | None) -> np.ndarray:
    x = np.asarray(x, float)
    if name == "M1":
        return np.exp(-x / tau)
    if name == "M2":
        return np.exp(-np.power(x / tau, beta))
    if name == "M3":
        return np.power(1.0 + x / tau, -beta)
    raise ValueError(f"no basis for model {name!r}")


def shape_of(name: str, tau: float, beta: float | None) -> dict[str, float | None]:
    return {"tau": float(tau), "beta": (float(beta) if name in ("M2", "M3") else None)}


def lin2(g: np.ndarray, y: np.ndarray) -> tuple[tuple[float, float], float]:
    """Closed-form y ~ c*1 + A*g, returns ((c, A), rss). Reference-checked
    against numpy least squares in the self-test."""
    g = np.asarray(g, float)
    y = np.asarray(y, float)
    n = g.size
    Sg, Sgg = float(g.sum()), float(np.dot(g, g))
    Sy, Sgy, Syy = float(y.sum()), float(np.dot(g, y)), float(np.dot(y, y))
    den = n * Sgg - Sg * Sg
    if den <= 1e-12 * max(1.0, n * Sgg):
        c = Sy / n
        return (c, 0.0), max(Syy - Sy * Sy / n, 0.0)
    c = (Sgg * Sy - Sg * Sgy) / den
    A = (n * Sgy - Sg * Sy) / den
    return (c, A), max(Syy - (c * Sy + A * Sgy), 0.0)


def total_rss(name: str, shape: Mapping[str, float], curves: Sequence[tuple]) -> float:
    tau, beta = shape["tau"], shape.get("beta")
    tot = 0.0
    for x, y in curves:
        g = basis(name, x, tau, beta)
        _, rss = lin2(g, y)
        tot += rss
    return tot


def fit_shared_shape(name: str, curves: Sequence[tuple]) -> tuple[dict, float]:
    """Shared shape parameter(s) minimising total RSS with per-curve (c,A)
    profiled out. Deterministic coarse grid + local Nelder-Mead refinement."""
    from scipy.optimize import minimize

    if len(curves) == 0:
        raise ValueError("fit_shared_shape needs at least one curve")
    if name == "M0":
        tot = 0.0
        for x, y in curves:
            r = np.asarray(y, float) - 1.0
            tot += float(np.dot(r, r))
        return {}, tot

    if name == "M1":
        taus = np.logspace(np.log10(0.5), np.log10(5000.0), 32)
        scored = sorted(((total_rss(name, {"tau": float(t)}, curves), float(t))
                         for t in taus), key=lambda z: z[0])
        best, best_tau = scored[0]
        for _, t0 in scored[:3]:
            def obj(z, t0=t0):
                tau = float(np.exp(z[0]))
                if tau < 0.1 or tau > 1e5:
                    return 1e18
                return total_rss(name, {"tau": tau}, curves)
            r = minimize(obj, [np.log(t0)], method="Nelder-Mead",
                         options={"xatol": 1e-5, "fatol": 1e-7, "maxiter": 4000})
            if r.fun < best:
                best, best_tau = float(r.fun), float(np.exp(r.x[0]))
        return {"tau": best_tau, "beta": None}, best

    taus = np.logspace(np.log10(0.5), np.log10(5000.0), 24)
    betas = np.linspace(0.1, 1.9, 11)
    cand = [(total_rss(name, {"tau": float(t), "beta": float(b)}, curves), float(t), float(b))
            for t in taus for b in betas]
    cand.sort(key=lambda z: z[0])
    best = cand[0][0]
    best_shape = {"tau": cand[0][1], "beta": cand[0][2]}
    for _, t0, b0 in cand[:4]:
        def obj(z, t0=t0, b0=b0):
            tau, beta = float(np.exp(z[0])), float(z[1])
            if tau < 0.1 or tau > 1e5 or beta < 0.05 or beta > 2.0:
                return 1e18
            return total_rss(name, {"tau": tau, "beta": beta}, curves)
        r = minimize(obj, [np.log(t0), b0], method="Nelder-Mead",
                     options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 6000})
        if r.fun < best:
            best = float(r.fun)
            best_shape = {"tau": float(np.exp(r.x[0])), "beta": float(r.x[1])}
    return best_shape, best


def eval_heldout(name: str, shape: Mapping[str, float], holds: Sequence[Mapping]) -> dict:
    """Held-out evaluation on normalized curves. Per-hold (c,A) are solved on the
    held-out hold itself (nuisance, not shared). Absolute RMSE is returned in
    both gf and mN (units explicit)."""
    per, rss_n, rss_gf, nn = [], 0.0, 0.0, 0
    for h in holds:
        x, y = np.asarray(h["x"], float), np.asarray(h["y_norm"], float)
        peak = float(h["peak_force_gf"])
        if name == "M0":
            c, A, pred = 1.0, 0.0, np.ones_like(y)
        else:
            g = basis(name, x, shape["tau"], shape.get("beta"))
            (c, A), _ = lin2(g, y)
            pred = c + A * g
        res = y - pred
        rn = float(np.dot(res, res))
        rg = float(np.dot(res * peak, res * peak))
        n = y.size
        per.append({"specimen": h.get("specimen"), "donor": h.get("donor"),
                    "step": h.get("step"), "peak_force_gf": peak,
                    "n_points": int(n),
                    "c": float(c), "A": float(A),
                    "rmse_rel": float(np.sqrt(rn / n)),
                    "rmse_gf": float(np.sqrt(rg / n)),
                    "rmse_mN": float(np.sqrt(rg / n) * GF_TO_N * 1e3)})
        rss_n += rn
        rss_gf += rg
        nn += n
    return {"n_holds": len(per), "n_points": nn,
            "pooled_rel_rmse": float(np.sqrt(rss_n / nn)) if nn else None,
            "pooled_abs_rmse_gf": float(np.sqrt(rss_gf / nn)) if nn else None,
            "pooled_abs_rmse_mN": (float(np.sqrt(rss_gf / nn) * GF_TO_N * 1e3)
                                   if nn else None),
            "per_hold": per}


def hit_bounds(name: str, popt: Sequence[float]) -> bool:
    if name == "M1":
        tau = popt[2]
        return bool(tau <= 1e-3 + 1e-9 or tau >= 1e4 - 1e-6)
    if name in ("M2", "M3"):
        tau, beta = popt[2], popt[3]
        return bool(tau <= 1e-3 + 1e-9 or tau >= 1e4 - 1e-6
                    or beta <= 0.05 + 1e-9 or beta >= 2.0 - 1e-9)
    return False


def shared_at_bound(name: str, shape: Mapping[str, float | None]) -> bool:
    """True if a shared shape sits on the optimizer/parameter boundary."""
    if name == "M1":
        tau = shape["tau"]
        return bool(tau is None or tau <= 0.1 + 1e-6 or tau >= 1e5 - 1e-3)
    if name in ("M2", "M3"):
        tau, beta = shape["tau"], shape.get("beta")
        return bool(tau is None or beta is None or tau <= 0.1 + 1e-6
                    or tau >= 1e5 - 1e-3 or beta <= 0.05 + 1e-6 or beta >= 2.0 - 1e-6)
    return False


def fit_hold_full(name: str, x: np.ndarray, y: np.ndarray) -> dict:
    """Independent per-hold fit in normalized space (for bound-hit diagnostics)."""
    from scipy.optimize import curve_fit

    if name == "M0":
        return {"params": [float(np.mean(y))], "rmse": float(np.std(y))}
    funcs = {"M1": (m1_exp, 3), "M2": (m2_stretch, 4), "M3": (m3_power, 4)}
    f, _ = funcs[name]
    c0 = float(np.median(y[-max(3, len(y) // 5):]))
    A0 = max(float(y[0] - c0), 1e-3)
    p0 = {"M1": [c0, A0, 20.0], "M2": [c0, A0, 20.0, 0.5],
          "M3": [c0, A0, 20.0, 0.3]}[name]
    bnd = {"M1": ([-np.inf, 0.0, 1e-3], [np.inf, np.inf, 1e4]),
           "M2": ([-np.inf, 0.0, 1e-3, 0.05], [np.inf, np.inf, 1e4, 2.0]),
           "M3": ([-np.inf, 0.0, 1e-3, 0.05], [np.inf, np.inf, 1e4, 2.0])}[name]
    try:
        popt, _ = curve_fit(f, x, y, p0=p0, maxfev=40000, bounds=bnd)
    except Exception as exc:  # noqa: BLE001
        return {"error": str(exc)}
    return {"params": [float(p) for p in popt],
            "rmse": float(np.sqrt(np.mean((y - f(x, *popt)) ** 2)))}


def shape_param_correlation(name: str, popt: Sequence[float], x: np.ndarray) -> float | None:
    """|corr| of the shape params from the numerical Jacobian. DIAGNOSTIC ONLY;
    never used as the sole identifiability gate."""
    if name not in ("M2", "M3"):
        return None
    funcs = {"M2": m2_stretch, "M3": m3_power}
    f = funcs[name]
    popt = np.asarray(popt, float)
    k = popt.size
    J = np.zeros((np.asarray(x).size, k))
    base = f(x, *popt)
    for j in range(k):
        h = 1e-6 * max(abs(popt[j]), 1.0)
        pj = popt.copy()
        pj[j] += h
        J[:, j] = (f(x, *pj) - base) / h
    try:
        cov = np.linalg.inv(J.T @ J)
        d = np.sqrt(np.diag(cov))
        corr = cov / np.outer(d, d)
    except np.linalg.LinAlgError:
        return None
    vals = [abs(float(corr[2, 3]))]
    return max(vals) if vals else None


# --------------------------------------------------------------- hold detector ----
def detect_holds_in_series(t: np.ndarray, force_signed_gf: np.ndarray, quant_gf: float,
                           *, specimen: str | None = None, donor: str | None = None,
                           grp: str | None = None, nmax: int = 300,
                           medfilt_s: float = 1.0, peak_dist_s: float = 30.0,
                           prominence_factor: float = 3.0, min_hold_s: float = 20.0,
                           hold_margin_lo_s: float = 2.0, hold_margin_hi_s: float = 3.0,
                           medfilt_hold_s: float = 0.5) -> list[dict]:
    """Recover step-and-hold segments from a force trace (protocol inference).

    Defaults reproduce the original detector exactly: sign-inverted force, 1 s median
    filter (``medfilt_s``), peak finding with a 3*quantisation prominence
    (``prominence_factor``) and 30 s minimum peak distance (``peak_dist_s``),
    hold window [t_peak+2 s, next_peak-3 s] with >= 50 samples and >= 20 s
    duration (``min_hold_s``), 0.5 s median filter on the hold
    (``medfilt_hold_s``), then uniform decimation to <= ``nmax`` samples.
    The filter/window knobs exist so a sensitivity sweep can perturb the
    detector; changing them is not a claim about the data. Writes nothing.
    """
    from scipy.signal import find_peaks, medfilt

    t = np.asarray(t, float)
    f = np.asarray(force_signed_gf, float)
    m = np.isfinite(t) & np.isfinite(f)
    t, f = t[m], f[m]
    holds: list[dict] = []
    if len(t) < 100:
        return holds
    T = -f
    dt = float(np.median(np.diff(t)))
    w = max(3, int(round(medfilt_s / dt)))
    w += (w % 2 == 0)
    Ts = medfilt(T, w)
    pk, _ = find_peaks(Ts, distance=max(10, int(round(peak_dist_s / dt))),
                       prominence=prominence_factor * quant_gf)
    for i, pi in enumerate(pk):
        tp = float(t[pi])
        sel = (t >= tp - 1.0) & (t <= tp + 5.0)
        if not np.any(sel):
            continue
        tpk = float(t[sel][np.argmax(T[sel])])
        Tpk = float(np.max(T[sel]))
        t_next = float(t[pk[i + 1]]) if i + 1 < len(pk) else float(t[-1])
        lo = tpk + hold_margin_lo_s
        hi = min(t_next - hold_margin_hi_s, float(t[-1]))
        if hi - lo < min_hold_s:
            continue
        selw = (t >= lo) & (t <= hi)
        if np.sum(selw) < 50:
            continue
        x = t[selw] - tpk
        y = T[selw]
        ww = max(3, int(round(medfilt_hold_s / dt)))
        ww += (ww % 2 == 0)
        ys = medfilt(y, ww)
        stride = max(1, int(np.ceil(x.size / max(1, int(nmax)))))
        idx = np.arange(0, x.size, stride)
        if idx.size == 0 or idx[-1] != x.size - 1:
            idx = np.append(idx, x.size - 1)
        holds.append({"grp": grp, "specimen": specimen, "donor": donor,
                      "peak_t_s": tpk, "peak_force_gf": Tpk,
                      "hold_duration_s": float(hi - lo), "n_full": int(np.sum(selw)),
                      "x": x[idx], "y_gf": ys[idx], "y_norm": (ys[idx] / Tpk)})
    return holds


def load_unipd_holds(npz_path: str, meta_path: str, *, nmax: int = 300,
                     **detector_kwargs
                     ) -> tuple[list[dict], dict]:
    """Data adapter over the parsed UniPD artifacts (read-only).

    Consumes ``unipd_series.npz`` + ``unipd_series_meta.json`` (paths supplied by the
    caller). Returns (holds, info). Holds are chronological per specimen,
    with a 1-based ``step`` index matching the original ordering so counts can be
    cross-checked.
    """
    import json

    Z = np.load(npz_path)
    with open(meta_path) as fh:
        meta = json.load(fh)

    diffs = []
    for sheet in ("SR flexor", "SR extensor"):
        for rec in meta["sheets"][sheet]:
            key = sheet.replace(" ", "_") + "/" + rec["id"]
            f = Z[key + "/grp0_Fz_gf"]
            f = f[np.isfinite(f)]
            d = np.abs(np.diff(f))
            d = d[d > 1e-9]
            if len(d):
                diffs.append(d)
    alld = np.concatenate(diffs) if diffs else np.array([np.nan])
    quant = float(np.median(alld[alld > 0]))

    holds: list[dict] = []
    specimens: dict[str, list[dict]] = {}
    for sheet in ("SR flexor", "SR extensor"):
        for rec in meta["sheets"][sheet]:
            sid = rec["id"]
            donor = sid.split("_")[0]
            key = sheet.replace(" ", "_") + "/" + sid
            spec_holds: list[dict] = []
            for g in ("grp0", "grp1"):
                tk, fk = key + "/" + g + "_time_s", key + "/" + g + "_Fz_gf"
                if tk not in Z.files:
                    continue
                spec_holds.extend(detect_holds_in_series(
                    Z[tk], Z[fk], quant, specimen=sid, donor=donor, grp=g,
                    nmax=nmax, **detector_kwargs))
            specimens[sid] = spec_holds
            for step, h in enumerate(spec_holds, start=1):
                h["step"] = step
                h["sheet"] = sheet
                holds.append(h)
    info = {"n_specimens": len(specimens), "n_holds": len(holds),
            "donors": sorted({h["donor"] for h in holds}),
            "n_donors": len({h["donor"] for h in holds}),
            "force_quantization_gf_median_min_nonzero": quant,
            "specimen_hold_counts": {k: len(v) for k, v in specimens.items()}}
    return holds, info


# ------------------------------------------------------- grouped validation ----
def make_folds(groups: Sequence[str], k: int | None = None) -> list[tuple[tuple[str, ...], tuple[str, ...]]]:
    """Deterministic grouped folds. k=None -> leave-one-group-out.

    Returns a list of (train_groups, test_groups) with a strict partition of the
    test groups (each group tested exactly once) and no group in both sides.
    """
    uniq = sorted(set(groups))
    if not uniq:
        raise ValueError("no groups supplied")
    k = len(uniq) if k is None else int(k)
    if k < 2 or k > len(uniq):
        raise ValueError(f"k must be in [2, {len(uniq)}]")
    buckets: list[list[str]] = [[] for _ in range(k)]
    for i, g in enumerate(uniq):
        buckets[i % k].append(g)
    folds = []
    for b in buckets:
        test = tuple(sorted(b))
        train = tuple(sorted(set(uniq) - set(b)))
        assert_no_leakage(train, test)
        folds.append((train, test))
    return folds


def assert_no_leakage(train_groups: Iterable[str], test_groups: Iterable[str]) -> None:
    overlap = set(train_groups) & set(test_groups)
    if overlap:
        raise ValueError(f"group leakage between calibration and validation: "
                         f"{sorted(overlap)}")


def grouped_validation(holds: Sequence[Mapping], models: Sequence[str],
                       group_key: str, k: int | None = None) -> dict:
    """Grouped cross-validation with exact point-weighted pooling.

    Calibration on holds whose ``group_key`` is not in the test groups;
    validation on the test groups only. Per-hold (c,A) are nuisance-solved on
    the held-out hold. Leakage is asserted per fold (calibration and validation
    never share a specimen or donor).
    """
    groups = [str(h[group_key]) for h in holds]
    folds = make_folds(groups, k)
    fold_rows = []
    acc = {name: {"rss_n": 0.0, "rss_gf": 0.0, "nn": 0} for name in models}
    for fi, (train_g, test_g) in enumerate(folds):
        assert_no_leakage(train_g, test_g)
        train = [h for h in holds if str(h[group_key]) in set(train_g)]
        test = [h for h in holds if str(h[group_key]) in set(test_g)]
        if not train or not test:
            raise ValueError(f"fold {fi} empty side: train={len(train)} test={len(test)}")
        row = {"fold": fi, "train_groups": list(train_g), "test_groups": list(test_g),
               "n_train_holds": len(train), "n_test_holds": len(test), "models": {}}
        for name in models:
            shape, rss = fit_shared_shape(name, [(h["x"], h["y_norm"]) for h in train])
            ev = eval_heldout(name, shape, test)
            acc[name]["rss_n"] += (ev["pooled_rel_rmse"] or 0.0) ** 2 * ev["n_points"]
            acc[name]["rss_gf"] += (ev["pooled_abs_rmse_gf"] or 0.0) ** 2 * ev["n_points"]
            acc[name]["nn"] += ev["n_points"]
            row["models"][name] = {
                "shared_shape": shape, "train_rss": rss,
                "pooled_rel_rmse": ev["pooled_rel_rmse"],
                "pooled_abs_rmse_mN": ev["pooled_abs_rmse_mN"]}
        fold_rows.append(row)
    pooled = {}
    for name in models:
        nn = acc[name]["nn"]
        pooled[name] = {
            "pooled_rel_rmse": float(np.sqrt(acc[name]["rss_n"] / nn)) if nn else None,
            "pooled_abs_rmse_mN": (float(np.sqrt(acc[name]["rss_gf"] / nn) * GF_TO_N * 1e3)
                                   if nn else None),
            "n_points": nn,
        }
    return {"group_key": group_key, "n_folds": len(folds), "folds": fold_rows,
            "pooled": pooled}


# --------------------------------------------------- identifiability / gates ----
def shape_param_bootstrap_ci(name: str, curves: Sequence[tuple], groups: Sequence[str],
                             n_boot: int = 200, seed: int = 20260916,
                             alpha: float = 0.05) -> dict:
    """Cluster (group) bootstrap confidence interval for the shared shape params.

    Resamples whole groups with replacement so uncertainty reflects between-unit
    variation, not within-hold pseudo-replication.
    """
    if name not in BASIS_MODELS:
        return {"model": name, "note": "no shape parameter"}
    uniq = sorted(set(groups))
    idx_by_group: dict[str, list[int]] = {g: [] for g in uniq}
    for i, g in enumerate(groups):
        idx_by_group[g].append(i)
    rng = np.random.default_rng(seed)
    taus, betas, failed = [], [], 0
    n_g = len(uniq)
    for _ in range(int(n_boot)):
        pick = rng.integers(0, n_g, size=n_g)
        idx: list[int] = []
        for p in pick:
            idx.extend(idx_by_group[uniq[int(p)]])
        if not idx:
            failed += 1
            continue
        try:
            shape, _ = fit_shared_shape(name, [curves[i] for i in idx])
        except Exception:  # noqa: BLE001
            failed += 1
            continue
        taus.append(shape["tau"])
        if name in ("M2", "M3"):
            betas.append(shape["beta"])

    def ci(vals):
        if not vals:
            return None
        lo, hi = np.percentile(vals, [100 * alpha / 2, 100 * (1 - alpha / 2)])
        return {"n": len(vals), "median": float(np.median(vals)),
                "lo": float(lo), "hi": float(hi),
                "rel_width": float((hi - lo) / np.median(vals)) if np.median(vals) else None}

    return {"model": name, "n_groups": n_g, "n_boot": int(n_boot), "n_failed": failed,
            "tau": ci(taus), "beta": ci(betas)}

def practically_identifiable(*, ci_rel_width: float | None, bound_fraction: float,
                             heldout_margin_frac: float | None,
                             max_abs_shape_corr: float | None = None,
                             max_ci_rel_width: float = 1.0,
                             max_bound_fraction: float = 0.25,
                             min_heldout_margin: float = 0.10) -> dict:
    """Practical identifiability verdict.

    A shape parameter set is reportable only when ALL hold:
      * held-out predictive margin over the simple alternative >= threshold
      * relative CI width <= max_ci_rel_width
      * bound-hit fraction <= max_bound_fraction
    The tau/beta correlation is reported for context but is NOT a gate; a high
    correlation may coexist with a usable predictive form, and a low correlation
    does not by itself make a parameter physically meaningful.
    """
    reasons = []
    if heldout_margin_frac is None or heldout_margin_frac < min_heldout_margin:
        reasons.append("heldout_margin_below_threshold")
    if ci_rel_width is None or ci_rel_width > max_ci_rel_width:
        reasons.append("ci_too_wide_or_unavailable")
    if bound_fraction > max_bound_fraction:
        reasons.append("too_many_bound_hits")
    return {
        "separately_reportable": not reasons,
        "reasons": reasons,
        "criteria": {"min_heldout_margin": min_heldout_margin,
                     "max_ci_rel_width": max_ci_rel_width,
                     "max_bound_fraction": max_bound_fraction},
        "diagnostic_max_abs_shape_corr": max_abs_shape_corr,
        "correlation_is_not_a_gate": True,
    }


PASS, FAIL, ABSTAIN = "PASS", "FAIL", "ABSTAIN"


def assess_geometry_prerequisites(*, area_mm2: Any = None,
                                  reference_length_mm: Any = None,
                                  strain_protocol_verified: bool | None = None,
                                  displacement_channel: Any = None) -> dict:
    """ABSTAIN gate for absolute material parameters.

    Absolute modulus/viscosity/relaxation time requires cross-sectional area,
    reference (gauge) length, a verified strain/displacement protocol and a
    displacement channel. UniPD has NONE of these. Returns an explicit verdict;
    never silently substitutes a default. Uses the installed framework
    ``CouplingEvidence`` ABSTAIN pattern when importable.
    """
    names = ("area_mm2", "reference_length_mm", "strain_protocol_verified",
             "displacement_channel")
    comparisons = {
        "area_mm2": True if area_mm2 is not None else None,
        "reference_length_mm": True if reference_length_mm is not None else None,
        "strain_protocol_verified": True if strain_protocol_verified else None,
        "displacement_channel": True if displacement_channel is not None else None,
    }
    if CouplingEvidence is not None:
        verdict = CouplingEvidence(names).assess(comparisons)
        state = verdict.state if verdict.state != "PASS" else "PASS"
        missing = tuple(verdict.missing)
        backend = "bodytwin.framework.coupling_evidence_v1.CouplingEvidence"
    else:  # standalone fallback with identical semantics
        missing = tuple(k for k in names if comparisons.get(k) is None)
        state = "ABSTAIN" if missing else ("PASS" if all(comparisons.values()) else "FAIL")
        backend = "local_fallback"
    return {"state": state, "missing": list(missing),
            "present": [k for k, v in comparisons.items() if v is True],
            "meaning": ("PASS here means the geometry/protocol PREREQUISITES exist; "
                        "it is not itself a material-parameter identification."),
            "backend": backend}


def _synthetic_check() -> dict:
    """Hermetic check for a standalone run: recover the analytic single-exponential
    reference on noise-free synthetic holds and show the absolute gate ABSTAINs
    without geometry. Synthetic numbers only; no tissue claim."""
    x = np.linspace(0.5, 250.0, 180)
    curves = []
    for s in range(4):
        c, A = 0.2 + 0.05 * s, 0.8 - 0.05 * s
        curves.append((x, c + A * np.exp(-x / TAU_REF_S)))
    shape, rss = fit_shared_shape("M1", curves)
    rel_err = abs(shape["tau"] - TAU_REF_S) / TAU_REF_S
    gate = assess_geometry_prerequisites()
    return {"tau_ref_s": TAU_REF_S, "tau_fit_s": float(shape["tau"]),
            "tau_rel_err": float(rel_err), "tau_within_tol": bool(rel_err < TAU_REL_TOL),
            "absolute_gate_state": gate["state"], "absolute_gate_missing": gate["missing"]}


if __name__ == "__main__":  # pragma: no cover - manual probe
    import json
    import sys as _sys
    _out = _synthetic_check()
    print(json.dumps(_out, indent=2))
    _sys.exit(0 if _out["tau_within_tol"] and _out["absolute_gate_state"] == "ABSTAIN" else 1)
