#!/usr/bin/env python3
"""Tests for the tendon-relaxation module.

Plain python3:  PYTHONPATH=src python tests/test_tendon_relaxation_v1.py
Also pytest-compatible (functions named test_*).

Covers the acceptance items that can be tested in isolation:
  * gf -> N conversion exactness and a declared-frame swap guard (failure test)
  * closed-form profiled LSQ vs numpy least squares (numerical reference)
  * analytic single-exponential recovery at a tolerance fixed in the module
  * constant-data healthy control (constant must win; no richer model invented)
  * grouped folds: disjoint calibration/validation, partition of test groups
  * explicit leakage detector refusal (failure test)
  * grouped validation never shares a specimen between calibration/validation
  * ABSTAIN gate when area/length/strain/displacement are missing
  * normalization is invariant to the force unit/frame (gf vs N)
"""
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
import numpy as np  # noqa: E402

from bodytwin.cells.musculoskeletal import tendon_relaxation_v1 as m  # noqa: E402


def _synthetic_holds(tau=13.7, betas=(0.6, 0.8), n=180, n_spec=4, constant=False,
                     seed=1):
    rng = np.random.default_rng(seed)
    x = np.linspace(0.5, 250.0, n)
    holds = []
    for s in range(n_spec):
        c = 0.2 + 0.05 * s
        A = 0.8 - 0.05 * s
        if constant:
            # a purely elastic (constant-displacement) response has r(x) = 1
            y = np.ones_like(x)
        else:
            beta = betas[s % len(betas)]
            y = c + A * np.exp(-np.power(x / tau, beta))
        holds.append({"specimen": f"sp{s}", "donor": f"d{s % 2}", "step": s + 1,
                      "x": x.copy(), "y_norm": y, "peak_force_gf": 100.0 + s})
    return holds


# --------------------------------------------------------------------- units ----
def test_gf_to_n_exact():
    assert m.GF_TO_N == 9.80665e-3
    got = m.gf_to_newton([1.0, 1000.0])
    assert abs(got[0] - 9.80665e-3) < 1e-15
    assert abs(got[1] - 9.80665) < 1e-12
    rt = m.newton_to_gf(got)
    assert np.allclose(rt, [1.0, 1000.0], rtol=0, atol=1e-12)


def test_unit_swap_rejected():
    # plausible gf data passes
    assert m.assert_declared_force_frame(np.array([500.0, 600.0]), "gf") == "gf"
    # newton-valued array declared as gf -> refused
    try:
        m.assert_declared_force_frame(np.array([5.0, 4.0]), "gf")
        raise AssertionError("expected refusal for newton-valued 'gf'")
    except ValueError:
        pass
    # gf-valued array declared as N -> refused
    try:
        m.assert_declared_force_frame(np.array([5000.0, 3000.0]), "N")
        raise AssertionError("expected refusal for gram-valued 'N'")
    except ValueError:
        pass
    # unknown frame refused
    try:
        m.assert_declared_force_frame(np.array([100.0]), "lbf")
        raise AssertionError("expected refusal for unknown unit")
    except ValueError:
        pass
    # magnitude of the swap for the record
    assert abs((1.0 / m.GF_TO_N) - 101.97162129779) < 1e-6


def test_unit_guard_asymmetry_documented_gap():
    # the magnitude-band guard is asymmetric and has a measured blind spot:
    # physiological newton loads mislabelled as gf are NOT caught. This test
    # locks the limitation so it cannot be silently forgotten.
    try:
        m.assert_declared_force_frame(np.array([5.0, 8.0]), "gf")
        small_newton_not_caught = True
    except ValueError:
        small_newton_not_caught = False
    assert small_newton_not_caught is False  # sub-10 newton values are caught
    try:
        m.assert_declared_force_frame(np.array([5.0, 50.0, 200.0, 1000.0]), "gf")
        phys_newton_not_caught = True
    except ValueError:
        phys_newton_not_caught = False
    assert phys_newton_not_caught is True  # documented blind spot: N -> gf at real loads


def test_short_window_noisy_recovery_degrades():
    # a short observation window with realistic noise leaves the fit residual small
    # but the parameter badly wrong; the recovery error, not the residual, detects it
    tau, beta, quant = m.TAU_REF_S, 0.70, 1.26312
    rng = np.random.default_rng(20260916)
    holds = []
    for s in range(6):
        x = np.linspace(0.5, 3.0, 260)
        P = 2000.0 * (1.0 + 0.15 * s)
        y = 0.36 + 0.64 * np.exp(-np.power(x / tau, beta))
        fg = np.round((P * y + rng.normal(scale=1.0, size=x.size)) / quant) * quant
        holds.append({"specimen": f"s{s}", "donor": f"d{s % 2}", "step": 1,
                      "x": x, "y_norm": fg / P, "peak_force_gf": P})
    curves = [(h["x"], h["y_norm"]) for h in holds]
    shape, _ = m.fit_shared_shape("M2", curves)
    ev = m.eval_heldout("M2", shape, holds)
    tau_err = abs(shape["tau"] - tau) / tau
    assert ev["pooled_rel_rmse"] < 0.01        # residual looks fine
    assert tau_err > 0.15                      # parameter is not


# ------------------------------------------------------- numerical reference ----
def test_lin2_matches_lstsq():
    rng = np.random.default_rng(7)
    for _ in range(20):
        g = rng.normal(size=50)
        y = 3.0 + 2.5 * g + rng.normal(scale=1e-3, size=50)
        (c, A), rss = m.lin2(g, y)
        M = np.column_stack([np.ones_like(g), g])
        coef, *_ = np.linalg.lstsq(M, y, rcond=None)
        ref = y - M @ coef
        assert abs(c - coef[0]) < m.LIN2_REL_TOL * max(1.0, abs(coef[0]))
        assert abs(A - coef[1]) < m.LIN2_REL_TOL * max(1.0, abs(coef[1]))
        assert abs(rss - float(ref @ ref)) < 1e-6 * max(1.0, float(y @ y))


def test_analytic_exp_recovery():
    holds = _synthetic_holds(constant=False, betas=(1.0,))
    curves = [(h["x"], h["y_norm"]) for h in holds]
    shape, _ = m.fit_shared_shape("M1", curves)
    assert abs(shape["tau"] - m.TAU_REF_S) / m.TAU_REF_S < m.TAU_REL_TOL, shape


def test_stretched_analytic_recovery():
    # known stretched form: M2 must recover both shape params tightly on clean data
    tau, beta = 13.7, 0.70
    holds = _synthetic_holds(tau=tau, betas=(beta,), n_spec=6)
    curves = [(h["x"], h["y_norm"]) for h in holds]
    shape, _ = m.fit_shared_shape("M2", curves)
    assert abs(shape["tau"] - tau) / tau < 0.05, shape
    assert abs(shape["beta"] - beta) / beta < 0.05, shape


def test_shared_at_bound_flags():
    assert m.shared_at_bound("M2", {"tau": 0.05, "beta": 0.7}) is True
    assert m.shared_at_bound("M2", {"tau": 5.0, "beta": 2.0}) is True
    assert m.shared_at_bound("M2", {"tau": 5.0, "beta": 0.7}) is False
    assert m.shared_at_bound("M1", {"tau": 5.0, "beta": None}) is False


def _synthetic_raw_trace(tau=13.7, beta=0.7, c_frac=0.35, dt=0.03, duration=410.0,
                         peaks=(10.0, 210.0), amps=(1000.0, 1500.0), seed=3):
    """Two step-and-hold peaks; stored with G's sign convention (tension = -Fz)."""
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, duration, dt)
    tension = np.full_like(t, c_frac * amps[0])
    for i, (tp, P) in enumerate(zip(peaks, amps)):
        t_end = peaks[i + 1] if i + 1 < len(peaks) else duration
        sel = (t >= tp) & (t < t_end)
        x = t[sel] - tp
        tension[sel] = c_frac * P + (1 - c_frac) * P * np.exp(-np.power(x / tau, beta))
    fz = -(tension + rng.normal(scale=0.5, size=t.size))  # signed so T = -fz
    return t, fz


def test_detector_finds_synthetic_holds():
    t, fz = _synthetic_raw_trace()
    holds = m.detect_holds_in_series(t, fz, 1.26312, specimen="syn", donor="d",
                                     grp="grp0")
    assert len(holds) == 2, [h["peak_force_gf"] for h in holds]
    # decimation bound respected
    assert all(h["x"].size <= 300 for h in holds)
    holds_lo = m.detect_holds_in_series(t, fz, 1.26312, specimen="syn", donor="d",
                                        grp="grp0", nmax=150)
    assert all(h["x"].size <= 150 for h in holds_lo)
    # normalized shape starts near 1 (median-filtered) and decays
    for h in holds:
        assert h["y_norm"][0] <= 1.0 + 1e-9
        assert h["y_norm"][0] > 0.8
        assert h["y_norm"][-1] < h["y_norm"][0]


def test_constant_data_control():
    holds = _synthetic_holds(constant=True)
    curves = [(h["x"], h["y_norm"]) for h in holds]
    _, rss0 = m.fit_shared_shape("M0", curves)
    ev0 = m.eval_heldout("M0", {}, holds)
    ev1 = m.eval_heldout("M1", m.fit_shared_shape("M1", curves)[0], holds)
    assert ev0["pooled_rel_rmse"] < 1e-9
    assert ev0["pooled_rel_rmse"] <= ev1["pooled_rel_rmse"] + 1e-9
    # healthy control: a richer model must not be declared justified here
    margin = ((ev1["pooled_rel_rmse"] - ev0["pooled_rel_rmse"])
              / ev1["pooled_rel_rmse"]) if ev1["pooled_rel_rmse"] else 0.0
    assert margin <= 0.0


# ---------------------------------------------------------- grouped validation ----
def test_grouped_folds_disjoint_and_partition():
    groups = ["a", "b", "c", "d", "e"]
    folds = m.make_folds(groups)
    assert len(folds) == 5
    tested = []
    for train, test in folds:
        assert not (set(train) & set(test))
        tested.extend(test)
    assert sorted(tested) == sorted(groups)


def test_leakage_detector_raises():
    try:
        m.assert_no_leakage(["a", "b"], ["b", "c"])
        raise AssertionError("expected leakage refusal")
    except ValueError:
        pass
    m.assert_no_leakage(["a"], ["b"])  # must not raise


def test_grouped_validation_never_shares_specimen():
    holds = _synthetic_holds(n_spec=4)
    out = m.grouped_validation(holds, ("M0", "M1"), group_key="specimen")
    assert out["n_folds"] == 4
    for fold in out["folds"]:
        assert not (set(fold["train_groups"]) & set(fold["test_groups"]))
        assert set(fold["test_groups"])


def test_grouped_validation_synthetic_recovers_richer():
    # curves generated with a stretched form: grouped held-out should prefer M2
    holds = _synthetic_holds(tau=13.7, betas=(0.55, 0.75), n_spec=6)
    out = m.grouped_validation(holds, ("M0", "M1", "M2"), group_key="specimen")
    rel = {k: out["pooled"][k]["pooled_rel_rmse"] for k in out["pooled"]}
    assert rel["M2"] < rel["M1"] < rel["M0"], rel


# -------------------------------------------------------------- abstain gate ----
def test_absolute_gate_abstains_without_geometry():
    v = m.assess_geometry_prerequisites()
    assert v["state"] == "ABSTAIN"
    assert set(v["missing"]) == {"area_mm2", "reference_length_mm",
                                 "strain_protocol_verified", "displacement_channel"}
    # full prerequisites -> PASS means only that the prerequisites exist
    v2 = m.assess_geometry_prerequisites(area_mm2=57.5, reference_length_mm=200.0,
                                         strain_protocol_verified=True,
                                         displacement_channel=100.0)
    assert v2["state"] == "PASS"
    assert "not itself" in v2["meaning"]
    # partial -> ABSTAIN, never PASS
    v3 = m.assess_geometry_prerequisites(area_mm2=57.5)
    assert v3["state"] == "ABSTAIN"


# ---------------------------------------------------- normalization invariance ----
def test_normalization_unit_invariance():
    holds = _synthetic_holds(n_spec=2)
    scaled = []
    for h in holds:
        h2 = dict(h)
        h2["y_norm"] = h["y_norm"]  # already normalized; rebuild from a unit change
        scaled.append(h2)
    # build raw gf traces, then show the normalized shape is unchanged when the
    # same trace is expressed in N
    raw_gf = np.array([500.0, 480.0, 470.0, 465.0])
    peak = raw_gf[0]
    r_gf = (-raw_gf) / (-peak)
    r_n = m.gf_to_newton(-raw_gf) / m.gf_to_newton(-peak)
    assert np.allclose(r_gf, r_n, rtol=0, atol=1e-15)


def test_reference_tolerances_declared():
    assert m.TAU_REL_TOL == 1e-3
    assert m.LIN2_REL_TOL == 1e-9
    assert m.TAU_REF_S > 0


def _run_all():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    results, failed = [], []
    for t in tests:
        try:
            t()
            results.append((t.__name__, "PASS"))
        except Exception as exc:  # noqa: BLE001
            results.append((t.__name__, f"FAIL: {exc}"))
            failed.append(t.__name__)
    for name, res in results:
        if res == "PASS":
            print(f"{res:60s} {name}")
        else:
            print(f"{name}: {res}")
    print(f"\n{len(tests) - len(failed)}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(_run_all())
