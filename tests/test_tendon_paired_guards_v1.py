"""Tests for the step-validity guards in ``tendon_paired_analysis_v1``.

Each guard gets a planted bad step that must return REJECT (pre-fit) or ABSTAIN
(post-fit) with the named reason, and no ``tau_s``. A clean synthetic
ramp-hold step must still PASS. Synthetic records use the module's own exact
single-tau SLS recurrence (known linear viscoelasticity, no new model).
"""
from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
from bodytwin.cells.musculoskeletal import tendon_paired_analysis_v1 as m  # noqa: E402

AREA_MM2 = 57.5
REF_MM = 200.0
E_INF = 0.8e9
E_BRANCH = 0.4e9


def _strain_ramp_hold(t, ramp_s, eps0):
    return np.where(t <= ramp_s, eps0 * t / ramp_s, eps0)


def _inputs(t, eps, stress):
    return m.PairedInputs(
        time_s=t, force_gf=stress * AREA_MM2 * 1e-6 / m.m.GF_TO_N,
        displacement_mm=eps * REF_MM, reference_length_mm=REF_MM,
        area_mm2=AREA_MM2, specimen="syn", donor="syn")


def _record(tau=30.0, end_s=300.0, n=3001, ramp_s=1.0, eps0=0.02):
    t = np.linspace(0.0, end_s, n)
    eps = _strain_ramp_hold(t, ramp_s, eps0)
    stress = m._sls_stress_from_strain(t, eps, E_INF, E_BRANCH, tau)
    return t, eps, stress


def _assert_rejected(fit, reason, state):
    assert fit["state"] == state, fit
    assert reason in fit["reasons"], fit["reasons"]
    assert "tau_s" not in fit, fit
    assert fit["guard_thresholds"] == m.GUARD_THRESHOLDS


# ---------------------------------------------------------------- clean PASS ----
def test_clean_step_passes_with_guard_metrics():
    t, eps, stress = _record()
    fit = m.fit_tau_from_pairs(_inputs(t, eps, stress))
    assert fit["state"] == "PASS", fit
    assert abs(fit["tau_s"] - 30.0) / 30.0 < 1e-6
    gm = fit["guard_metrics"]
    assert gm["retraction_frac"] == 0.0
    assert gm["hold_over_tau"] > m.GUARD_MIN_HOLD_OVER_TAU
    assert gm["ci_rel_width"] < m.GUARD_MAX_CI_REL_WIDTH
    assert gm["rmse_rel"] < m.GUARD_MAX_RMSE_REL


def test_drifting_step_still_passes():
    # +/-2 % hold drift is a protocol flag, not a reason to reject the step
    for frac in (-0.02, 0.02):
        t = np.linspace(0.0, 300.0, 3001)
        eps = _strain_ramp_hold(t, 1.0, 0.02)
        hold = t > 1.0
        eps[hold] = 0.02 * (1.0 + frac * (t[hold] - 1.0) / 299.0)
        stress = m._sls_stress_from_strain(t, eps, E_INF, E_BRANCH, 30.0)
        fit = m.fit_tau_from_pairs(_inputs(t, eps, stress))
        assert fit["state"] == "PASS", (frac, fit)
        assert abs(fit["tau_s"] - 30.0) / 30.0 < 1e-6


def test_guards_false_restores_unguarded_fit():
    t, eps, stress = _record(end_s=60.0, n=601)
    eps = np.where(t > 40.0, 0.02 * np.clip((60.0 - t) / 20.0, 0, 1), eps)
    stress = m._sls_stress_from_strain(t, eps, E_INF, E_BRANCH, 30.0)
    inp = _inputs(t, eps, stress)
    assert m.fit_tau_from_pairs(inp)["state"] == "REJECT"
    assert m.fit_tau_from_pairs(inp, guards=False)["state"] == "PASS"


# ------------------------------------------------------------- pre-fit REJECT ----
def test_reject_too_few_samples():
    t, eps, stress = _record(end_s=300.0, n=15)
    _assert_rejected(m.fit_tau_from_pairs(_inputs(t, eps, stress)),
                     "too_few_samples", "REJECT")


def test_reject_time_not_strictly_increasing():
    t, eps, stress = _record()
    t = t.copy()
    t[1500] = t[1499]  # duplicated time stamp
    _assert_rejected(m.fit_tau_from_pairs(_inputs(t, eps, stress)),
                     "time_not_strictly_increasing", "REJECT")


def test_reject_no_tensile_loading():
    # compression-only step: strain never rises above its first sample
    t, eps, stress = _record()
    _assert_rejected(m.fit_tau_from_pairs(_inputs(t, -eps, -stress)),
                     "no_tensile_loading", "REJECT")


def test_reject_unloading_step():
    # ramp, hold, then full unload (the catch-matrix fixture shape)
    t = np.linspace(0.0, 60.0, 6001)
    eps = _strain_ramp_hold(t, 1.0, 0.02)
    eps = np.where(t > 40.0, 0.02 * np.clip((60.0 - t) / 20.0, 0, 1), eps)
    stress = m._sls_stress_from_strain(t, eps, E_INF, E_BRANCH, 30.0)
    fit = m.fit_tau_from_pairs(_inputs(t, eps, stress))
    _assert_rejected(fit, "unloading_after_peak", "REJECT")
    assert fit["guard_metrics"]["retraction_frac"] > m.GUARD_MAX_RETRACTION_FRAC


def test_unloading_threshold_boundary():
    # 5 % retraction (above every real accepted step, 3.25 %) still passes;
    # 15 % is rejected
    for frac, state in ((0.05, "PASS"), (0.15, "REJECT")):
        t = np.linspace(0.0, 300.0, 3001)
        eps = _strain_ramp_hold(t, 1.0, 0.02)
        eps = np.where(t > 250.0, 0.02 * (1.0 - frac * np.clip((t - 250.0) / 50.0, 0, 1)), eps)
        stress = m._sls_stress_from_strain(t, eps, E_INF, E_BRANCH, 30.0)
        fit = m.fit_tau_from_pairs(_inputs(t, eps, stress))
        assert fit["state"] == state, (frac, fit.get("reasons"))


def test_reject_no_hold_pure_ramp():
    # loading ramp to the last sample, no hold at all
    t = np.linspace(0.0, 25.0, 2501)
    eps = 0.04 * t / 25.0
    stress = m._sls_stress_from_strain(t, eps, E_INF, E_BRANCH, 30.0)
    _assert_rejected(m.fit_tau_from_pairs(_inputs(t, eps, stress)),
                     "no_hold", "REJECT")


def test_reject_force_collapse_below_step_start():
    # hold whose force falls far below the step-start force (rupture-like)
    t, eps, stress = _record()
    stress = stress.copy()
    peak = float(np.max(stress))
    stress[t > 200.0] = -0.5 * peak
    fit = m.fit_tau_from_pairs(_inputs(t, eps, stress))
    _assert_rejected(fit, "force_collapse_below_step_start", "REJECT")
    assert fit["guard_metrics"]["end_force_frac"] < m.GUARD_MIN_END_FORCE_FRAC


# ------------------------------------------------------------ post-fit ABSTAIN ----
def test_abstain_hold_shorter_than_two_tau():
    # tau = 100 s but only a 120 s hold
    t, eps, stress = _record(tau=100.0, end_s=121.0, n=1211)
    fit = m.fit_tau_from_pairs(_inputs(t, eps, stress))
    _assert_rejected(fit, "hold_shorter_than_two_tau", "ABSTAIN")
    assert abs(fit["unaccepted_fit"]["tau_s"] - 100.0) / 100.0 < 1e-3


def test_abstain_tau_below_sampling_resolution():
    # tau = 1 s sampled every 0.5 s: less than 5 samples per tau
    t, eps, stress = _record(tau=1.0, end_s=300.0, n=601, ramp_s=0.5)
    fit = m.fit_tau_from_pairs(_inputs(t, eps, stress))
    _assert_rejected(fit, "tau_below_sampling_resolution", "ABSTAIN")


def _noisy(t, stress, frac, seed):
    inc = float(np.max(stress) - stress[0])
    out = stress + np.random.default_rng(seed).normal(0.0, frac * inc, size=stress.size)
    out[0] = stress[0]
    return out


def test_abstain_single_tau_misfit_end_to_end():
    # 15 % stress noise on a long, well-sampled hold: tau CI stays narrow but the
    # residual exceeds 10 % of the stress increment
    t, eps, stress = _record(tau=30.0)
    fit = m.fit_tau_from_pairs(_inputs(t, eps, _noisy(t, stress, 0.15, 3)))
    _assert_rejected(fit, "single_tau_sls_misfit", "ABSTAIN")
    assert fit["guard_metrics"]["rmse_rel"] > m.GUARD_MAX_RMSE_REL


def test_abstain_tau_ci_too_wide_end_to_end():
    # 71 samples, 8 % noise: residual is acceptable but tau is not identified
    t, eps, stress = _record(tau=30.0, end_s=70.0, n=71)
    fit = m.fit_tau_from_pairs(_inputs(t, eps, _noisy(t, stress, 0.08, 1)))
    _assert_rejected(fit, "tau_ci_too_wide", "ABSTAIN")
    assert fit["guard_metrics"]["ci_rel_width"] > m.GUARD_MAX_CI_REL_WIDTH


def test_tau_identifiability_each_reason_direct():
    ok = dict(tau=30.0, ci=[29.0, 31.0], rmse=1.0, stress_increment=100.0,
              median_dt_s=0.1, hold_duration_s=300.0)
    assert m.tau_identifiability(**ok) == []
    cases = {
        "tau_at_fit_bound": dict(tau=1e-9, hold_duration_s=300.0),
        "tau_ci_not_finite": dict(ci=None),
        "tau_ci_too_wide": dict(ci=[1.0, 45.0]),
        "single_tau_sls_misfit": dict(rmse=11.0),
        "tau_below_sampling_resolution": dict(median_dt_s=7.0),
        "hold_shorter_than_two_tau": dict(hold_duration_s=59.0),
    }
    for reason, patch in cases.items():
        kw = dict(ok)
        kw.update(patch)
        assert reason in m.tau_identifiability(**kw), (reason, kw)
    nan_ci = dict(ok, ci=[float("nan"), 31.0])
    assert "tau_ci_not_finite" in m.tau_identifiability(**nan_ci)


def test_real_data_worst_values_inside_thresholds():
    w = m.GUARD_REAL_DATA_WORST
    assert w["retraction_frac_max"] < m.GUARD_MAX_RETRACTION_FRAC
    assert w["end_force_frac_min"] > m.GUARD_MIN_END_FORCE_FRAC
    assert w["ci_rel_width_max"] < m.GUARD_MAX_CI_REL_WIDTH
    assert w["rmse_rel_max"] < m.GUARD_MAX_RMSE_REL
    assert w["tau_over_dt_min"] > m.GUARD_MIN_TAU_OVER_DT
    assert w["hold_over_tau_min"] > m.GUARD_MIN_HOLD_OVER_TAU
    assert w["hold_over_loading_min"] > m.GUARD_MIN_HOLD_OVER_LOADING


def test_hold_start_frac_robust_to_positive_drift():
    t = np.linspace(0.0, 300.0, 3001)
    eps = _strain_ramp_hold(t, 1.0, 0.02)
    eps[t > 1.0] = 0.02 * (1.0 + 0.02 * (t[t > 1.0] - 1.0) / 299.0)
    h = m.hold_start_index_frac(eps)
    assert t[h] <= 1.0, t[h]


def test_guards_scale_invariant_to_units():
    # x1000 displacement / force (unit-frame error) must not change the verdict
    t, eps, stress = _record()
    for dscale, fscale in ((1e3, 1.0), (1.0, 1e3)):
        inp = m.PairedInputs(
            time_s=t, force_gf=stress * fscale * AREA_MM2 * 1e-6 / m.m.GF_TO_N,
            displacement_mm=eps * REF_MM * dscale, reference_length_mm=REF_MM,
            area_mm2=AREA_MM2)
        assert m.fit_tau_from_pairs(inp)["state"] == "PASS"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    bad = 0
    for fn in fns:
        try:
            fn()
            print("PASS", fn.__name__)
        except Exception as exc:  # noqa: BLE001
            bad += 1
            print("FAIL", fn.__name__, exc)
    print(f"{len(fns) - bad}/{len(fns)} passed")
    sys.exit(1 if bad else 0)
