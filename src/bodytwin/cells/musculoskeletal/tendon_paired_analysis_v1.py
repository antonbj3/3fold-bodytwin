"""Paired force + displacement tendon analysis (data-anchored late mechanics).

Why this module exists
----------------------
The normalized relaxation analysis in ``tendon_relaxation_v1`` (force + time only) cannot separate
material time-dependence from an unverified hold-displacement protocol: any
nonzero hold-strain drift makes a single-tau standard-linear-solid (SLS) look as
if it needed a two-parameter stretched form (a synthetic protocol-bias check during
development; that evidence file is not part of this repository). The decisive missing anchor is a PAIRED
force + displacement channel (plus specimen geometry: cross-sectional area and a
reference length).

This module is the analysis path that becomes usable the moment such paired data
exists. It is deliberately honest about what is identifiable:

* No paired displacement / no geometry  ->  absolute parameters ABSTAIN.
* Force + displacement but no geometry   ->  normalized shape + protocol
  diagnostics only; no absolute stress.
* Force + displacement + geometry        ->  absolute stress and the SLS
  relaxation time ``tau`` become recoverable, because the ACTUAL observed strain
  history (not an assumed constant-hold) is propagated through the exact SLS.

The point of paired data is that ``tau`` is recoverable even while the hold
strain drifts: the drift is observed, hence it is a known input, not a confound.

Design and honesty boundary
---------------------------
* Pure Python / NumPy / SciPy. Performs no I/O and no writes; the caller supplies
  arrays explicitly.
* Units are explicit and pinned to the sibling analyzer's exact factor:
  ``GF_TO_N`` is reused from the sibling module ``tendon_relaxation_v1``
  (package import; a standalone script run adds the src root to ``sys.path``).
* ``PASS`` at a prerequisite gate means the prerequisites EXIST; it is never, by
  itself, a material-parameter identification.
* No extrapolation beyond the observed strain range.

Attribution: the hold detector and the geometry-prerequisite gate are reused
read-only from the sibling module ``tendon_relaxation_v1.py``. The single-tau
SLS propagation is standard linear viscoelasticity (the same exact per-step
recurrence as an SLS reference script used during development, which is not part of
this repository); it is known mathematics, not a new tissue model. The paired input
layout follows a synthetic paired-dataset description from the same development work.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from dataclasses import dataclass
from typing import Any

import numpy as np

# ------------------------------------------------------------- sibling module ----
# Reuse the relaxation analyzer (read-only). As a package module the sibling is
# imported normally; when this file is run as a standalone script (the repo's
# cell runner does that) the src root is put on sys.path and the package import
# is used as well, so both paths share one module object.
_HERE = os.path.dirname(os.path.abspath(__file__))
if __package__:
    from . import tendon_relaxation_v1 as m
else:  # standalone script run
    _SRC_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
    if _SRC_ROOT not in sys.path:
        sys.path.insert(0, _SRC_ROOT)
    from bodytwin.cells.musculoskeletal import tendon_relaxation_v1 as m

__all__ = [
    "MM_PER_M", "MM2_PER_M2",
    "gf_to_n", "mm_to_m", "mm2_to_m2", "strain_from_displacement",
    "PairedInputs", "verify_constant_displacement", "protocol_diagnostic",
    "normalized_relaxation_holds", "absolute_parameter_gate", "apparent_modulus",
    "fit_tau_from_pairs", "step_validity", "tau_identifiability", "self_test",
    "REJECT", "GUARD_THRESHOLDS",
]

# ------------------------------------------------------------------- units ----
# 1 gram-force is, by definition, the weight of 1 gram under standard gravity;
# the exact factor lives in the sibling analyzer and is reused, never re-derived.
MM_PER_M = 1e-3
MM2_PER_M2 = 1e-6

PASS, FAIL, ABSTAIN = "PASS", "FAIL", "ABSTAIN"
REJECT = "REJECT"

# ------------------------------------------------------ step-validity guards ----
# ``fit_tau_from_pairs`` only returns PASS for a loading-then-hold relaxation step
# whose single tau is identified. Two layers:
#   * pre-fit  -> REJECT  : the record is not a relaxation step at all;
#   * post-fit -> ABSTAIN : a relaxation step, but tau is not identified.
# Thresholds are modelling choices (assumptions), not validated physiological
# limits. They were set against the real animal anchors that currently pass
# (rat Zenodo 10.5281/zenodo.7371107, 19 usable steps; porcine Leeds
# 10.5518/1802, 18 specimens), so a threshold can only be as general as those
# two datasets. The observed worst real value for each metric is recorded in
# ``GUARD_REAL_DATA_WORST``. Those margins were measured on 2026-09-17 during development;
# the measurement is not reproducible from this repository, so treat them as recorded
# assumptions, not as verified values.
GUARD_MIN_SAMPLES = 20            # 3 free parameters; >= ~6 samples per parameter
GUARD_HOLD_START_FRAC = 0.95      # hold starts where strain first reaches 95 % of its peak increment
GUARD_MIN_HOLD_SAMPLES = 10       # a hold needs samples to show relaxation at all
GUARD_MIN_HOLD_OVER_LOADING = 0.5 # hold must last >= half the time spent reaching it
GUARD_MAX_RETRACTION_FRAC = 0.10  # strain may fall at most 10 % of its increment after the peak
GUARD_MIN_END_FORCE_FRAC = -0.20  # end force may undershoot step-start force by <= 20 % of the increment
GUARD_MAX_CI_REL_WIDTH = 1.0      # 95 % CI width / tau; > 1 means the CI is wider than tau itself
GUARD_MAX_RMSE_REL = 0.10         # rmse / incremental stress range (same limit as the import wrappers)
GUARD_MIN_TAU_OVER_DT = 5.0       # tau must span >= 5 median sample intervals
GUARD_MIN_HOLD_OVER_TAU = 2.0     # hold >= 2 tau: 86.5 % of the branch stress has relaxed
GUARD_TAU_BOUNDS = (1e-9, 1e7)    # curve_fit bounds on tau; a fit at a bound is not identified

GUARD_THRESHOLDS = {
    "min_samples": GUARD_MIN_SAMPLES,
    "hold_start_strain_frac": GUARD_HOLD_START_FRAC,
    "min_hold_samples": GUARD_MIN_HOLD_SAMPLES,
    "min_hold_over_loading": GUARD_MIN_HOLD_OVER_LOADING,
    "max_retraction_frac": GUARD_MAX_RETRACTION_FRAC,
    "min_end_force_frac": GUARD_MIN_END_FORCE_FRAC,
    "max_ci_rel_width": GUARD_MAX_CI_REL_WIDTH,
    "max_rmse_rel": GUARD_MAX_RMSE_REL,
    "min_tau_over_dt": GUARD_MIN_TAU_OVER_DT,
    "min_hold_over_tau": GUARD_MIN_HOLD_OVER_TAU,
    "tau_bounds_s": list(GUARD_TAU_BOUNDS),
}

# Worst value among the 37 real steps that pass (19 rat + 18 porcine), as recorded during
# development before the guards were added (not reproducible from this repository).
# Used only for documentation and tests.
GUARD_REAL_DATA_WORST = {
    "retraction_frac_max": 0.0325,      # rat anchor (limit 0.10)
    "end_force_frac_min": -0.0866,      # rat anchor (limit -0.20)
    "ci_rel_width_max": 0.152,          # rat anchor (limit 1.0)
    "rmse_rel_max": 0.0662,             # rat anchor (limit 0.10)
    "tau_over_dt_min": 178.0,           # rat anchor, 17.8 s / 0.1 s (limit 5)
    "hold_over_tau_min": 2.37,          # porcine anchor, 324.8 s / 137.2 s (limit 2)
    "hold_over_loading_min": 2.28,      # porcine anchor, 300.6 s / 131.8 s (limit 0.5)
}


def gf_to_n(force_gf: Any) -> np.ndarray:
    """Convert gram-force to newton via the sibling analyzer's exact factor."""
    return np.asarray(force_gf, dtype=float) * m.GF_TO_N


def mm_to_m(length_mm: Any) -> np.ndarray:
    """Convert millimetre to metre (exact factor 1e-3 m/mm)."""
    return np.asarray(length_mm, dtype=float) * MM_PER_M


def mm2_to_m2(area_mm2: Any) -> np.ndarray:
    """Convert square millimetre to square metre (exact factor 1e-6 m^2/mm^2)."""
    return np.asarray(area_mm2, dtype=float) * MM2_PER_M2


def strain_from_displacement(displacement_mm: Any,
                             reference_length_mm: float | None) -> np.ndarray:
    """Engineering strain = displacement / reference (gauge) length.

    Both arguments are in millimetre; the ratio is dimensionless. Refuses a
    nonpositive or non-finite reference length instead of substituting a default.
    """
    disp = np.asarray(displacement_mm, dtype=float)
    if not np.all(np.isfinite(disp)):
        raise ValueError("displacement_mm contains non-finite values")
    if reference_length_mm is None:
        raise ValueError("reference_length_mm is required to form strain")
    length = float(reference_length_mm)
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError(
            "reference_length_mm must be finite and > 0, got "
            f"{reference_length_mm!r}")
    return disp / length


# ------------------------------------------------------------------- inputs ----
@dataclass(frozen=True)
class PairedInputs:
    """Immutable paired force + displacement record (one specimen/step).

    ``time_s`` and ``force_gf`` are the channels the old UniPD analysis had.
    ``displacement_mm`` is the decisive added anchor. ``reference_length_mm`` and
    ``area_mm2`` are the specimen geometry; ``None`` means "not measured" and is
    never replaced by a default. ``provenance`` (optional) may carry an explicit
    ``strain_protocol_verified`` flag; see :func:`absolute_parameter_gate`.

    The constructor validates shape/finiteness and copies the arrays, so the
    caller's inputs are never mutated.
    """

    time_s: np.ndarray
    force_gf: np.ndarray
    displacement_mm: np.ndarray
    reference_length_mm: float | None = None
    area_mm2: float | None = None
    specimen: str = ""
    donor: str = ""
    provenance: dict | None = None

    def __post_init__(self) -> None:
        t = np.asarray(self.time_s, dtype=float)
        f = np.asarray(self.force_gf, dtype=float)
        d = np.asarray(self.displacement_mm, dtype=float)
        if t.ndim != 1 or f.ndim != 1 or d.ndim != 1:
            raise ValueError("time_s, force_gf and displacement_mm must be 1-D")
        if not (t.size == f.size == d.size):
            raise ValueError(
                "time_s, force_gf and displacement_mm must have equal length, got "
                f"{t.size}, {f.size}, {d.size}")
        if t.size == 0:
            raise ValueError("paired inputs must be non-empty")
        for name, arr in (("time_s", t), ("force_gf", f), ("displacement_mm", d)):
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} contains non-finite values")
        if self.reference_length_mm is not None:
            length = float(self.reference_length_mm)
            if not np.isfinite(length) or length <= 0.0:
                raise ValueError("reference_length_mm must be finite and > 0")
        if self.area_mm2 is not None:
            area = float(self.area_mm2)
            if not np.isfinite(area) or area <= 0.0:
                raise ValueError("area_mm2 must be finite and > 0")
        prov = dict(self.provenance) if isinstance(self.provenance, dict) else self.provenance
        # Store copies so we neither alias nor mutate the caller's arrays/objects.
        object.__setattr__(self, "time_s", t.copy())
        object.__setattr__(self, "force_gf", f.copy())
        object.__setattr__(self, "displacement_mm", d.copy())
        object.__setattr__(self, "reference_length_mm",
                           None if self.reference_length_mm is None
                           else float(self.reference_length_mm))
        object.__setattr__(self, "area_mm2",
                           None if self.area_mm2 is None else float(self.area_mm2))
        object.__setattr__(self, "provenance", prov)


# ------------------------------------------------- protocol verification ----
def verify_constant_displacement(displacement_mm: Any, time_s: Any,
                                 hold_start_s: float, hold_end_s: float,
                                 rel_tol: float = 0.01) -> dict:
    """Decide whether the hold really is (near) constant displacement.

    This is the diagnostic that was impossible with force-only data. It works on
    a displacement **or** a strain series (both are dimensionless-ratio friendly).

    Classification (priority order):
      * ``nonmonotonic``  the hold is not monotone within the tolerance;
      * ``creep``         monotone and the relative drift is >= +``rel_tol``;
      * ``slip``          monotone and the relative drift is <= -``rel_tol``;
      * ``constant``      otherwise.

    ``drift_abs_mm`` is end-minus-start over the window and ``drift_rel`` is that
    difference relative to the first hold sample; ``slope_mm_per_s`` is the
    ordinary least-squares slope. When applied to strain the unit labels remain
    "mm" by name only; they are the same dimensionless ratio either way.
    """
    d = np.asarray(displacement_mm, dtype=float)
    t = np.asarray(time_s, dtype=float)
    if d.ndim != 1 or t.ndim != 1 or d.size != t.size:
        raise ValueError("displacement_mm and time_s must be equal-length 1-D arrays")
    if not np.all(np.isfinite(d)) or not np.all(np.isfinite(t)):
        raise ValueError("displacement_mm/time_s contain non-finite values")
    lo, hi = float(hold_start_s), float(hold_end_s)
    if not hi > lo:
        raise ValueError("hold_end_s must be greater than hold_start_s")
    tol = float(rel_tol)
    if not np.isfinite(tol) or tol < 0.0:
        raise ValueError("rel_tol must be finite and >= 0")
    sel = (t >= lo) & (t <= hi)
    if np.count_nonzero(sel) < 2:
        raise ValueError("hold window contains fewer than two samples")
    th, dh = t[sel], d[sel]
    slope = float(np.polyfit(th, dh, 1)[0])
    drift_abs = float(dh[-1] - dh[0])
    base = float(dh[0])
    drift_rel = float(drift_abs / base) if base != 0.0 else 0.0
    diffs = np.diff(dh)
    # a monotone check tolerating only numerical noise: a tenth of the absolute
    # drift that would change the classification at all (rel_tol * |base|).
    mono_tol = 0.1 * tol * max(abs(base), 1e-15)
    mono_up = bool(np.all(diffs >= -mono_tol))
    mono_dn = bool(np.all(diffs <= mono_tol))
    # the small epsilon keeps a protocol defined at exactly rel_tol (e.g. +1%)
    # from being classified "constant" by floating-point round-off
    eff_tol = max(tol - 1e-12, 0.0)
    if not (mono_up or mono_dn):
        classification = "nonmonotonic"
    elif drift_rel > 0.0 and drift_rel >= eff_tol:
        classification = "creep"
    elif drift_rel < 0.0 and drift_rel <= -eff_tol:
        classification = "slip"
    else:
        classification = "constant"
    return {
        "constant": classification == "constant",
        "drift_abs_mm": drift_abs,
        "drift_rel": drift_rel,
        "slope_mm_per_s": slope,
        "classification": classification,
        "hold_start_s": lo,
        "hold_end_s": hi,
        "n_hold_samples": int(np.count_nonzero(sel)),
        "rel_tol": tol,
    }


def protocol_diagnostic(inputs: PairedInputs, ramp_end_s: float,
                        hold_end_s: float, rel_tol: float = 0.01) -> dict:
    """Apply the constant-displacement check to a paired record's hold.

    When a reference length is present the check runs on STRAIN; otherwise it
    runs on the raw displacement in millimetre and says so explicitly (the
    normalized shape can still be analyzed, but no strain claim is made).
    """
    if inputs.reference_length_mm is not None:
        series = strain_from_displacement(inputs.displacement_mm,
                                          inputs.reference_length_mm)
        channel = "strain"
        used_strain = True
    else:
        series = inputs.displacement_mm
        channel = "displacement_mm"
        used_strain = False
    verdict = verify_constant_displacement(series, inputs.time_s, ramp_end_s,
                                           hold_end_s, rel_tol=rel_tol)
    out = dict(verdict)
    out.update({
        "channel": channel,
        "used_strain": used_strain,
        "reference_length_mm": inputs.reference_length_mm,
        "specimen": inputs.specimen,
        "donor": inputs.donor,
        "note": ("constant-displacement check on "
                 + ("strain (reference length present)" if used_strain
                    else "raw displacement in mm (no reference length present)")),
    })
    return out


# --------------------------------------------------- normalized holds + pairs ----
def _force_quantization_gf(force_gf: np.ndarray) -> float:
    """Median absolute non-zero first difference; a robust detector prominence."""
    diffs = np.abs(np.diff(np.asarray(force_gf, dtype=float)))
    diffs = diffs[diffs > 1e-12]
    return float(np.median(diffs)) if diffs.size else 1e-9


def _hold_window_indices(t: np.ndarray, peak_t_s: float, hold_duration_s: float,
                         hold_margin_lo_s: float = 2.0, nmax: int = 300) -> np.ndarray:
    """Replicate the detector's hold window so paired samples line up exactly."""
    lo = peak_t_s + hold_margin_lo_s
    hi = lo + hold_duration_s
    idx = np.flatnonzero((t >= lo) & (t <= hi))
    if idx.size == 0:
        return idx
    stride = max(1, int(np.ceil(idx.size / max(1, int(nmax)))))
    keep = idx[::stride]
    if keep.size == 0 or keep[-1] != idx[-1]:
        keep = np.append(keep, idx[-1])
    return keep


def normalized_relaxation_holds(inputs: PairedInputs, n_peaks_max: int = 5,
                                min_hold_s: float = 20.0) -> list[dict]:
    """Recover step-hold segments from the FORCE trace, with paired displacement.

    The hold detector is the sibling analyzer's ``detect_holds_in_series``
    (read-only), applied to the signed force in its expected convention
    (tension-negative). Each returned hold carries the normalized force shape
    AND the paired strain/displacement at the same sample indices:

      ``x_s``, ``y_norm``, ``peak_force_gf``, ``peak_strain``,
      ``strain_hold_drift_rel`` (plus the displacement drift and diagnostics).

    If the record is already a single hold (the detector finds no step peaks), the
    post-peak span is returned as one hold instead.
    """
    t = inputs.time_s
    f = inputs.force_gf
    length = inputs.reference_length_mm
    raw = (inputs.displacement_mm / length if length is not None else None)
    quant = _force_quantization_gf(f)
    signed = -f  # tension-positive input -> detector's T = -f convention
    detected = m.detect_holds_in_series(
        t, signed, quant, specimen=inputs.specimen, donor=inputs.donor,
        nmax=300, min_hold_s=float(min_hold_s))

    out: list[dict] = []
    if detected:
        for h in detected[:int(n_peaks_max)]:
            tpk = float(h["peak_t_s"])
            duration = float(h["hold_duration_s"])
            idx = _hold_window_indices(t, tpk, duration)
            pidx = int(np.argmin(np.abs(t - tpk)))
            peak_strain = float(raw[pidx]) if raw is not None else None
            strain_drift_rel = None
            disp_drift_rel = None
            if idx.size >= 2:
                dh = inputs.displacement_mm[idx]
                if dh[0] != 0.0:
                    disp_drift_rel = float((dh[-1] - dh[0]) / abs(dh[0]))
                if raw is not None and raw[idx][0] != 0.0:
                    strain_drift_rel = float((raw[idx][-1] - raw[idx][0])
                                             / abs(raw[idx][0]))
            out.append({
                "x_s": np.asarray(h["x"], dtype=float),
                "y_norm": np.asarray(h["y_norm"], dtype=float),
                "peak_force_gf": float(h["peak_force_gf"]),
                "peak_strain": peak_strain,
                "strain_hold_drift_rel": strain_drift_rel,
                "displacement_hold_drift_rel": disp_drift_rel,
                "peak_t_s": tpk,
                "hold_duration_s": duration,
                "n_points": int(np.asarray(h["x"]).size),
                "single_hold_fallback": False,
            })
        return out

    # single-hold fallback: no step peaks were detected, treat post-peak as the hold
    if t.size < 2:
        return []
    pidx = int(np.argmax(np.abs(f)))
    x = t - float(t[pidx])
    sel = x >= 0.0
    if np.count_nonzero(sel) < 2:
        return []
    xs = x[sel]
    fs = f[sel]
    peak = float(np.max(np.abs(fs)))
    if peak == 0.0:
        peak = 1.0
    rawsel = raw[sel] if raw is not None else None
    peak_strain = float(raw[pidx]) if raw is not None else None
    strain_drift_rel = None
    disp_drift_rel = None
    dsel = inputs.displacement_mm[sel]
    if dsel[0] != 0.0:
        disp_drift_rel = float((dsel[-1] - dsel[0]) / abs(dsel[0]))
    if rawsel is not None and rawsel[0] != 0.0:
        strain_drift_rel = float((rawsel[-1] - rawsel[0]) / abs(rawsel[0]))
    out.append({
        "x_s": xs,
        "y_norm": fs / peak,
        "peak_force_gf": peak,
        "peak_strain": peak_strain,
        "strain_hold_drift_rel": strain_drift_rel,
        "displacement_hold_drift_rel": disp_drift_rel,
        "peak_t_s": float(t[pidx]),
        "hold_duration_s": float(t[-1] - t[pidx]),
        "n_points": int(xs.size),
        "single_hold_fallback": True,
    })
    return out


# ------------------------------------------------------------- parameter gate ----
def absolute_parameter_gate(inputs: PairedInputs) -> dict:
    """ABSTAIN gate for absolute material parameters from a paired record.

    Explicitly maps the four prerequisites of
    ``m.assess_geometry_prerequisites`` from the record:

      * ``area_mm2``                 -> ``inputs.area_mm2`` (or None);
      * ``reference_length_mm``      -> ``inputs.reference_length_mm`` (or None);
      * ``strain_protocol_verified`` -> explicit ``provenance`` flag if supplied,
        else derived True only when both a reference length and a finite
        displacement channel exist, else None;
      * ``displacement_channel``     -> the observed ``inputs.displacement_mm``
        channel (or None).

    ``PASS`` here means the prerequisites EXIST; it is never a material-parameter
    identification.
    """
    area = inputs.area_mm2
    ref = inputs.reference_length_mm
    disp = inputs.displacement_mm
    prov = inputs.provenance if isinstance(inputs.provenance, dict) else {}
    if "strain_protocol_verified" in prov:
        spv: bool | None = bool(prov["strain_protocol_verified"])
    elif ref is not None and disp.size and np.all(np.isfinite(disp)):
        spv = True
    else:
        spv = None
    verdict = dict(m.assess_geometry_prerequisites(
        area_mm2=area,
        reference_length_mm=ref,
        strain_protocol_verified=spv,
        displacement_channel=disp if disp.size else None,
    ))
    verdict.update({
        "inputs_specimen": inputs.specimen,
        "inputs_donor": inputs.donor,
        "inputs_area_mm2": area,
        "inputs_reference_length_mm": ref,
        "inputs_strain_protocol_verified": spv,
        "inputs_displacement_channel": ("present" if disp.size else None),
    })
    return verdict


def apparent_modulus(inputs: PairedInputs) -> dict:
    """Apparent (specimen-level) modulus proxy, only when geometry is present.

    Requires cross-sectional area and reference length. Returns the peak stress
    ``force_N / area_m2`` and an equilibrium stress estimate (median of the last
    10% of the post-peak hold). No value is extrapolated beyond the observed
    strain. Without geometry it returns ``ABSTAIN`` and no number.
    """
    missing: list[str] = []
    if inputs.area_mm2 is None:
        missing.append("area_mm2")
    if inputs.reference_length_mm is None:
        missing.append("reference_length_mm")
    if missing:
        return {
            "state": ABSTAIN,
            "missing": missing,
            "meaning": ("absolute stress/modulus requires geometry "
                        "(cross-sectional area + reference length); force alone "
                        "is not enough."),
        }
    area_m2 = float(mm2_to_m2(inputs.area_mm2))
    strain = strain_from_displacement(inputs.displacement_mm,
                                      inputs.reference_length_mm)
    stress_pa = gf_to_n(inputs.force_gf) / area_m2
    peak_idx = int(np.argmax(strain))
    n_hold = int(strain.size - peak_idx)
    n_tail = max(1, int(np.ceil(0.1 * n_hold)))
    peak_stress = float(np.max(stress_pa))
    eq_stress = float(np.median(stress_pa[-n_tail:]))
    return {
        "state": PASS,
        "peak_stress_Pa": peak_stress,
        "equilibrium_stress_Pa": eq_stress,
        "equilibrium_fraction_of_peak": (eq_stress / peak_stress
                                         if peak_stress != 0.0 else None),
        "peak_strain": float(strain[peak_idx]),
        "strain_end": float(strain[-1]),
        "area_mm2": inputs.area_mm2,
        "area_m2": area_m2,
        "reference_length_mm": inputs.reference_length_mm,
        "specimen": inputs.specimen,
        "donor": inputs.donor,
        "provenance": inputs.provenance,
        "n_points": int(strain.size),
        "n_hold_tail": n_tail,
        "provenance_note": ("peak = max(force_N/area_m2); equilibrium = median of "
                            "the last 10% of the post-peak hold; no extrapolation "
                            "beyond observed strain."),
    }


# --------------------------------------------------------------- tau from pairs ----
def _sls_stress_from_strain(t: np.ndarray, eps: np.ndarray, e_inf: float,
                            e_branch: float, tau: float) -> np.ndarray:
    """Exact single-tau SLS stress for an arbitrary observed strain history.

    ``sigma = E_inf*eps + q``, ``q' = E_branch*eps' - q/tau``. The per-step
    recurrence is the exact closed-form update (known linear-viscoelastic
    mathematics), so no ODE error is introduced when the strain path is
    piecewise linear. Starts from ``q = 0`` at the first sample.
    """
    n = eps.size
    sig = np.empty(n, dtype=float)
    sig[0] = e_inf * eps[0]
    q = 0.0
    for i in range(n - 1):
        dt = float(t[i + 1] - t[i])
        if dt <= 0.0:
            sig[i + 1] = e_inf * eps[i + 1] + q
            continue
        v = (eps[i + 1] - eps[i]) / dt
        a = -np.expm1(-dt / tau)
        ss = e_branch * tau * v
        b = q - ss
        q = q * (1.0 - a) + ss * a
        sig[i + 1] = e_inf * eps[i + 1] + q
    return sig


def _hold_start_index(eps: np.ndarray) -> int:
    """Index where the fast loading ramp ends and the (possibly drifting) hold begins.

    The ramp is where |d eps/dt| is at least a tenth of its maximum; the hold
    begins just after the last such sample. Using ``argmax(eps)`` here would put
    the hold start at the final sample for a positively drifting hold and collapse
    the time-constant guess.
    """
    d = np.abs(np.diff(eps))
    if d.size and float(np.max(d)) > 0.0:
        big = np.flatnonzero(d >= 0.1 * float(np.max(d)))
        if big.size:
            return int(min(int(big[-1]) + 1, eps.size - 1))
    return int(np.argmax(np.abs(eps)))


def _auto_p0(t: np.ndarray, eps: np.ndarray, stress: np.ndarray) -> np.ndarray:
    """Deterministic starting guess for (E_inf, E_branch, tau)."""
    eps_pk = float(np.max(np.abs(eps)))
    sig_pk = float(np.max(np.abs(stress)))
    hold_idx = _hold_start_index(eps)
    n_tail = max(1, int(np.ceil(0.1 * max(1, stress.size - hold_idx))))
    sig_eq = float(np.median(stress[-n_tail:]))
    eps_eq = float(np.abs(eps[-1])) if abs(eps[-1]) > 1e-15 else eps_pk
    e_inf0 = max(abs(sig_eq) / eps_eq, 1e3) if eps_eq != 0.0 else 1e6
    e_branch0 = max((sig_pk - abs(sig_eq)) / eps_pk, 1e3) if eps_pk != 0.0 else 1e6
    tau0 = float(np.clip(0.3 * (t[-1] - t[hold_idx]), 1e-3, 1e5))
    return np.array([e_inf0, e_branch0, tau0], dtype=float)


def _p0_candidates(t: np.ndarray, eps: np.ndarray, stress: np.ndarray) -> list[np.ndarray]:
    """Deterministic multi-start around the automatic guess.

    A positive hold-strain drift can put ``curve_fit`` in a shallow local minimum;
    a small deterministic set of time-constant scales (and stiffer/softer branch
    splits) makes the recovery robust without introducing randomness.
    """
    base = _auto_p0(t, eps, stress)
    span = max(float(t[-1] - t[0]), 1e-9)
    cands = [base]
    for frac in (0.1, 0.25, 0.5):
        c = base.copy()
        c[2] = max(frac * span, 1e-3)
        cands.append(c)
    c = base.copy()
    c[0] *= 0.5
    c[1] *= 2.0
    cands.append(c)
    c = base.copy()
    c[0] *= 1.5
    c[1] *= 0.5
    cands.append(c)
    return cands


def step_validity(t: np.ndarray, eps: np.ndarray, stress: np.ndarray) -> dict:
    """Pre-fit check that a record is a loading-then-hold relaxation step.

    Works on time, engineering strain and stress (any consistent units; every
    metric is a ratio, so it is invariant to L0/area/unit scale). Returns
    ``{"ok", "reasons", "metrics"}``. Reasons (each a REJECT cause):

    * ``too_few_samples``            fewer than ``GUARD_MIN_SAMPLES`` samples;
    * ``time_not_strictly_increasing`` a duplicated or reversed time stamp;
    * ``no_tensile_loading``          strain or stress never rises above its
      first sample (nothing was loaded, so there is nothing to relax);
    * ``unloading_after_peak``        strain falls by more than
      ``GUARD_MAX_RETRACTION_FRAC`` of its increment after the peak (unloading /
      retraction step; a relaxation hold keeps strain near its peak);
    * ``no_hold``                     fewer than ``GUARD_MIN_HOLD_SAMPLES``
      samples after strain first reaches ``GUARD_HOLD_START_FRAC`` of its peak
      increment, or that segment lasts less than ``GUARD_MIN_HOLD_OVER_LOADING``
      of the time spent reaching it (a pure loading or failure ramp: the hold is
      the measurement of a relaxation test and is never a short tail of the ramp);
    * ``force_collapse_below_step_start`` final stress lies more than
      ``|GUARD_MIN_END_FORCE_FRAC|`` of the step's stress increment below the
      step-start stress: the specimen carries clearly less load at a higher strain
      than before the step (loss of load-bearing capacity, not relaxation of
      this step).
    """
    t = np.asarray(t, dtype=float)
    eps = np.asarray(eps, dtype=float)
    stress = np.asarray(stress, dtype=float)
    n = int(t.size)
    reasons: list[str] = []
    metrics: dict[str, Any] = {"n_samples": n}
    if n < GUARD_MIN_SAMPLES:
        reasons.append("too_few_samples")
        return {"ok": False, "reasons": reasons, "metrics": metrics}
    dts = np.diff(t)
    metrics["min_dt_s"] = float(np.min(dts))
    metrics["median_dt_s"] = float(np.median(dts))
    if not np.all(dts > 0.0):
        reasons.append("time_not_strictly_increasing")
    pk = int(np.argmax(eps))
    eps_inc = float(eps[pk] - eps[0])
    sig_inc = float(np.max(stress) - stress[0])
    metrics["strain_increment"] = eps_inc
    metrics["stress_increment"] = sig_inc
    if not (eps_inc > 0.0 and sig_inc > 0.0):
        reasons.append("no_tensile_loading")
        return {"ok": False, "reasons": reasons, "metrics": metrics}
    retraction = float((eps[pk] - np.min(eps[pk:])) / eps_inc)
    metrics["retraction_frac"] = retraction
    if retraction > GUARD_MAX_RETRACTION_FRAC:
        reasons.append("unloading_after_peak")
    h = hold_start_index_frac(eps)
    metrics["hold_start_index"] = h
    metrics["hold_start_s"] = float(t[h])
    metrics["hold_duration_s"] = float(t[-1] - t[h])
    metrics["hold_samples"] = int(n - 1 - h)
    loading_s = float(t[h] - t[0])
    metrics["loading_duration_s"] = loading_s
    metrics["hold_over_loading"] = (float((t[-1] - t[h]) / loading_s)
                                    if loading_s > 0.0 else None)
    if (n - 1 - h < GUARD_MIN_HOLD_SAMPLES
            or (loading_s > 0.0
                and (t[-1] - t[h]) < GUARD_MIN_HOLD_OVER_LOADING * loading_s)):
        reasons.append("no_hold")
    end_frac = float((stress[-1] - stress[0]) / sig_inc)
    metrics["end_force_frac"] = end_frac
    if end_frac < GUARD_MIN_END_FORCE_FRAC:
        reasons.append("force_collapse_below_step_start")
    return {"ok": not reasons, "reasons": reasons, "metrics": metrics}


def hold_start_index_frac(eps: np.ndarray, frac: float = GUARD_HOLD_START_FRAC) -> int:
    """First index where strain reaches ``frac`` of its peak increment.

    Used by the guards only (the fit's starting guess keeps
    :func:`_hold_start_index`). The 95 % crossing is robust to actuator jitter
    during the hold (which can move the 10 %-of-max-slope rule far into the hold;
    one rat anchor step put it at 232 s) and to a positive drifting hold (whose
    peak is the last sample).
    """
    eps = np.asarray(eps, dtype=float)
    inc = float(np.max(eps) - eps[0])
    if inc <= 0.0:
        return int(eps.size - 1)
    return int(np.argmax(eps - eps[0] >= frac * inc))


def tau_identifiability(tau: float, ci: list | None, rmse: float,
                        stress_increment: float, median_dt_s: float,
                        hold_duration_s: float) -> list[str]:
    """Post-fit reasons why a fitted tau is not an identified time constant.

    Each reason is an ABSTAIN cause:

    * ``tau_at_fit_bound``            tau within 0.1 % of a curve_fit bound;
    * ``tau_ci_not_finite``           covariance not available/finite;
    * ``tau_ci_too_wide``             95 % CI width / tau > ``GUARD_MAX_CI_REL_WIDTH``;
    * ``single_tau_sls_misfit``       rmse / stress increment > ``GUARD_MAX_RMSE_REL``;
    * ``tau_below_sampling_resolution`` tau < ``GUARD_MIN_TAU_OVER_DT`` x median dt;
    * ``hold_shorter_than_two_tau``   hold < ``GUARD_MIN_HOLD_OVER_TAU`` x tau
      (the relaxed fraction is then too small to separate tau from E_inf).
    """
    reasons: list[str] = []
    lo, hi = GUARD_TAU_BOUNDS
    if tau <= lo * 1.001 or tau >= hi * 0.999:
        reasons.append("tau_at_fit_bound")
    if ci is None or not all(np.isfinite(ci)):
        reasons.append("tau_ci_not_finite")
    elif tau > 0.0 and (ci[1] - ci[0]) / tau > GUARD_MAX_CI_REL_WIDTH:
        reasons.append("tau_ci_too_wide")
    if not (stress_increment > 0.0) or rmse / stress_increment > GUARD_MAX_RMSE_REL:
        reasons.append("single_tau_sls_misfit")
    if tau < GUARD_MIN_TAU_OVER_DT * median_dt_s:
        reasons.append("tau_below_sampling_resolution")
    if hold_duration_s < GUARD_MIN_HOLD_OVER_TAU * tau:
        reasons.append("hold_shorter_than_two_tau")
    return reasons


def fit_tau_from_pairs(inputs: PairedInputs, p0: Any = None,
                       guards: bool = True) -> dict:
    """Fit a single-tau SLS to the paired record and return ``tau`` with a 95% CI.

    Requires geometry (area + reference length). The EXACT observed strain
    history is propagated through the standard SLS
    (``scipy.optimize.curve_fit`` over ``E_inf``, ``E_branch``, ``tau``), so a
    drifting hold does not masquerade as material complexity: the drift is a
    known input. Without geometry it returns ``ABSTAIN`` (a physical time
    constant is not identifiable from force alone).

    Step-validity guards (``guards=True``, the default):

    * :func:`step_validity` runs before the fit. A record that is not a
      loading-then-hold relaxation step (unloading, no hold, force collapse,
      bad time base, no loading) returns ``state == "REJECT"`` with
      ``reasons`` and no fit.
    * :func:`tau_identifiability` runs after the fit. A tau that is not
      identified returns ``state == "ABSTAIN"`` with ``reasons``. The fitted
      numbers are kept only under ``unaccepted_fit``, never as ``tau_s``.

    ``guards=False`` restores the unguarded behaviour for diagnostics only.
    """
    missing: list[str] = []
    if inputs.area_mm2 is None:
        missing.append("area_mm2")
    if inputs.reference_length_mm is None:
        missing.append("reference_length_mm")
    if missing:
        return {
            "state": ABSTAIN,
            "missing": missing,
            "meaning": ("tau (a physical relaxation time) is not identifiable from "
                        "force alone; paired displacement + geometry are required."),
        }
    from scipy.optimize import curve_fit

    t = inputs.time_s
    eps = strain_from_displacement(inputs.displacement_mm, inputs.reference_length_mm)
    area_m2 = float(mm2_to_m2(inputs.area_mm2))
    stress = gf_to_n(inputs.force_gf) / area_m2

    validity = step_validity(t, eps, stress) if guards else None
    if validity is not None and not validity["ok"]:
        return {
            "state": REJECT,
            "reasons": list(validity["reasons"]),
            "reason": "; ".join(validity["reasons"]),
            "guard_metrics": validity["metrics"],
            "guard_thresholds": dict(GUARD_THRESHOLDS),
            "specimen": inputs.specimen,
            "donor": inputs.donor,
            "n_points": int(t.size),
            "meaning": ("not a loading-then-hold relaxation step; no tau is fitted "
                        "or reported."),
        }

    if p0 is not None:
        p0v = np.asarray(p0, dtype=float).ravel()
        if p0v.size != 3:
            raise ValueError("p0 must supply three values: (E_inf, E_branch, tau)")
        starts = [p0v]
    else:
        starts = _p0_candidates(t, eps, stress)

    def model(tt: np.ndarray, e_inf: float, e_branch: float, tau: float) -> np.ndarray:
        return _sls_stress_from_strain(tt, eps, e_inf, e_branch, tau)

    best: tuple[float, np.ndarray, np.ndarray] | None = None
    last_exc: str | None = None
    for start in starts:
        try:
            popt, pcov = curve_fit(
                model, t, stress, p0=start,
                bounds=([0.0, 0.0, 1e-9], [np.inf, np.inf, 1e7]),
                maxfev=40000)
        except Exception as exc:  # noqa: BLE001 - try the next start
            last_exc = str(exc)
            continue
        rmse_i = float(np.sqrt(np.mean((stress - model(t, *popt)) ** 2)))
        if best is None or rmse_i < best[0]:
            best = (rmse_i, popt, pcov)
    if best is None:
        return {
            "state": FAIL,
            "reason": last_exc or "curve_fit failed for every starting guess",
            "p0_attempts": [[float(x) for x in s] for s in starts],
            "n_points": int(t.size),
        }

    rmse, popt, pcov = best
    tau = float(popt[2])
    ci: list[float] | None = None
    ci_rel: float | None = None
    if np.all(np.isfinite(pcov)):
        perr = np.sqrt(np.diag(pcov))
        if np.all(np.isfinite(perr)):
            half = 1.96 * float(perr[2])
            ci = [tau - half, tau + half]
            ci_rel = float((ci[1] - ci[0]) / tau) if tau != 0.0 else None
    if validity is not None:
        vm = validity["metrics"]
        id_reasons = tau_identifiability(
            tau, ci, rmse, vm["stress_increment"], vm["median_dt_s"],
            vm["hold_duration_s"])
        vm = dict(vm)
        vm.update({
            "ci_rel_width": ci_rel,
            "rmse_rel": (float(rmse / vm["stress_increment"])
                         if vm["stress_increment"] > 0.0 else None),
            "tau_over_dt": float(tau / vm["median_dt_s"]),
            "hold_over_tau": (float(vm["hold_duration_s"] / tau)
                              if tau > 0.0 else None),
        })
        if id_reasons:
            return {
                "state": ABSTAIN,
                "reasons": id_reasons,
                "reason": "; ".join(id_reasons),
                "guard_metrics": vm,
                "guard_thresholds": dict(GUARD_THRESHOLDS),
                "unaccepted_fit": {
                    "tau_s": tau, "tau_ci95_s": ci, "tau_ci95_rel_width": ci_rel,
                    "e_inf_Pa": float(popt[0]), "e_branch_Pa": float(popt[1]),
                    "rmse_Pa": rmse,
                },
                "specimen": inputs.specimen,
                "donor": inputs.donor,
                "n_points": int(t.size),
                "meaning": ("a relaxation step, but the single-tau SLS time constant "
                            "is not identified; the fitted numbers are diagnostic only."),
            }
    else:
        vm = None
    return {
        "state": PASS,
        "guard_metrics": vm,
        "tau_s": tau,
        "tau_ci95_s": ci,
        "tau_ci95_rel_width": ci_rel,
        "e_inf_Pa": float(popt[0]),
        "e_branch_Pa": float(popt[1]),
        "rmse_Pa": rmse,
        "area_m2": area_m2,
        "reference_length_mm": inputs.reference_length_mm,
        "specimen": inputs.specimen,
        "donor": inputs.donor,
        "n_points": int(t.size),
        "model": ("single-tau SLS with the exact observed strain history "
                  "propagated (E_inf, E_branch, tau)"),
        "provenance_note": ("tau is recovered from paired force+displacement; a "
                            "drifting hold is observed input, not a confound."),
    }


# ------------------------------------------------------------------ self-test ----
def self_test() -> dict:
    """Hermetic self-test. No I/O, no writes."""
    out: dict[str, bool] = {}

    # exact unit conversions
    out["gf_to_n_exact_factor"] = bool(
        float(gf_to_n(1.0)) == m.GF_TO_N == 9.80665e-3
        and abs(float(gf_to_n(1000.0)) - 9.80665) < 1e-12)
    out["mm_and_mm2_conversions_exact"] = bool(
        float(mm_to_m(1.0)) == 1e-3 and float(mm2_to_m2(1.0)) == 1e-6
        and float(mm_to_m(250.0)) == 0.25)

    # strain formula + refusal
    s = strain_from_displacement([1.0, 2.0, 5.0], 1000.0)
    out["strain_formula"] = bool(np.allclose(s, [1e-3, 2e-3, 5e-3], rtol=0, atol=1e-15))
    try:
        strain_from_displacement([1.0], 0.0)
        out["strain_refuses_nonpositive_length"] = False
    except ValueError:
        out["strain_refuses_nonpositive_length"] = True
    try:
        strain_from_displacement([1.0], float("nan"))
        out["strain_refuses_nonfinite_length"] = False
    except ValueError:
        out["strain_refuses_nonfinite_length"] = True

    # constant-displacement detector: constant pass, drift fail
    t = np.linspace(0.0, 100.0, 201)
    const = np.full_like(t, 5.0)
    v_const = verify_constant_displacement(const, t, 10.0, 90.0)
    out["constant_detector_true_for_constant"] = bool(
        v_const["classification"] == "constant" and v_const["constant"] is True)
    drift = 5.0 * (1.0 + 0.02 * (t - 10.0) / 80.0)
    v_drift = verify_constant_displacement(drift, t, 10.0, 90.0)
    out["constant_detector_false_for_drift"] = bool(
        v_drift["classification"] == "creep" and v_drift["constant"] is False)
    nonmono = 5.0 + 0.2 * np.sin(np.linspace(0.0, 6.0 * np.pi, t.size))
    v_nm = verify_constant_displacement(nonmono, t, 10.0, 90.0)
    out["constant_detector_nonmonotonic"] = bool(
        v_nm["classification"] == "nonmonotonic")

    # absolute gate ABSTAIN without geometry
    t_syn, eps_syn, stress_syn = _synthetic_ramp_hold()
    area_mm2 = 57.5
    ref_mm = 200.0
    force_gf = stress_syn * (mm2_to_m2(area_mm2)) / m.GF_TO_N
    disp_mm = eps_syn * ref_mm
    no_geo = PairedInputs(time_s=t_syn, force_gf=force_gf, displacement_mm=disp_mm,
                          specimen="syn", donor="d")
    gate = absolute_parameter_gate(no_geo)
    # the displacement channel is present (it is the whole point of paired
    # inputs); the three geometry/protocol prerequisites are missing -> ABSTAIN.
    out["absolute_gate_abstains_without_geometry"] = bool(
        gate["state"] == ABSTAIN
        and {"area_mm2", "reference_length_mm"} <= set(gate["missing"])
        and "displacement_channel" not in gate["missing"])

    with_geo = PairedInputs(time_s=t_syn, force_gf=force_gf, displacement_mm=disp_mm,
                            reference_length_mm=ref_mm, area_mm2=area_mm2,
                            specimen="syn", donor="d")
    gate2 = absolute_parameter_gate(with_geo)
    out["absolute_gate_passes_with_all_four"] = bool(gate2["state"] == PASS)

    # tau recovery on clean paired data
    fit = fit_tau_from_pairs(with_geo)
    out["fit_tau_recovers_known_sls"] = bool(
        fit["state"] == PASS
        and abs(fit["tau_s"] - _SELF_TEST_TAU) / _SELF_TEST_TAU < 0.05)
    out["fit_tau_abstains_without_geometry"] = bool(
        fit_tau_from_pairs(no_geo)["state"] == ABSTAIN)

    # guard: an unloading step is rejected, not fitted
    eps_un = np.where(t_syn > 200.0, eps_syn * np.clip((300.0 - t_syn) / 100.0, 0, 1),
                      eps_syn)
    stress_un = _sls_stress_from_strain(t_syn, eps_un, _SELF_TEST_E_INF,
                                        _SELF_TEST_E_BRANCH, _SELF_TEST_TAU)
    unload = PairedInputs(time_s=t_syn,
                          force_gf=stress_un * mm2_to_m2(area_mm2) / m.GF_TO_N,
                          displacement_mm=eps_un * ref_mm,
                          reference_length_mm=ref_mm, area_mm2=area_mm2)
    fu = fit_tau_from_pairs(unload)
    out["fit_tau_rejects_unloading_step"] = bool(
        fu["state"] == REJECT and "unloading_after_peak" in fu["reasons"]
        and "tau_s" not in fu)
    return out


# Known synthetic truth used only by self_test and its helper.
_SELF_TEST_TAU = 100.0
_SELF_TEST_E_INF = 0.8e9
_SELF_TEST_E_BRANCH = 0.4e9


def _synthetic_ramp_hold(ramp_s: float = 0.5, end_s: float = 300.0,
                         n: int = 601, eps0: float = 0.02) -> tuple:
    """Clean ramp-then-constant-hold strain protocol + exact SLS stress."""
    t = np.linspace(0.0, end_s, n)
    eps = np.where(t <= ramp_s, eps0 * t / ramp_s, eps0)
    stress = _sls_stress_from_strain(t, eps, _SELF_TEST_E_INF,
                                     _SELF_TEST_E_BRANCH, _SELF_TEST_TAU)
    return t, eps, stress


if __name__ == "__main__":  # pragma: no cover - manual probe
    import json
    print(json.dumps(self_test(), indent=2))
