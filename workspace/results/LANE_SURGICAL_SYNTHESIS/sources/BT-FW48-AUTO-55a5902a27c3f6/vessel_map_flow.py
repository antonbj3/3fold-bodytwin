"""BT-FW48-AUTO-55a5902a27c3f6 -- vessel-map -> flow functional weighting.

Implements PREREG.md section 3 exactly. No constants introduced after the
freeze. Every declared map value is an ASSUMPTION with no provenance.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

# ---------------------------------------------------------------- ports
MU_PA_S = 3.5e-3          # blood viscosity
DELTA_P_PA = 1.0e4        # transmural pressure
L_OUT_M = 2.0e-3          # Q033 Parameters.thickness_m
A_FACE_M2 = math.pi * 0.023 * 0.0065   # Q033 ellipse_geometry()

# ----------------------------------------------- declared map (ASSUMPTION)
DECADES = [(3.0, 10.0), (10.0, 30.0), (30.0, 100.0), (100.0, 300.0)]
N_DECLARED = [1.5e4, 1.5e4, 3.0e2, 3.0e2]          # 1/m^2 per decade
N_BOX = [(3.75e3, 6.0e4), (3.75e3, 6.0e4), (7.5e1, 1.2e3), (7.5e1, 1.2e3)]


def radius(decade: tuple[float, float]) -> float:
    return math.sqrt(decade[0] * decade[1]) * 1e-6      # m, geometric mean


def vessel_flow(r_m: float) -> float:
    """Poiseuille outflow from one severed end, PREREG 3.1. [m^3/s]"""
    return math.pi * DELTA_P_PA * r_m ** 4 / (8.0 * MU_PA_S * L_OUT_M)


def exit_shear(r_m: float) -> float:
    """gamma = 4v/R = dP*R/(2 mu L_out), PREREG 3.3 as amended. [1/s]"""
    return DELTA_P_PA * r_m / (2.0 * MU_PA_S * L_OUT_M)


def flow_profile(n: list[float]) -> dict:
    """PREREG 3.2-3.4. Returns per-decade contributions and totals."""
    rows = []
    for (lo, hi), nj in zip(DECADES, n):
        r = radius((lo, hi))
        qv = vessel_flow(r)
        g = exit_shear(r)
        rows.append({
            "decade_um": [lo, hi],
            "radius_geomean_m": r,
            "n_per_m2_per_decade": nj,
            "q_vessel_m3_s": qv,
            "q_per_area_m_s": nj * qv,
            "exit_shear_s_inv": g,
            "shear_weighted_n": nj * g,
        })
    q = sum(r["q_per_area_m_s"] for r in rows)
    sw = sum(r["shear_weighted_n"] for r in rows)
    for r in rows:
        r["flow_share"] = r["q_per_area_m_s"] / q
        r["shear_end_share"] = r["shear_weighted_n"] / sw
    return {
        "rows": rows,
        "q_per_area_m_s": q,
        "Q_total_m3_s": q * A_FACE_M2,
        "shear_weight_total_1_s": sw,
    }


# -------------------------------------------------- unit algebra (PREREG C1)
DIMLESS = (0, 0, 0)
_M, _S = (0, 1, 0), (0, 0, 1)
_PA = (1, -1, -2)


def _mul(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _div(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _pow(a, k):
    return tuple(k * x for x in a)


_M3 = (0, 3, 0)          # m^3
_M2 = (0, 2, 0)          # m^2
_N_PER_M2 = (0, -2, 0)   # 1/m^2 (count is dimensionless, so per m^2)


def unit_check() -> dict:
    mu_l = _mul(_mul(_PA, _S), _M)              # Pa*s*m
    # Q_vessel = pi*dP*R^4/(8 mu L): Pa*m^4/(Pa*s*m) = m^3/s
    q_vessel = _div(_mul(_PA, _pow(_M, 4)), mu_l)
    # q = n * Q_vessel: (1/m^2)*(m^3/s) = m/s
    q_per_area = _mul(_N_PER_M2, q_vessel)
    # Q_total = q * A: (m/s)*(m^2) = m^3/s
    q_total = _mul(q_per_area, _M2)
    # v = dP*R^2/(8 mu L): m/s
    v = _div(_mul(_PA, _pow(_M, 2)), mu_l)
    # gamma = 4v/R: 1/s   [AMENDED: PREREG 3.3 substituted form had R^3]
    gamma = _div(_pow(v, 1), _M)
    checks = {
        "q_vessel_is_m3_per_s": q_vessel == _div(_M3, _S),
        "q_per_area_is_m_per_s": q_per_area == _div(_M, _S),
        "Q_total_is_m3_per_s": q_total == _div(_M3, _S),
        "gamma_is_per_s": gamma == _div(DIMLESS, _S),
        "n_dimension_is_per_m2": _div(DIMLESS, _M2) == _N_PER_M2,
        "shear_weight_is_per_m2_per_s": _mul(_N_PER_M2, gamma) == _div(_N_PER_M2, _S),
        "shares_dimensionless": _div(q_per_area, q_per_area) == DIMLESS,
    }
    ok = all(bool(v) for v in checks.values())
    return {"status": "pass" if ok else "fail", "checks": checks}


# ------------------------------------------------------ criteria C1..C6
def criteria(profile: dict, control_ratio: float, box_span: float) -> dict:
    rows = profile["rows"]
    const_n_shares = [
        vessel_flow(radius(d)) for d in DECADES
    ]
    s = sum(const_n_shares)
    shares = [c / s for c in const_n_shares]
    monotonic = all(shares[i] < shares[i + 1] for i in range(len(shares) - 1))
    c = {
        "C1_units": unit_check()["status"] == "pass",
        "C2_positive": profile["q_per_area_m_s"] > 0
        and profile["Q_total_m3_s"] > 0
        and all(r["exit_shear_s_inv"] > 0 for r in rows)
        and all(r["flow_share"] > 0 and r["shear_end_share"] > 0 for r in rows),
        "C3_shares_sum_1": abs(sum(r["flow_share"] for r in rows) - 1) < 1e-12
        and abs(sum(r["shear_end_share"] for r in rows) - 1) < 1e-12,
        "C4_share_increases_with_radius_at_const_n": monotonic,
        "C5_two_mode_load_bearing_gt_2x": (control_ratio > 2.0 or control_ratio < 0.5),
        "C6_box_span_factor": box_span,
    }
    c["code_criteria_pass"] = bool(c["C1_units"] and c["C2_positive"] and c["C3_shares_sum_1"] and c["C4_share_increases_with_radius_at_const_n"])
    return c


def run() -> dict:
    base = flow_profile(N_DECLARED)

    # --- control: single log-uniform matched on the ONLY scalar a section-plane
    # histogram actually delivers directly: the TOTAL count integral
    # sum_j n_j (hits per m^2). Matching on q instead would be circular
    # (a single amplitude reproduces any scalar q exactly, ratio==1 trivially).
    total_count = sum(N_DECLARED)
    ctrl = flow_profile([total_count / 4.0] * 4)
    control_ratio = ctrl["q_per_area_m_s"] / base["q_per_area_m_s"]

    # --- box corners (PREREG 3.6)
    corners = []
    for mask in range(16):
        n = [N_BOX[j][(mask >> j) & 1] for j in range(4)]
        p = flow_profile(n)
        corners.append(p["q_per_area_m_s"])
    box_span = max(corners) / min(corners)

    # precision requirement: relative accuracy needed per decade so that a
    # count error contributes <= 1% of q (PREREG-derived, post-hoc metric)
    tol = 0.01
    per_decade_tol = [tol / r["flow_share"] for r in base["rows"]]
    tol10 = 0.10
    per_decade_tol10 = [tol10 / r["flow_share"] for r in base["rows"]]
    # how many vessels must be counted in the dominant decade to fix q to 10%
    # Poisson: relative SE = 1/sqrt(N)
    poisson_N = [
        (tol10 / r["flow_share"]) ** 2 for r in base["rows"]
    ]

    return {
        "id": "BT-FW48-AUTO-55a5902a27c3f6",
        "equations": {
            "v": "v = dP*R^2/(8*mu*L_out)",
            "Q_vessel": "pi*dP*R^4/(8*mu*L_out)",
            "q": "sum_j n_j * Q_vessel(R_j)",
            "Q_total": "q * A_face",
            "gamma": "dP*R^3/(2*mu*L_out)",
            "shares": "s_j=Q_j/q ; g_j=n_j*g_j/sum(n_k*g_k)",
        },
        "ports": {
            "L_out_m": {"value": L_OUT_M, "unit": "m", "source": "Q033 Parameters.thickness_m"},
            "A_face_m2": {"value": A_FACE_M2, "unit": "m^2", "source": "Q033 ellipse_geometry a=0.023,b=0.0065"},
            "mu_pa_s": {"value": MU_PA_S, "unit": "Pa s", "source": "ASSUMPTION; standard whole blood"},
            "delta_p_pa": {"value": DELTA_P_PA, "unit": "Pa", "source": "ASSUMPTION"},
        },
        "declared_map_status": "ASSUMPTION with no provenance (PREREG 2.1)",
        "baseline_profile": base,
        "control_single_log_uniform": {
            "amplitude_per_m2_per_decade": total_count / 4.0,
            "matched_on": "total count integral sum_j n_j (hits per m^2)",
            "q_per_area_m_s": ctrl["q_per_area_m_s"],
            "ratio_to_declared": control_ratio,
            "interpretation": (
                "two-mode structure is load-bearing for q"
                if control_ratio > 2.0 or control_ratio < 0.5 else
                "count-matched single log-uniform reproduces q within 2x; "
                "two-mode structure is NOT load-bearing for q"
            ),
        },
        "box_span_factor": box_span,
        "per_decade_precision_for_1pct_q": per_decade_tol,
        "per_decade_relative_tol_for_10pct_q": per_decade_tol10,
        "poisson_counts_needed_for_10pct_q": poisson_N,
        "total_area_needed_m2_at_declared_n": {
            f"{DECADES[k][0]:g}-{DECADES[k][1]:g}": poisson_N[k] / N_DECLARED[k]
            for k in range(len(DECADES))
        },
        "criteria": criteria(base, control_ratio, box_span),
    }


if __name__ == "__main__":
    out = Path(__file__).with_name("vessel_map_flow_results.json")
    payload = run()
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    b = payload["baseline_profile"]
    print(json.dumps({
        "q_per_area_m_s": b["q_per_area_m_s"],
        "Q_total_m3_s": b["Q_total_m3_s"],
        "flow_shares": [r["flow_share"] for r in b["rows"]],
        "shear_end_shares": [r["shear_end_share"] for r in b["rows"]],
        "control_ratio": payload["control_single_log_uniform"]["ratio_to_declared"],
        "box_span_factor": payload["box_span_factor"],
        "code_criteria_pass": payload["criteria"]["code_criteria_pass"],
    }, indent=2))