#!/usr/bin/env python3
"""CLOUD-W4-A194-ARMS synthetic moment-arm audit.

Compares three moment-arm estimators on original synthetic muscle paths:
  VW : virtual work, r = -dL/dq, with dL/dq computed exactly by forward-mode
       automatic differentiation (dual numbers) through the path-length code.
  FL : force line, r = axis . ((Q - j) x u), the moment of a unit tension
       applied at the distal attachment Q along the path's last segment u.
       Uses only positions, no derivative of L.
  FD : central finite difference, r = -(L(q+h) - L(q-h)) / (2h).

No external anatomy data is used. Units: metres and radians in code; results
are reported in mm. Run: python3 arms_audit.py  (writes results_synthetic.json)
"""
import json
import math
import sys

import numpy as np

# ---- pre-declared audit settings (fixed before running) ----
H_FD = 1e-5                      # rad, central-difference step
SWITCH_EXCLUSION = math.radians(0.5)  # rad, band around switches excluded from PASS test
REL_TOL = 1e-3                   # 0.1 % (PREREG synthetic criterion)
ABS_FLOOR = 1e-6                 # m (1 um); rel. error uses max(|r|, floor)
Q_GRID = np.radians(np.linspace(-120.0, 120.0, 2401))  # 0.1 deg grid


# ---------------- minimal forward-mode AD ----------------
class D:
    __slots__ = ("v", "d")

    def __init__(self, v, d=0.0):
        self.v, self.d = float(v), float(d)

    def __add__(s, o):
        o = o if isinstance(o, D) else D(o)
        return D(s.v + o.v, s.d + o.d)
    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, D) else D(o)
        return D(s.v - o.v, s.d - o.d)

    def __rsub__(s, o):
        return D(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, D) else D(o)
        return D(s.v * o.v, s.d * o.v + s.v * o.d)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, D) else D(o)
        return D(s.v / o.v, (s.d * o.v - s.v * o.d) / (o.v * o.v))

    def __rtruediv__(s, o):
        return D(o) / s

    def __neg__(s):
        return D(-s.v, -s.d)


def val(x):
    return x.v if isinstance(x, D) else float(x)


def dsqrt(x):
    x = x if isinstance(x, D) else D(x)
    r = math.sqrt(x.v)
    return D(r, x.d / (2 * r))


def dsin(x):
    x = x if isinstance(x, D) else D(x)
    return D(math.sin(x.v), math.cos(x.v) * x.d)


def dcos(x):
    x = x if isinstance(x, D) else D(x)
    return D(math.cos(x.v), -math.sin(x.v) * x.d)


def dacos(x):
    x = x if isinstance(x, D) else D(x)
    return D(math.acos(x.v), -x.d / math.sqrt(1 - x.v * x.v))


def datan2(y, x):
    y = y if isinstance(y, D) else D(y)
    x = x if isinstance(x, D) else D(x)
    den = x.v * x.v + y.v * y.v
    return D(math.atan2(y.v, x.v), (x.v * y.d - y.v * x.d) / den)


def norm(vec):
    return dsqrt(sum(c * c for c in vec))


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


# ---------------- kinematics ----------------
def rot_z(q, p):
    c, s = dcos(q), dsin(q)
    return [c * p[0] - s * p[1], s * p[0] + c * p[1], D(0.0) + p[2]]


def rodrigues(axis, q, p):
    """Rotate p about unit axis through origin by q (all AD-compatible)."""
    k = axis
    c, s = dcos(q), dsin(q)
    kxp = [k[1] * p[2] - k[2] * p[1], k[2] * p[0] - k[0] * p[2], k[0] * p[1] - k[1] * p[0]]
    kdp = k[0] * p[0] + k[1] * p[1] + k[2] * p[2]
    return [p[i] * c + kxp[i] * s + k[i] * kdp * (1 - c) for i in range(3)]


# ---------------- path models ----------------
# Each model(q) returns (L, Q, A, axis, joint, switch_state):
#   Q = distal attachment (moving with q), A = point Q pulls toward.
def straight_planar(P, Q0):
    def model(q):
        Q = rot_z(q, Q0)
        return norm(sub(Q, P)), Q, P, (0, 0, 1), (0, 0, 0), 0
    return model


def via_distal(P, V0, Q0):
    """Via point V fixed on distal body; segment V-Q lies in distal body."""
    def model(q):
        V, Q = rot_z(q, V0), rot_z(q, Q0)
        L = norm(sub(V, P)) + norm(sub(Q, V))
        # only P-V crosses the joint: force acts on distal body at V toward P
        return L, V, P, (0, 0, 1), (0, 0, 0), 0
    return model


def cylinder_wrap(P, Q0, c, R, side):
    """Planar geodesic around a circle (cylinder) fixed on the proximal body.
    side=+1: path passes the cylinder counter-clockwise from P to Q."""
    def model(q):
        Q = rot_z(q, Q0)
        cx, cy = c
        dPx, dPy = P[0] - cx, P[1] - cy
        dQx, dQy = Q[0] - cx, Q[1] - cy
        dP, dQ = dsqrt(dPx * dPx + dPy * dPy), dsqrt(dQx * dQx + dQy * dQy)
        phP, phQ = datan2(dPy, dPx), datan2(dQy, dQx)
        aP, aQ = dacos(R / dP), dacos(R / dQ)
        tP = phP + side * aP          # angle of tangent point from P
        tQ = phQ - side * aQ          # angle of tangent point from Q
        th = side * (tQ - tP)
        k = math.floor((val(th) + math.pi) / (2 * math.pi))
        th = th - 2 * math.pi * k     # wrap to (-pi, pi]
        if val(th) > 0.0:             # wrapping active
            L = dsqrt(dP * dP - R * R) + R * th + dsqrt(dQ * dQ - R * R)
            TQ = [cx + R * dcos(tQ), cy + R * dsin(tQ), D(0.0)]
            return L, Q, TQ, (0, 0, 1), (0, 0, 0), 1
        return norm(sub(Q, P)), Q, P, (0, 0, 1), (0, 0, 0), 0
    return model


def straight_oblique3d(P, Q0, axis, joint):
    ax = np.asarray(axis, float)
    ax = ax / np.linalg.norm(ax)
    j = np.asarray(joint, float)

    def model(q):
        Qr = rodrigues(list(ax), q, list(np.asarray(Q0) - j))
        Q = [Qr[i] + j[i] for i in range(3)]
        return norm(sub(Q, P)), Q, P, tuple(ax), tuple(j), 0
    return model


def conditional_via(P, V, Q0, qlo, qhi):
    """Via point on proximal body active only for qlo <= q <= qhi
    (OpenSim-style conditional path point): L and r jump at the switches."""
    def model(q):
        Q = rot_z(q, Q0)
        if qlo <= val(q) <= qhi:
            return norm(sub(V, P)) + norm(sub(Q, V)), Q, V, (0, 0, 1), (0, 0, 0), 1
        return norm(sub(Q, P)), Q, P, (0, 0, 1), (0, 0, 0), 0
    return model


# ---------------- estimators ----------------
def r_vw(model, q):
    return -model(D(q, 1.0))[0].d


def r_fl(model, q):
    _, Q, A, axis, j, _ = model(D(q))
    Qv = np.array([val(x) for x in Q]) - np.array(j)
    Av = np.array([val(x) for x in A]) - np.array(j)
    u = (Av - Qv) / np.linalg.norm(Av - Qv)
    return float(np.dot(axis, np.cross(Qv, u)))


def r_fd(model, q, h=H_FD):
    return -(val(model(D(q + h))[0]) - val(model(D(q - h))[0])) / (2 * h)


def state(model, q):
    return model(D(q))[5]


def switch_angles(model, grid):
    s = np.array([state(model, q) for q in grid])
    idx = np.nonzero(np.diff(s))[0]
    out = []
    for i in idx:  # bisection to 1e-12 rad
        a, b = grid[i], grid[i + 1]
        sa = state(model, a)
        for _ in range(60):
            m = 0.5 * (a + b)
            if state(model, m) == sa:
                a = m
            else:
                b = m
        out.append(0.5 * (a + b))
    return out


def relerr(x, ref):
    return abs(x - ref) / max(abs(ref), ABS_FLOOR)


# ---------------- synthetic "muscle groups" ----------------
def attachment_sets(base_P, base_Q, rng, n=5, jitter=0.01):
    """Controlled attachment perturbations: base + n-1 seeded jitters (+-1 cm)."""
    sets = [(np.array(base_P, float), np.array(base_Q, float))]
    for _ in range(n - 1):
        sets.append((np.array(base_P) + rng.uniform(-jitter, jitter, 3) * [1, 1, 0],
                     np.array(base_Q) + rng.uniform(-jitter, jitter, 3) * [1, 1, 0]))
    return sets


def build_groups(rng):
    G = {}
    G["G1_straight_planar"] = [straight_planar(list(P), list(Q))
                               for P, Q in attachment_sets([-0.15, 0.03, 0], [0.05, 0.02, 0], rng)]
    G["G2_via_point_distal"] = [via_distal(list(P), [0.03, 0.035, 0.0], list(Q))
                                for P, Q in attachment_sets([-0.15, 0.03, 0], [0.12, 0.01, 0], rng)]
    G["G3_cylinder_centred"] = [cylinder_wrap(list(P), list(Q), (0.0, 0.0), 0.025, -1)
                                for P, Q in attachment_sets([-0.20, 0.035, 0], [0.18, 0.030, 0], rng)]
    G["G4_cylinder_offset"] = [cylinder_wrap(list(P), list(Q), (-0.01, 0.012), 0.018, -1)
                               for P, Q in attachment_sets([-0.20, 0.04, 0], [0.18, 0.030, 0], rng)]
    G["G5_oblique_3d"] = [straight_oblique3d(list(P + [0, 0, 0.02]), list(Q + [0, 0, -0.01]),
                                             (0.2, -0.3, 1.0), (0.004, -0.003, 0.01))
                          for P, Q in attachment_sets([-0.15, 0.04, 0], [0.06, 0.02, 0], rng)]
    return G


def audit_model(model, grid, closed_form=None):
    sw = switch_angles(model, grid)
    rows = []
    for q in grid:
        near = any(abs(q - s) < SWITCH_EXCLUSION for s in sw)
        vw, fl, fd = r_vw(model, q), r_fl(model, q), r_fd(model, q)
        rows.append((q, near, vw, fl, fd, closed_form(q) if closed_form else None))
    away = [r for r in rows if not r[1]]
    nearr = [r for r in rows if r[1]]
    res = {
        "n_angles": len(rows), "n_away": len(away), "switches_deg": [math.degrees(s) for s in sw],
        "r_range_mm": [1e3 * min(r[2] for r in rows), 1e3 * max(r[2] for r in rows)],
        "max_rel_vw_vs_fl_away": max(relerr(r[2], r[3]) for r in away),
        "max_rel_vw_vs_fd_away": max(relerr(r[2], r[4]) for r in away),
        "max_abs_vw_vs_fl_away_mm": 1e3 * max(abs(r[2] - r[3]) for r in away),
        "max_abs_vw_vs_fd_away_mm": 1e3 * max(abs(r[2] - r[4]) for r in away),
    }
    if closed_form:
        res["max_abs_vw_vs_closed_form_away_mm"] = 1e3 * max(
            abs(r[2] - r[5]) for r in away if r[5] is not None)
    if nearr:
        res["max_abs_vw_vs_fd_near_switch_mm"] = 1e3 * max(abs(r[2] - r[4]) for r in nearr)
        res["max_abs_vw_vs_fl_near_switch_mm"] = 1e3 * max(abs(r[2] - r[3]) for r in nearr)
    return res


def fd_step_study(model, q):
    ref = r_vw(model, q)
    return {f"{h:g}": 1e3 * abs(r_fd(model, q, h) - ref) for h in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)}


def attachment_sensitivity(Q0, P, q=math.radians(30), d=1e-3):
    """Change in moment arm (mm) per 1 mm shift of the distal attachment."""
    base = r_vw(straight_planar(P, Q0), q)
    out = {}
    for name, dv in (("x", [d, 0, 0]), ("y", [0, d, 0])):
        Qs = [Q0[i] + dv[i] for i in range(3)]
        out[name] = 1e3 * (r_vw(straight_planar(P, Qs), q) - base)
    return out


def main():
    rng = np.random.default_rng(194)
    groups = build_groups(rng)
    out = {"settings": {"h_fd_rad": H_FD, "switch_exclusion_deg": math.degrees(SWITCH_EXCLUSION),
                        "rel_tol": REL_TOL, "abs_floor_m": ABS_FLOOR,
                        "q_grid_deg": [-120, 120, 2401], "seed": 194,
                        "python": sys.version.split()[0], "numpy": np.__version__},
           "groups": {}}
    for gname, models in groups.items():
        per = []
        for m in models:
            cf = None
            if gname == "G3_cylinder_centred":  # wrapped: r = +R (corrected sign, see RESULTS.md)
                cf = (lambda mm: (lambda q: (0.025 if state(mm, q) == 1 else None)))(m)
            per.append(audit_model(m, Q_GRID, cf))
        worst_rel = max(max(p["max_rel_vw_vs_fl_away"], p["max_rel_vw_vs_fd_away"]) for p in per)
        out["groups"][gname] = {"configs": per, "worst_rel_err_away": worst_rel,
                                "synthetic_pass": worst_rel < REL_TOL}
    # counterexample: conditional via point (discontinuous switch)
    cv = conditional_via([-0.15, 0.03, 0], [-0.01, 0.03, 0], [0.05, 0.02, 0],
                         math.radians(-30), math.radians(40))
    out["counterexample_conditional_via"] = audit_model(cv, Q_GRID)
    qs = math.radians(40) - 0.5 * H_FD
    out["counterexample_conditional_via"]["fd_straddling_switch_mm"] = {
        "q_deg": math.degrees(qs), "vw_mm": 1e3 * r_vw(cv, qs), "fd_mm": 1e3 * r_fd(cv, qs)}
    # FD step study on offset cylinder, away from and at lift-off
    m4 = groups["G4_cylinder_offset"][0]
    sw4 = switch_angles(m4, Q_GRID)
    out["fd_step_study_mm"] = {"G4_q=0deg": fd_step_study(m4, 0.0)}
    if sw4:
        out["fd_step_study_mm"][f"G4_q=switch+1e-4rad({math.degrees(sw4[0]):.3f}deg)"] = \
            fd_step_study(m4, sw4[0] + 1e-4)
    out["attachment_sensitivity_mm_per_mm"] = attachment_sensitivity(
        [0.05, 0.02, 0.0], [-0.15, 0.03, 0.0])
    out["synthetic_method_pass"] = all(g["synthetic_pass"] for g in out["groups"].values())
    out["n_groups_pass"] = sum(g["synthetic_pass"] for g in out["groups"].values())
    with open("results_synthetic.json", "w") as f:
        json.dump(out, f, indent=1)
    for g, v in out["groups"].items():
        print(f"{g}: worst rel err away from switches = {v['worst_rel_err_away']:.3e} "
              f"pass={v['synthetic_pass']} switches={[round(s,3) for p in v['configs'] for s in p['switches_deg']]}")
    print("counterexample:", json.dumps(out["counterexample_conditional_via"], indent=1))
    print("fd step:", json.dumps(out["fd_step_study_mm"], indent=1))
    print("sens:", out["attachment_sensitivity_mm_per_mm"])
    print("SYNTHETIC METHOD PASS:", out["synthetic_method_pass"])


if __name__ == "__main__":
    main()
