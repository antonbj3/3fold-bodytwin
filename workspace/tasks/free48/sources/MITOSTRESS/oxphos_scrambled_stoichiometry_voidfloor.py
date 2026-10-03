#!/usr/bin/env python3
"""
Scrambled-stoichiometry void-floor adversary for MODEL-OXIDATIVE-PHOSPHORYLATION.
Re-executes the model's own two mechanisms (pure stoichiometric arithmetic, and the
dynamical ODE steady-state solve) but with deliberately WRONG proton-per-electron /
proton-per-ATP integers, and asks: does the model STILL converge to something that
looks like the textbook P/O (~2.3-2.5 NADH), i.e. is "landing near the measured value"
a property of the arithmetic/dynamics structure rather than evidence the specific
stoichiometry (CI=4,CIII=4,CIV=2,c_ring=8) is correct?

PRE-REGISTERED falsifier: if >=20% of a WIDE grid of deliberately-wrong-but-textbook-
scale integer stoichiometries ALSO land within the model's own tolerance band of the
measured P/O (gate g08: |P/O-2.5|<=0.10, and the looser textbook band [2.2,2.8]), the
"convergence" claim is void-floor-consistent (hollow) -- the arithmetic's own bounded
small-integer structure does the work, not the specific stoichiometry.
"""
import sys
import json
import itertools

import os
from pathlib import Path

# PENDING_INDEPENDENT_REVIEW -- path portability only; no model change.
_HERE = Path(__file__).resolve().parent
_DEP = "oxphosswing_9f2c_reimpl_steadystate.py"
_DEP_DEFAULT = "source_repository/data/body_twin/agent_scratch_preserved"


def _dep_dir():
    """Locate the reused reimpl module: env override, then beside this file, then
    its canonical location. Lets the cell run from any working directory."""
    cands = []
    env = os.environ.get("BODYTWIN_SCRATCH_DIR")
    if env:
        cands.append(Path(env))
    cands += [_HERE, Path(_DEP_DEFAULT)]
    for c in cands:
        if (c / _DEP).exists():
            return c
    raise FileNotFoundError(f"{_DEP} not found; set BODYTWIN_SCRATCH_DIR")


OUT_DIR = Path(os.environ.get("CELL_OUT_DIR", str(_HERE)))
OUT_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(_dep_dir()))
import oxphosswing_9f2c_reimpl_steadystate as ox  # noqa: E402

MEASURED_NADH = 2.5
GATE_TOL = 0.10          # the model's own g08 tolerance
LOOSE_TEXTBOOK_BAND = (2.2, 2.8)

# ============================================================================
# PART 1 -- pure stoichiometric-arithmetic model (mitochondrial_oxphos.py logic).
# P/O(NADH) = (CI + CIII + CIV) / (c_ring/3 + 1)
# Real model uses CI in {3,4}, CIII=4 fixed, CIV=2 fixed, c_ring in {8,10}.
# Adversary: sweep a WIDE, deliberately-scrambled grid across the same "textbook small
# integer" scale (no absurd 1000x values -- that would be an unfair adversary; these are
# all plausible-LOOKING integers a careless stoichiometry table could contain).
# ============================================================================
CI_RANGE = range(1, 8)      # real: 3 or 4
CIII_RANGE = range(1, 8)    # real: 4 (Q-cycle)
CIV_RANGE = range(1, 6)     # real: 2
CRING_RANGE = range(6, 16)  # real: 8 (mammal) or 10 (yeast); known biological c-ring span is 8-15

grid_results = []
for ci, ciii, civ, cring in itertools.product(CI_RANGE, CIII_RANGE, CIV_RANGE, CRING_RANGE):
    h2e = ci + ciii + civ
    hatp = cring / 3.0 + 1.0
    po = h2e / hatp
    grid_results.append((ci, ciii, civ, cring, h2e, hatp, po))

total_n = len(grid_results)
n_within_gate_tol = sum(1 for r in grid_results if abs(r[6] - MEASURED_NADH) <= GATE_TOL)
n_within_loose_band = sum(1 for r in grid_results
                           if LOOSE_TEXTBOOK_BAND[0] <= r[6] <= LOOSE_TEXTBOOK_BAND[1])
frac_gate_tol = n_within_gate_tol / total_n
frac_loose_band = n_within_loose_band / total_n

# Real-stoichiometry-only subgrid membership check (sanity: real values must be IN the swept grid)
real_ci_values = {3, 4}
real_hits = [r for r in grid_results if r[0] in real_ci_values and r[1] == 4 and r[2] == 2
             and r[3] in (8, 10)]

# ============================================================================
# PART 2 -- dynamical ODE steady-state model (the KINETIC mechanism, not pure ratio
# arithmetic). Scramble nH_ETC (real=10 NADH path) and nH_ATP (real=3.667-4.333) to
# deliberately wrong integers/values spanning the same order of magnitude, and ask:
# (a) does a steady-state root still EXIST (the convergence claim itself)?
# (b) what P/O = J_syn/J_ETC results at state3 (ADP=1.5mM)?
# ============================================================================
def po_at_state3(nh_etc, nh_atp, adp=1.5):
    dg0 = ox.P["dG0_ETC_NADH"]
    roots, _ = ox.solve_steady(adp, nH_ETC=nh_etc, dG0_ETC=dg0, nH_ATP=nh_atp)
    cands = [r for r in roots if 0 < r < 300]
    if not cands:
        return None, None
    dpsi = max(cands)
    _, nadh, jetc, jsyn, jleak = ox.charge_residual(dpsi, adp, nh_etc, dg0, nH_ATP=nh_atp)
    if jetc <= 0:
        return dpsi, None
    return dpsi, jsyn / jetc


# Deliberately WRONG (scrambled) stoichiometry set -- same textbook integer scale, swapped/
# permuted/replaced values that do NOT correspond to any real biological reading.
SCRAMBLED_ODE_CASES = [
    ("real_nH_ETC10_nH_ATP3.667", 10, 3.6667),
    ("scrambled_swap_nH_ETC=4_nH_ATP=10", 4, 10.0),          # swap the two numbers' roles
    ("scrambled_nH_ETC=6_nH_ATP=6", 6, 6.0),                 # arbitrary equal wrong pair
    ("scrambled_nH_ETC=14_nH_ATP=2", 14, 2.0),               # implausibly high H+/2e-, low H+/ATP
    ("scrambled_nH_ETC=7_nH_ATP=5", 7, 5.0),                 # mid-range wrong pair
    ("scrambled_nH_ETC=3_nH_ATP=1", 3, 1.0),                 # both far too low
    ("scrambled_nH_ETC=20_nH_ATP=8", 20, 8.0),               # both far too high, ratio~2.5 by design check
]

ode_results = []
for label, nh_etc, nh_atp in SCRAMBLED_ODE_CASES:
    dpsi, po = po_at_state3(nh_etc, nh_atp)
    converged = dpsi is not None
    ode_results.append(dict(label=label, nH_ETC=nh_etc, nH_ATP=nh_atp,
                             converged=converged, dPsi=dpsi, PO=po,
                             within_gate_tol=(po is not None and abs(po - MEASURED_NADH) <= GATE_TOL),
                             within_loose_band=(po is not None and LOOSE_TEXTBOOK_BAND[0] <= po <= LOOSE_TEXTBOOK_BAND[1])))

n_ode_scrambled = len(SCRAMBLED_ODE_CASES) - 1  # exclude the "real" row
n_ode_scrambled_converged = sum(1 for r in ode_results[1:] if r["converged"])
n_ode_scrambled_within_loose = sum(1 for r in ode_results[1:] if r["within_loose_band"])

out = dict(
    part1_arithmetic_grid=dict(
        total_n=total_n,
        n_within_gate_tol_0p10=n_within_gate_tol,
        frac_within_gate_tol=frac_gate_tol,
        n_within_loose_band_2p2_2p8=n_within_loose_band,
        frac_within_loose_band=frac_loose_band,
        n_real_stoichiometry_combos_in_grid=len(real_hits),
        note="P/O(NADH)=(CI+CIII+CIV)/(c_ring/3+1) swept over CI in [1,7], CIII in [1,7], "
             "CIV in [1,5], c_ring in [6,15] -- textbook-scale-plausible but including many "
             "biologically-wrong combinations.",
    ),
    part2_ode_scrambled=ode_results,
    part2_summary=dict(
        n_scrambled_cases=n_ode_scrambled,
        n_scrambled_converged=n_ode_scrambled_converged,
        n_scrambled_within_loose_band=n_ode_scrambled_within_loose,
        frac_scrambled_converged=n_ode_scrambled_converged / n_ode_scrambled,
        frac_scrambled_within_loose_band=n_ode_scrambled_within_loose / n_ode_scrambled,
    ),
)

print(json.dumps(out, indent=2, default=str))
with open(OUT_DIR / "voidfloor_oxphos_results.json", "w") as f:
    json.dump(out, f, indent=2, default=str)
