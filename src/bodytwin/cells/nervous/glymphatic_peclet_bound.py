"""Peclet bound for convective vs diffusive solute transport in brain parenchyma.

PRE-REGISTERED CLAIM (C):  For a SMALL solute (TMA+, Rh ~ 0.35 nm) the parenchymal
Peclet number Pe = u_ECS * L / D_eff is < 0.05 under EVERY physiologically-anchored
interstitial-flow scenario, INCLUDING the steelman scenario in which 100% of total CSF
production is forced to transit the parenchyma.  ~C would be Pe >= 0.05 for any such
scenario.  Threshold 0.05 fixed BEFORE running (5% convective contribution = the level
below which convection cannot be the dominant clearance route for that solute).

SECOND PRE-REGISTERED CLAIM (C2, size-crossover): the hydrodynamic radius Rh* at which
Pe = 1 exceeds the measured ECS gap half-width (~20 nm, gaps ~40 nm, Nicholson &
Hrabetova 2017 PMID28755756) under the CENTRAL flow scenario -- i.e. any species large
enough for convection to dominate is sterically hindered/excluded by the ECS itself.
~C2 = Rh* < 20 nm under the central scenario.

THIRD PRE-REGISTERED CLAIM (C3, sleep): the Xie2013-measured sleep ECS expansion
(alpha 0.14 -> 0.23) LOWERS Pe at fixed volumetric ISF flow (u_ECS = u_sup/alpha), by
>= 30%.  ~C3 = Pe ratio >= 1.0.  This makes "sleep enhances clearance" mechanism-neutral.

ALL INPUTS ARE EXTERNAL AND CITED (see PARAMS).  No fitted parameters.
EXTERNAL ANCHOR (not used as an input, held out): Hladky & Barrand 2019 PMID31299992
independently derive, by a completely different route (perivenous-outflow clearance
~1 uL/min/g divided by Ray et al's 40 cm^2/g interfacial area -> 0.25 um/min superficial,
vs 50 um/min needed for diffusion-parity for TMA+), a convection:diffusion ratio of
1/200 = 0.005 for TMA+.  This script must land within a factor of 3 of 0.005 in its
matched scenario or the computation is considered UNANCHORED (machine-checked below).

Reads: nothing (all inputs embedded in PARAMS).
Writes: glymphatic_peclet_evidence.json.
Gate: verdicts C (Pe < 0.05 for TMA+), C2 (crossover radius above the ECS half-gap),
C3 (sleep expansion lowers Pe by >= 30%), plus the held-out Hladky & Barrand anchor check.
Run:  python3 glymphatic_peclet_bound.py           (prints + writes evidence)
      python3 glymphatic_peclet_bound.py --selftest (unit tests, exit 1 on fail)
"""

from __future__ import annotations

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import argparse
import json
import math
import os
import sys

# ---------------------------------------------------------------- constants
K_B = 1.380649e-23  # J/K (SI defined)
T_BODY = 310.15     # K (37 C)
ETA_37 = 0.6913e-3  # Pa.s, water viscosity at 37 C (CRC)

PARAMS = {
    "alpha_awake": (0.14, "Xie et al 2013 Science PMID24136970 (awake ECS volume fraction)"),
    "alpha_sleep": (0.23, "Xie et al 2013 Science PMID24136970 (sleep/anesth ECS volume fraction; +~60%)"),
    "alpha_nominal": (0.20, "Nicholson & Hrabetova 2017 Biophys J PMID28755756 (ECS ~20% of tissue)"),
    "lambda_small": (1.581, "Nicholson & Hrabetova 2017: D_eff = 2/5 D_free for small molecules -> lambda = sqrt(2.5)"),
    "lambda_large": (2.3, "Thorne/Nicholson-class measurements: lambda rises with size for MDa dextrans (sensitivity 1.6-2.5 swept)"),
    "ecs_gap_nm": (40.0, "Nicholson & Hrabetova 2017 PMID28755756: ECS gaps predominantly ~40 nm"),
    "L_um": (200.0, "peri-arterial to peri-venous path length (inter-vessel spacing), mouse/rat cortex; swept 100-500"),
}

# Volumetric interstitial-flow scenarios, expressed as Q/V in uL per gram tissue per min.
# Tissue density taken as 1 g/cm^3 so uL/g/min == uL/cm^3/min.
FLOW_SCENARIOS = {
    # Cserr-lineage classical ISF bulk-flow production rate
    "cserr_low": (0.10, "classical ISF bulk flow 0.1-0.3 uL/g/min (Cserr; reviewed Abbott 2004)"),
    "cserr_high": (0.30, "classical ISF bulk flow upper end"),
    # Hladky & Barrand 2019's glymphatic-favourable estimate from perivenous outflow
    "hladky_perivenous": (1.0, "Hladky & Barrand 2019 PMID31299992: reported clearances ~1 uL/min/g"),
    # STEELMAN / void-floor: ALL mouse CSF production forced through parenchyma
    "steelman_all_csf": (0.325 / 0.4, "mouse CSF production 0.325 uL/min over 0.4 g brain = 0.81 uL/g/min, 100% forced transparenchymal"),
    # Ray, Iliff & Heys 2019 PMID30836968 UPPER LIMIT (not an estimate) -- handled separately
}

# Solutes: (Rh in nm, label, citation)
SOLUTES = [
    (0.35, "TMA+ (74 Da, RTI probe)", "Nicholson RTI standard; D_free 1.3e-5 cm2/s at 37C"),
    (0.5, "mannitol (182 Da)", "Stokes-Einstein"),
    (1.8, "amyloid-beta 1-42 monomer (4.5 kDa)", "Stokes-Einstein, Rh ~1.8 nm"),
    (2.0, "10 kDa dextran", "Smith et al 2017 PMID28826498 / Hadjiev & Amsden 2015"),
    (3.5, "albumin (66 kDa)", "Stokes-Einstein, Rh ~3.5 nm"),
    (5.0, "70 kDa dextran", "Smith et al 2017 PMID28826498"),
    (12.0, "2000 kDa dextran", "Smith et al 2017 PMID28826498 (largest ECS-mobile tracer used)"),
]


# ---------------------------------------------------------------- physics
def stokes_einstein_D(rh_nm: float, T: float = T_BODY, eta: float = ETA_37) -> float:
    """Free-solution diffusion coefficient [cm^2/s] from hydrodynamic radius [nm]."""
    if rh_nm <= 0:
        raise ValueError("rh_nm must be > 0")
    d_si = K_B * T / (6.0 * math.pi * eta * rh_nm * 1e-9)  # m^2/s
    return d_si * 1e4  # cm^2/s


def d_eff(d_free: float, lam: float) -> float:
    """Effective ECS diffusion coefficient: D_eff = D_free / lambda^2."""
    if lam < 1.0:
        raise ValueError("tortuosity lambda must be >= 1")
    return d_free / (lam * lam)


def superficial_velocity(q_ul_per_g_per_min: float, L_um: float) -> float:
    """Slab geometry: fluid produced in a tissue volume must cross path length L.
    u_sup [cm/s] = (Q/V)[1/s] * L[cm].  Q/V in uL/(cm^3 min) -> 1/s."""
    qv_per_s = q_ul_per_g_per_min * 1e-3 / 60.0  # cm^3/(cm^3 s) = 1/s
    return qv_per_s * (L_um * 1e-4)


def peclet(u_ecs_cm_s: float, L_um: float, d_eff_cm2_s: float) -> float:
    return u_ecs_cm_s * (L_um * 1e-4) / d_eff_cm2_s


def lam_for_rh(rh_nm: float) -> float:
    """Tortuosity model: 1.581 for small solutes, rising to 2.3 for Rh >= 12 nm.
    Linear interpolation in Rh between 0.5 and 12 nm.  Larger lambda => larger Pe,
    i.e. this model is GENEROUS to the convection hypothesis for big solutes."""
    lo_r, hi_r = 0.5, 12.0
    lo_l, hi_l = PARAMS["lambda_small"][0], PARAMS["lambda_large"][0]
    if rh_nm <= lo_r:
        return lo_l
    if rh_nm >= hi_r:
        return hi_l
    f = (rh_nm - lo_r) / (hi_r - lo_r)
    return lo_l + f * (hi_l - lo_l)


def pe_for(rh_nm: float, q: float, L_um: float, alpha: float, lam: float | None = None) -> float:
    lam = lam_for_rh(rh_nm) if lam is None else lam
    de = d_eff(stokes_einstein_D(rh_nm), lam)
    u_sup = superficial_velocity(q, L_um)
    return peclet(u_sup / alpha, L_um, de)


def crossover_rh(q: float, L_um: float, alpha: float, target_pe: float = 1.0) -> float:
    """Bisect for the Rh [nm] at which Pe == target_pe."""
    lo, hi = 0.1, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if pe_for(mid, q, L_um, alpha) < target_pe:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


# ---------------------------------------------------------------- selftest
def selftest() -> int:
    fails = []

    def chk(name, cond, detail=""):
        if not cond:
            fails.append(f"{name}: {detail}")

    # 1. Stokes-Einstein must reproduce the textbook small-ion D within 2x.
    d_tma = stokes_einstein_D(0.35)
    chk("SE_TMA", 8e-6 < d_tma < 2.5e-5, f"D(TMA+)={d_tma:.3e} cm2/s vs literature 1.3e-5")
    # 2. Nicholson 2/5 rule must be recovered by the tortuosity used.
    ratio = d_eff(1.0, PARAMS["lambda_small"][0])
    chk("NICHOLSON_2_5", abs(ratio - 0.4) < 0.005, f"D_eff/D_free={ratio:.4f} vs 0.400")
    # 3. Pe must scale as L^2 (u_sup ~ L, Pe ~ u*L).
    p1 = pe_for(0.35, 1.0, 100.0, 0.20)
    p2 = pe_for(0.35, 1.0, 200.0, 0.20)
    chk("PE_L2", abs(p2 / p1 - 4.0) < 1e-6, f"ratio={p2/p1:.4f} vs 4.0")
    # 4. Pe must be linear in Q.
    chk("PE_LINQ", abs(pe_for(0.35, 2.0, 200.0, 0.2) / pe_for(0.35, 1.0, 200.0, 0.2) - 2.0) < 1e-9, "")
    # 5. Pe must rise with solute size (convection's relative advantage grows with size).
    seq = [pe_for(r, 1.0, 200.0, 0.20) for r, _, _ in SOLUTES]
    chk("PE_MONOTONE_SIZE", all(b > a for a, b in zip(seq, seq[1:])), f"{seq}")
    # 6. Crossover bisection must invert pe_for.
    rstar = crossover_rh(1.0, 200.0, 0.20)
    chk("CROSSOVER_INVERTS", abs(pe_for(rstar, 1.0, 200.0, 0.20) - 1.0) < 1e-3, f"Pe(Rh*)={pe_for(rstar,1.0,200.0,0.20):.4f}")
    # 7. EXTERNAL ANCHOR (held out, never used as input): Hladky & Barrand 2019 give
    #    convection:diffusion = 1/200 = 0.005 for TMA+ at their ~1 uL/min/g clearance.
    pe_matched = pe_for(0.35, 1.0, 200.0, 0.20)
    chk("EXTERNAL_ANCHOR_HLADKY", 0.005 / 3.0 <= pe_matched <= 0.005 * 3.0,
        f"Pe(TMA+, 1 uL/g/min, L=200um)={pe_matched:.5f} vs external 0.005 (factor-3 band)")
    # 8. Guard against a degenerate/void result: a 200x flow must break the bound.
    chk("VOID_FLOOR_BREAKS", pe_for(0.35, 200.0, 200.0, 0.20) > 1.0,
        "at Ray/Iliff/Heys upper-limit flow the bound MUST fail, else the test is vacuous")
    # 9. Sanity: alpha increase lowers Pe at fixed Q.
    chk("SLEEP_LOWERS_PE", pe_for(0.35, 1.0, 200.0, 0.23) < pe_for(0.35, 1.0, 200.0, 0.14), "")

    for f in fails:
        print("SELFTEST FAIL " + f)
    print(f"SELFTEST {'PASS' if not fails else 'FAIL'} ({9-len(fails)}/9)")
    return 1 if fails else 0


# ---------------------------------------------------------------- main
def run() -> dict:
    L = PARAMS["L_um"][0]
    alpha = PARAMS["alpha_nominal"][0]
    out = {"params": {k: {"value": v[0], "source": v[1]} for k, v in PARAMS.items()},
           "constants": {"k_B": K_B, "T_K": T_BODY, "eta_Pa_s": ETA_37},
           "scenarios": {}, "sweeps": {}, "verdicts": {}}

    for sname, (q, cite) in FLOW_SCENARIOS.items():
        u_sup = superficial_velocity(q, L)
        rows = []
        for rh, label, scite in SOLUTES:
            lam = lam_for_rh(rh)
            df = stokes_einstein_D(rh)
            rows.append({"solute": label, "Rh_nm": rh, "lambda": round(lam, 3),
                         "D_free_cm2_s": df, "D_eff_cm2_s": d_eff(df, lam),
                         "Pe": pe_for(rh, q, L, alpha), "source": scite})
        out["scenarios"][sname] = {
            "Q_uL_per_g_per_min": q, "source": cite,
            "u_superficial_um_per_min": u_sup * 1e4 * 60.0,
            "u_ECS_um_per_min": (u_sup / alpha) * 1e4 * 60.0,
            "solutes": rows,
            "crossover_Rh_nm_Pe1": crossover_rh(q, L, alpha),
        }

    # Ray/Iliff/Heys 2019 UPPER LIMIT expressed as a superficial velocity, not a Q.
    u_ray = 50.0 / 60.0 * 1e-4  # 50 um/min -> cm/s
    out["scenarios"]["ray_iliff_heys_UPPER_LIMIT"] = {
        "u_superficial_um_per_min": 50.0,
        "source": "Ray, Iliff & Heys 2019 PMID30836968: simulations give an UPPER LIMIT ~50 um/min superficial",
        "note": "this is an upper limit stated by the glymphatic-camp modellers themselves, not an estimate",
        "solutes": [{"solute": lab, "Rh_nm": rh,
                     "Pe": peclet(u_ray / alpha, L, d_eff(stokes_einstein_D(rh), lam_for_rh(rh)))}
                    for rh, lab, _ in SOLUTES],
    }

    # Sensitivity sweep on L and lambda for the small solute + the largest tracer.
    out["sweeps"]["L_um"] = {
        str(int(Lx)): {"Pe_TMA_hladky": pe_for(0.35, 1.0, Lx, alpha),
                       "Pe_2MDa_hladky": pe_for(12.0, 1.0, Lx, alpha)}
        for Lx in (100.0, 200.0, 350.0, 500.0)}
    out["sweeps"]["lambda_large_solute"] = {
        str(lm): pe_for(12.0, 1.0, L, alpha, lam=lm) for lm in (1.6, 2.0, 2.3, 2.5)}

    # Sleep leg (C3): fixed Q, alpha awake -> sleep.
    pe_wake = pe_for(0.35, 1.0, L, PARAMS["alpha_awake"][0])
    pe_sleep = pe_for(0.35, 1.0, L, PARAMS["alpha_sleep"][0])
    out["sleep_leg"] = {"Pe_awake": pe_wake, "Pe_sleep": pe_sleep,
                        "Pe_ratio_sleep_over_wake": pe_sleep / pe_wake,
                        "alpha_awake": PARAMS["alpha_awake"][0], "alpha_sleep": PARAMS["alpha_sleep"][0],
                        "note": "at fixed volumetric ISF flow, ECS expansion lowers interstitial velocity as 1/alpha and raises D_eff; both lower Pe"}

    # PERIVASCULAR-SPACE CONTRAST LEG: same Pe formalism, measured PVS velocity.
    # Mestre et al 2018 Nat Commun PMID30451853: particle-tracking velocimetry, n=7 mice,
    # typical PVS flow speed 18.7 um/s (95% CI 9.4-28).  PVS is an open fluid channel:
    # no ECS tortuosity (lambda=1) and no alpha division.
    v_pvs = 18.7e-4  # cm/s
    out["pvs_contrast"] = {
        "v_PVS_um_per_s": 18.7, "ci95_um_per_s": [9.4, 28.0], "n_mice": 7,
        "source": "Mestre et al 2018 Nat Commun 9:4878 PMID30451853 (particle tracking velocimetry)",
        "ratio_to_parenchymal_u_ECS_hladky": 18.7 * 60.0 / out["scenarios"]["hladky_perivenous"]["u_ECS_um_per_min"],
        "Pe_PVS": {lab: {str(int(Lp)): peclet(v_pvs, Lp, stokes_einstein_D(rh))
                         for Lp in (200.0, 1000.0)} for rh, lab, _ in SOLUTES},
    }

    # ---- pre-registered verdicts (machine-evaluated)
    small_pes = {s: v["solutes"][0]["Pe"] for s, v in out["scenarios"].items() if s != "ray_iliff_heys_UPPER_LIMIT"}
    out["verdicts"]["C_small_solute_Pe_lt_0.05"] = {
        "threshold": 0.05, "per_scenario": small_pes,
        "max": max(small_pes.values()),
        "PASS": max(small_pes.values()) < 0.05}
    rstar = out["scenarios"]["hladky_perivenous"]["crossover_Rh_nm_Pe1"]
    out["verdicts"]["C2_crossover_exceeds_ECS_half_gap"] = {
        "Rh_star_nm": rstar, "ECS_half_gap_nm": PARAMS["ecs_gap_nm"][0] / 2.0,
        "PASS": rstar > PARAMS["ecs_gap_nm"][0] / 2.0}
    out["verdicts"]["C3_sleep_lowers_Pe_by_ge_30pct"] = {
        "ratio": pe_sleep / pe_wake, "PASS": (pe_sleep / pe_wake) <= 0.70}
    # honest counter-verdict: is the LARGE-solute claim bounded too?
    pe_large = out["scenarios"]["hladky_perivenous"]["solutes"][-1]["Pe"]
    out["verdicts"]["LARGE_SOLUTE_NOT_BOUNDED"] = {
        "Pe_2MDa_dextran": pe_large,
        "OPEN": pe_large > 0.1,
        "note": "Pe of order 0.1-1 for the largest ECS-mobile tracer: the bound does NOT close the large-solute case"}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--out", default=os.path.join(OUT_ROOT, "glymphatic_peclet_bound",
                                                  "glymphatic_peclet_evidence.json"))
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if selftest() != 0:
        print("ABORT: selftest failed")
        return 1
    res = run()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(os.path.abspath(a.out), "w") as fh:
        json.dump(res, fh, indent=1)
    for s, v in res["scenarios"].items():
        print(f"\n[{s}]  u_sup={v.get('u_superficial_um_per_min'):.4g} um/min")
        for r in v["solutes"]:
            print(f"   Pe={r['Pe']:.4g}  {r['solute']}")
        if "crossover_Rh_nm_Pe1" in v:
            print(f"   Rh*(Pe=1) = {v['crossover_Rh_nm_Pe1']:.3g} nm")
    print("\nVERDICTS:", json.dumps(res["verdicts"], indent=1))
    print("SLEEP:", json.dumps(res["sleep_leg"], indent=1))
    print("\nwrote", os.path.abspath(a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
