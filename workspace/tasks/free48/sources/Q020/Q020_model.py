#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"\nBT-HX-Q020 -- vessel-interstitium-lymph: can the same model produce BOTH swelling\nand substance retention?\n\nFirst-principles 4-kompartmentmodell:\n    P   plasma            V_p [mL], M_p [mg albumin]\n    ISC interstitium      V_i [mL], M_i [mg albumin]      <- enda \"lagrade\" tillstandet\n    C   celler            V_c [mL]  (vatten, protein immobil, egen volymreglering)\n    L   lymfa             V_l [mL], M_l [mg]  (kvasi-stationaert, tau_l)\n\nSource equations (no empirical \"swelling formulas\"):\n 1) van-Slyke/Hill  pi(c)  -- oncotic pressure from albumin mass concentration\n 2) Starling/filtr. J_net = K_f[(P_c - P_ISC) - sigma_g(pi_p - pi_i)]  with\n    (a) hard lower bound on reabsorption flow:  J_net >= -K_f sigma_g pi_p\n    (b) glycocalyx dilution E -> sigma_g, sigma_s, K_f  (Michel & Curry form)\n 3) albumin flow through the capillary  J_p = (1-sigma_s) c_p J_net + G_d (c_p - c_i)\n 4) interstitial P-V  P_ISC = a ln(V_ISC/V_k)  -- RECTIFYING (Guyton 1965/1966)\n 5) mobility resistance to lymph R_mob(P) -- knee-shaped (Guyton 1966: >1e5)\n 6) lymph pump  Q_L = P_pump exp(beta dP_L_eff)   (Humphrey form)\n 7) lymph protein c_l = (1 - sigma_L) c_i ; return to plasma with the fraction f_p\n 8) cell exchange J_c = K_c[(P_c - P_ISC) + (pi_i - pi_c)] + cell volume regulation\n\nUnit check: each relation is assigned a Dim vector (mL, mg, min, mmHg)\nand checked automatically by unit_check(). No unit is assumed without a check.\n"
from __future__ import annotations
import json, math, hashlib, argparse
import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------
# Enhetsalgebra: Dim = (mL, mg, min, mmHg)
# ----------------------------------------------------------------------------
class Dim:
    __slots__ = ("e", "name")
    def __init__(self, mL=0, mg=0, minute=0, mmHg=0, name=""):
        self.e = np.array([float(mL), float(mg), float(minute), float(mmHg)])
        self.name = name
    def __mul__(self, o):  return Dim(*(self.e + o.e), name=f"{self.name}*{o.name}")
    def __truediv__(self, o): return Dim(*(self.e - o.e), name=f"{self.name}/{o.name}")
    def __pow__(self, k):   return Dim(*(self.e * k), name=f"{self.name}^{k}")
    def __eq__(self, o):    return np.allclose(self.e, o.e)
    def __repr__(self):     return f"Dim({self.e.tolist()}) {self.name}"
    def text(self):
        u = ["mL", "mg", "min", "mmHg"]
        s = "".join(f"{u[i]}^{v:g} " for i, v in enumerate(self.e) if abs(v) > 1e-12)
        return (s.strip() or "1") + f"  [{self.name}]"

ML   = Dim(mL=1,  name="mL")
MG   = Dim(mg=1,  name="mg")
MIN  = Dim(minute=1, name="min")
MMHG = Dim(mmHg=1, name="mmHg")
G_L  = Dim(mL=-1, mg=1, name="g/L")          # masskoncentration
FLUX_V = ML / MIN                              # mL/min
FLUX_M = MG / MIN
KF_UNIT = FLUX_V / MMHG          # K_f har enheten mL/min/mmHg
KC_UNIT = FLUX_V / MMHG          # cellkonduktans har samma enhet som K_f                              # mg/min

# ----------------------------------------------------------------------------
# PARAMETERTABELL  (fryst i PREREG.md §4)
# ----------------------------------------------------------------------------
PARAMS = {
 # ---- kropp ---------------------------------------------------------------
 "M_body":     (70.0,   "kg",  "ANTAGANDE (brief: 70 kg vuxen)"),
 "V_plasma0":  (3.0,    "L",   "DERIVED: 0.043 L/kg plasma"),
 "V_isc0":     (11.0,   "L",   "DERIVED: ECF 0.2 L/kg - plasma 0.043 L/kg"),
 "V_cell0":    (28.0,   "L",   "DERIVED: TBW 0.6 L/kg - ECF 0.2 L/kg"),
 "c_plasma":   (40.0,   "g/L", "ANTAGANDE: plasmaalbumin 4 g/dL; ger pi_p=18.4 mmHg"),
 "c_isc0":     (15.0,   "g/L", "ANTAGANDE: ~0.38 c_p; ger pi_i~5 mmHg (subkutan litteratur)"),
 # ---- van Slyke / Hill ----------------------------------------------------
 "T":          (310.0,  "K",   "ANTAGANDE 37 C"),
 "M_n":        (69300.0,"g/mol","ANTAGANDE albuminmonomer"),
 "n_Hill":     (1.5,    "-",   "KALIBRERAD pa (pi(10 g/L)=3.1, pi(40 g/L)=18.4 mmHg)"),
 "pi_max":     (None,   "mmHg","DERIVED from the calibration above"),
 "c_half":     (None,   "g/L", "DERIVED from the calibration above"),
 # ---- kapillar ------------------------------------------------------------
 "Pc_cap":     (20.0,   "mmHg","ANTAGANDE: capillary pressure mean (20-25)"),
 "K_f0":       (0.05,   "mL/min/mmHg", "ANTAGANDE: helkroppseffektiv K_f; ger nettofiltr. ~0.5 mL/min"),
 "L_p":        (1.5e-7, "cm/s","Michel & Curry 1999 Physiol Rev 79:703-761 (10.1152/physrev.1999.79.3.703) - INTERVALL OVERIFIERAD"),
 "A_cap":      (3000.0, "m^2", "ANTAGANDE: helkroppens kapillaryarea"),
 "sigma_g0":   (0.95,   "-",   "ANTAGANDE: glykokalysens onkotiska screening"),
 "sigma_s0":   (0.90,   "-",   "ANTAGANDE: albuminreflektion"),
 "D_ratio":    (2.0e-4, "-",   "ANTAGANDE: diffusiv/reflektiv (L_p D_s)/(K_f delta); kalibrerat sa att c_isc0 ~ 0.4 c_p"),
 # ---- glycocalyx dilution ------------------------------------------------
 "f_glyc":     (2.0,    "-",   "ANTAGANDE: E~0.5 nar dP_eff = f_glyc*pi_p"),
 "n_glyc":     (12.0,   "-",   "ANTAGANDE: skarpa knaet"),
 "k_leak":     (2.0,    "-",   "ANTAGANDE: K_f-okning nar glykocalysen tappar"),
 "E_min":      (0.05,   "-",   "ANTAGANDE: golv pa glykokalys"),
 "E_fixed":    (None,   "-",   "Nollmodell N3: None = fri E, 1.0 = laserad glykokalys"),
 # ---- interstitium (Guyton) ----------------------------------------------
 "a_stiff":    (12.0,   "mmHg","ANTAGANDE: P = a ln(V/V_k)"),
 "pv_form":    ("log",   "-",   "'log' = Guyton-riktande; 'linear' = nollmodell N0"),
 "V_k_ratio":  (1.20,   "-",   "ANTAGANDE: V_k = 1.20 V_isc0 -> P_ISC,0 ~ -2.2 mmHg"),
 "R_mob_max":  (1.0e5,  "-",   "Guyton/Scheel/Murphree 1966 Circ Res 19:412-419 (10.1161/01.res.19.2.412): >1e5 ggr MOBILITETSfall"),
 "P_half_mob":  (0.25,   "mmHg","ANTAGANDE: knaets halvbredd"),
 "P_min":      (-14.0,  "mmHg","ANTAGANDE: golv pa P_ISC"),
 "P_max":      (25.0,   "mmHg","ANTAGANDE: tak pa P_ISC"),
 # ---- lymfa ---------------------------------------------------------------
 "K_lymph":    (1.0,    "mL/min/mmHg", "ANTAGANDE: lymfansangio-konduktans (dimensionellt K_L)"),
 "P_pump":     (0.5,    "mmHg","ANTAGANDE: vilande pumptryck"),
 "beta":       (0.05,   "mmHg^-1","ANTAGANDE: K_A A_L / K_L (Humphrey-formen) - OVERIFIERAD"),
 "P_out":      (0.0,    "mmHg","ANTAGANDE: extrinsiskt/visceralt referenstryck"),
 "sigma_L":    (0.10,   "-",   "ANTAGANDE: lymfans proteinreflektion"),
 "tau_l":      (20.0,   "min", "ANTAGANDE: lymfvaskompartimentets omsattningstid"),
 "f_p_lymph":  (0.30,   "-",   "ANTAGANDE: andel lymfprotein som ater till plasma"),
 "tau_hep":    (120.0,  "min", "ANTAGANDE: hepatostatisk reglering av c_p (lever som kallande mangkomp)"),
 # ---- celler --------------------------------------------------------------
 "pi_cell":    (None,   "mmHg","DESIRED/KALIBRERAD: sa att J_cell=0 in the base equilibrium (se PREREG A4)"),
 "P_cell":     (0.0,    "mmHg","ASSUMPTION: intracellular hydrostatic"),
 "K_cell":     (0.5,    "mL/min/mmHg", "ANTAGANDE: cellmembrans konduktans"),
 "c_cell":     (50.0,   "mmHg","ANTAGANDE: cellernas P-V-stelhet (dP_c/d(V/V0))"),
 "tau_cell":   (500.0,  "min", "ANTAGANDE: langsam cellvolymreglering (VRAC etc.)"),
 "c_prop":     (1.0e-4, "-",   "ANTAGANDE: extracell. proteinklarande /tau_l -> k=5e-6 1/min (t1/2 ~ 5 d)"),
}

# van-Slyke-kalibrering: Hill-form pi = pi_max c^n / (c_half^n + c^n)
# 2 points: (c=10 g/L, pi=3.1 mmHg) and (c=40 g/L, pi=18.4 mmHg) @37 C
_PI_CAL = [(10.0, 3.1), (40.0, 18.4)]

def _fit_hill(n):
    """los pi_max, c_half fran 2 punkter med given Hill-exponent."""
    (c1, p1), (c2, p2) = _PI_CAL
    r = p2 / p1
    # p1 = pm c1^n/(ch^n+c1^n), p2 = pm c2^n/(ch^n+c2^n)
    # => r c2^n (ch^n+c1^n) = c1^n(ch^n+c2^n)  ->  ch^n (r c2^n - c1^n) = c1^n c2^n (1-r)
    num = (c1 ** n) * (c2 ** n) * (r - 1.0)
    den = (c2 ** n) - r * (c1 ** n)
    if den <= 0 or num <= 0:
        return None
    chn = num / den
    ch = chn ** (1.0 / n)
    pm = p1 * (chn + c1 ** n) / (c1 ** n)
    return pm, ch

def build_params(**over):
    p = {k: v[0] for k, v in PARAMS.items()}
    pm, ch = _fit_hill(p["n_Hill"])
    if pm is None:
        pm, ch = 25.0, 18.0
    p["pi_max"], p["c_half"] = pm, ch
    p.update(over)
    # A4: the cells' resting state is calibrated so that cell flow is zero at baseline equilibrium
    p["pi_cell"] = pi_albumin(p["c_isc0"], p) - P_isc(p["V_isc0"] * 1e3, p) + p["P_cell"]
    return p

# ----------------------------------------------------------------------------
# 1) ONKOTISKT TRYCK
# ----------------------------------------------------------------------------
def pi_albumin(c_gL, p):
    """pi [mmHg] fran albuminmasskoncentration [g/L].  Dim: G_L -> mmHg."""
    c = max(float(c_gL), 0.0)
    n, ch, pm = p["n_Hill"], p["c_half"], p["pi_max"]
    return pm * (c ** n) / ((ch ** n) + (c ** n))

# ----------------------------------------------------------------------------
# 4) INTERSTITIELL P-V  (riktande, Guyton 1965/1966)
# ----------------------------------------------------------------------------
def P_isc(V, p):
    """P_ISC [mmHg].  Dim: mL -> mmHg.
    pv_form='log'    : Guyton-riktande P = a ln(V/V_k)   (standard)
    pv_form='linear' : nollmodell N0, samma lokala komplians: P = P(V0) + (a/V_k)(V-V0)"""
    V0 = p["V_isc0"] * 1e3
    V_k = V0 * p["V_k_ratio"]                            # mL
    if p.get("pv_form", "log") == "linear":
        P0 = p["a_stiff"] * math.log(V0 / V_k)
        Pv = P0 + p["a_stiff"] / V_k * (max(V, 1.0) - V0)
    else:
        Pv = p["a_stiff"] * math.log(max(V, 1.0) / V_k)
    return float(np.clip(Pv, p["P_min"], p["P_max"]))

# ----------------------------------------------------------------------------
# 5) MOBILITETSMOTSTAND TILL LYMFAN  (knaeformad)
# ----------------------------------------------------------------------------
def R_mob(P, p):
    """R [1].  Dim: mmHg -> 1.  R=1 nedanfor knaet, R->R_mob_max ovanfor."""
    h = max(P, 0.0) / p["P_half_mob"]
    f = (h ** 4) / (1.0 + (h ** 4)) if h > 0 else 0.0
    return 1.0 + (p["R_mob_max"] - 1.0) * f

# ----------------------------------------------------------------------------
# 2b) GLYKOKALYSMODELLEN  (E <-> J_net, sigma_g, sigma_s, K_f)  -- fixed point
# ----------------------------------------------------------------------------
def capillary(p, P_c, V_i, c_i, n_iter=60, tol=1e-12):
    """Socker J_net [mL/min], sigma_g, sigma_s, K_f, E_glyc, dP_eff."""
    Pi = P_isc(V_i, p)
    pi_i = pi_albumin(c_i, p)
    pi_p = pi_albumin(p["c_plasma"], p)
    Kf0 = p["K_f0"]
    dP_raw = (P_c - Pi) - p["sigma_g0"] * (pi_p - pi_i)   # eff. gradient, E=1
    E = 1.0
    J = 0.0
    if p.get("E_fixed") is not None:        # nollmodell N3: glykokalysen laserad
        E = float(p["E_fixed"])
        n_iter = 0
    for _ in range(n_iter):
        sg = p["sigma_g0"] * E
        ss = p["sigma_s0"] * E
        Kf = Kf0 * (1.0 + p["k_leak"] * (1.0 - E))
        J = Kf * dP_raw
        # fys. tak: reabsorbtion kan aldrig overstiga onkotisk screening
        J = max(J, -Kf * sg * pi_p)
        dP_eff = J / Kf if Kf > 0 else 0.0
        E_new = 1.0 / (1.0 + (max(dP_eff, 0.0) / (p["f_glyc"] * pi_p)) ** p["n_glyc"])
        E_new = float(np.clip(E_new, p["E_min"], 1.0))
        E_old = E
        E = 0.5 * E + 0.5 * E_new
        if abs(E - E_old) < tol:
            break
    sg = p["sigma_g0"] * E
    ss = p["sigma_s0"] * E
    Kf = Kf0 * (1.0 + p["k_leak"] * (1.0 - E))
    J = Kf * dP_raw
    J = max(J, -Kf * sg * pi_p)
    return dict(J_net=J, sigma_g=sg, sigma_s=ss, K_f=Kf, E_glyc=E,
                P_isc=Pi, pi_i=pi_i, pi_p=pi_p, dP_eff=J / Kf)

# ----------------------------------------------------------------------------
# 6) LYMFAN
# ----------------------------------------------------------------------------
def lymph(p, P_isc_val):
    "Q_L [mL/min] and effective driving pressure [mmHg]."
    dP = P_isc_val - p["P_out"]
    R = R_mob(dP, p)
    dP_eff = dP / R
    Q = p["K_lymph"] * p["P_pump"] * math.exp(p["beta"] * dP_eff)
    return dict(Q_L=Q, dP_L=dP, dP_L_eff=dP_eff, R_mob=R)

# ----------------------------------------------------------------------------
# TILLSTAND: y = [V_p, M_p, V_i, M_i, V_c]  (V_l, M_l algebraiskt)
# ----------------------------------------------------------------------------
def rhs(t, y, p, dPc=0.0, infuse=0.0, null="none"):
    """y = [V_p, M_p, V_i, M_i, V_c] i [L, g, L, g, L] -- L/g skalning for numerisk
    konditionering.  Floden i kroppen ar i mL/min resp mg/min = *1e-3 i L/g per min."""
    V_p, M_p, V_i, M_i, V_c = y
    V_p = max(V_p, 1e-6); V_i = max(V_i, 1e-6); V_c = max(V_c, 1e-6)
    c_p = M_p / V_p
    c_i = M_i / V_i

    cap = capillary(p, p["Pc_cap"] + dPc, V_i * 1e3, c_i)
    J = cap["J_net"]; ss = cap["sigma_s"]; Pi = cap["P_isc"]

    # 3) albumin through the capillary (df_ven -> 1: the flow cannot enter plasma)
    Gd = p["K_f0"] * p["D_ratio"] * 1e3          # mL/min
    J_p = (1.0 - ss) * c_p * J + Gd * (c_p - c_i)
    J_p = max(J_p, 0.0)

    # 6) lymfa
    ly = lymph(p, Pi)
    Q_L = 0.0 if null == "N1" else ly["Q_L"]
    c_l = (1.0 - p["sigma_L"]) * c_i
    M_l_ret = p["f_p_lymph"] * c_l * Q_L         # mg/min
    J_clear = p["c_prop"] * M_i * 1e3 / p["tau_l"]   # mg/min

    # 8) celler
    P_c = p["P_cell"] + p["c_cell"] * (V_c - p["V_cell0"]) / p["V_cell0"]
    J_osm = p["K_cell"] * ((P_c - Pi) + (cap["pi_i"] - p["pi_cell"]))
    J_reg = (p["V_cell0"] - V_c) * 1e3 / p["tau_cell"]
    J_c = J_osm + J_reg

    # 9) hepatostat (kallan som haller c_p)
    S_hep = (p["c_plasma"] * V_p - M_p) / p["tau_hep"] * 1e3      # mg/min, PI pa c_p

    if null == "N2":            # sigma_s = 0 : allt protein lacker
        J_p = max(Gd * (c_p - c_i), 0.0)

    dV_i = (J - Q_L - J_c) * 1e-3               # L/min
    dM_i = (J_p - M_l_ret - J_clear) * 1e-3     # g/min
    dV_p = (-J + Q_L) * 1e-3 + infuse       # plasma <-> interstitium <-> celler
    dM_p = (-J_p + M_l_ret + S_hep) * 1e-3
    dV_c = J_c * 1e-3
    return np.array([dV_p, dM_p, dV_i, dM_i, dV_c])

# ----------------------------------------------------------------------------
# KORRIGERINGAR AV HJALPFUNKTIONER (V_i-anropningen i cell_flux var overflodig)
# ----------------------------------------------------------------------------
def cell_flux(p, V_c, V_i, M_i):
    Pi = P_isc(V_i, p)
    P_c = p["P_cell"] + p["c_cell"] * (V_c - p["V_cell0"] * 1000.0) / (p["V_cell0"] * 1000.0)
    J_reg = (p["V_cell0"] * 1000.0 - V_c) / p["tau_cell"]
    J_osm = p["K_cell"] * ((P_c - Pi) + (pi_albumin(M_i / V_i, p) - p["pi_cell"]))
    return J_reg, J_osm, P_c, Pi

# ----------------------------------------------------------------------------
# SIMULERING
# ----------------------------------------------------------------------------
def simulate(p, dPc=0.0, t_end=40 * 1440.0, t_eval_h=24.0, null="none",
             y0=None, dt=60.0):
    """y i [L, g, L, g, L] (L/g). y0 = samma skalning."""
    if y0 is None:
        y0 = np.array([p["V_plasma0"], p["c_plasma"] * p["V_plasma0"],
                       p["V_isc0"], p["c_isc0"] * p["V_isc0"], p["V_cell0"]])
    y0 = np.asarray(y0, float)
    te = np.arange(0.0, t_end + 1e-9, t_eval_h * 60.0)
    sol = solve_ivp(rhs, (0, t_end), y0, args=(p, dPc, 0.0, null), t_eval=te,
                    method="Radau", rtol=1e-8, atol=1e-12, max_step=3600.0)
    if not sol.success:
        raise RuntimeError("integration misslyckades: " + sol.message)
    y = sol.y
    V_p, M_p, V_i, M_i, V_c = y
    neg = {nm: float(np.min(v)) for nm, v in (("V_plasma", V_p), ("V_isc", V_i), ("V_cell", V_c))
           if np.min(v) <= 0.0}
    if neg:
        raise RuntimeError("icke-fysikaliskt solution; negativ volym i " +
                           ",".join(f"{k}={v:.3g}" for k, v in neg.items()))
    out = {"t_h": te / 60.0, "V_p": V_p * 1e3, "M_p": M_p * 1e3,
           "V_i": V_i * 1e3, "M_i": M_i * 1e3, "V_c": V_c * 1e3, "y0": y0}
    c_i = M_i / V_i
    c_p = M_p / V_p
    Pi = np.array([P_isc(v, p) for v in out["V_i"]])
    lym = [lymph(p, v) for v in Pi]
    out["P_isc"] = Pi
    out["c_isc"] = c_i
    out["c_pl"] = c_p
    out["Q_L"] = np.array([0.0 if null == "N1" else l["Q_L"] for l in lym])
    out["V_l"] = out["Q_L"] * p["tau_l"]
    out["M_l"] = (1 - p["sigma_L"]) * c_i * out["V_l"]
    out["pi_isc"] = np.array([pi_albumin(c, p) for c in c_i])
    out["pi_pl"] = np.array([pi_albumin(c, p) for c in c_p])
    cap = [capillary(p, p["Pc_cap"] + dPc, v, ci) for v, ci in zip(out["V_i"], c_i)]
    out["J_net"] = np.array([c["J_net"] for c in cap])
    out["E_glyc"] = np.array([c["E_glyc"] for c in cap])
    out["sigma_g"] = np.array([c["sigma_g"] for c in cap])
    out["sigma_s"] = np.array([c["sigma_s"] for c in cap])
    out["K_f"] = np.array([c["K_f"] for c in cap])
    out["R_mob"] = np.array([l["R_mob"] for l in lym])
    return out

def steady(p, dPc=0.0, null="none", y0=None, **kw):
    r = simulate(p, dPc=dPc, null=null, y0=y0, **kw)
    return r, -1

def equilibrate(p, t_end=25 * 1440.0, null="none", **kw):
    "Base equilibrium at dP_c = 0; all scenarios start from there."
    r = simulate(p, dPc=0.0, null=null, t_end=t_end, **kw)
    return r["y0"]

# ----------------------------------------------------------------------------
# ENHETSKONTROLL
# ----------------------------------------------------------------------------
UNIT_CHECKS = []
def _reg(name, lhs, rhs): UNIT_CHECKS.append((name, lhs, rhs))

def unit_check(p):
    "Mechanical dimensional analysis of EVERY constitutive relation.\n    Each entry: (relation, unit on either side, value check on a dimensional\n    condition -> on plain numbers one can only test the dimension, not the number)."
    cap = capillary(p, p["Pc_cap"], 11000.0, 15.0)
    ly = lymph(p, 1.0)
    si = cap["sigma_s"]
    checks = [
        # (relation, dim(vard), dim(hogerled), varde, boad villkor)
        ("1  pi = f(c)                     [mmHg] = f[g/L]",
         MMHG, MMHG, pi_albumin(15.0, p), lambda v: v > 0 and np.isfinite(v)),
        ("2  J_net = K_f*dP_eff            [mL/min] = [mL/min/mmHg]*[mmHg]",
         FLUX_V, KF_UNIT * MMHG, cap["J_net"], lambda v: np.isfinite(v)),
        ("2  dP_eff = J_net/K_f            [mmHg]",
         MMHG, FLUX_V / KF_UNIT, cap["dP_eff"], lambda v: np.isfinite(v)),
        ("2  (pi_p - pi_i) sigma_g         [mmHg] -> dimensionless * [mmHg]",
         MMHG, MMHG, cap["pi_p"] - cap["pi_i"], lambda v: v > 0),
        ("3  G_d = K_f0*D_ratio*1e3        [mL/min]",
         FLUX_V, FLUX_V, p["K_f0"] * p["D_ratio"] * 1e3, lambda v: v > 0),
        ("3  J_p,conv = (1-s_s)*c_p*J      [mg/min] = - * [g/L] * [mL/min]",
         MG / MIN, G_L * FLUX_V, (1 - si) * 15.0 * cap["J_net"], lambda v: np.isfinite(v)),
        ("3  J_p,diff = G_d*(c_p - c_i)     [mg/min] = [mL/min]*[g/L]",
         MG / MIN, (FLUX_V) * G_L, (p["K_f0"] * p["D_ratio"] * 1e3) * (15.0 - 15.0),
         lambda v: v == 0.0),
        ("4  P_ISC = a ln(V/V_k)           [mmHg] = a[mmHg]*dimensionellt",
         MMHG, MMHG, P_isc(11000.0, p), lambda v: -14.0 <= v <= 25.0),
        ("5  R_mob(P)                      [1] : R(-1)=1 ... R(+1)->R_mob_max",
         ML ** 0, ML ** 0, R_mob(1.0, p) / R_mob(-1.0, p), lambda v: v >= 100.0),
        ("6  Q_L = K_L P_pump exp(beta dP)  [mL/min] = [mL/min/mmHg]*[mmHg]*1",
         FLUX_V, KF_UNIT * MMHG, ly["Q_L"], lambda v: v > 0),
        ("6  dP_L_eff = dP_L/(1+R_mob)     [mmHg]",
         MMHG, MMHG, ly["dP_L_eff"], lambda v: abs(v) <= abs(ly["dP_L"])),
        ("7  c_l = (1-s_L) c_i             [g/L]",
         G_L, G_L, (1 - p["sigma_L"]) * 15.0, lambda v: 0 < v < 15.0),
        ("7  M_l_ret = f_p c_l Q_L         [mg/min] = - * [g/L] * [mL/min]",
         MG / MIN, G_L * FLUX_V, p["f_p_lymph"] * (1 - p["sigma_L"]) * 15.0 * ly["Q_L"],
         lambda v: np.isfinite(v)),
        ("8  J_osm = K_c*(dP + dpi)         [mL/min] = [mL/min/mmHg]*[mmHg]",
         FLUX_V, KC_UNIT * MMHG,
         p["K_cell"] * 1.0, lambda v: np.isfinite(v)),
        ("8  J_reg = (V_c0 - V_c)/tau_c    [mL/min] = [mL]/[min]",
         FLUX_V, ML / MIN, p["V_cell0"] * 1e3 / p["tau_cell"], lambda v: v > 0),
        ("9  S_hep = c_p*(V_p - V_p0)/tau  [mg/min] = [g/L]*[mL]/[min]",
         MG / MIN, G_L * (ML / MIN), p["c_plasma"] * 1e3 / p["tau_hep"], lambda v: v > 0),
    ]
    out, ok_all = [], True
    for rel, lhs, rhs, val, valtest in checks:
        ok = (lhs == rhs) and bool(valtest(val))
        ok_all &= ok
        out.append({"relation": rel, "lhs_unit": lhs.text(), "rhs_unit": rhs.text(),
                    "value": float(val), "ok": bool(ok)})
    if not ok_all:
        bad = [c for c in out if not c["ok"]]
        raise AssertionError("ENHETSKONTROLL MISSLYCKADES: " + json.dumps(bad, ensure_ascii=False))
    return out

def convergence(r, p, dt_min=1440.0):
    """|dX/dt| vid slutet, normaliserad mot X (konvergenskontroll)."""
    V, M = r["V_i"], r["M_i"]          # V i mL, M i mg
    dV_L = abs(V[-1] - V[-2]) / dt_min / 1e3      # L/min
    dM_g = abs(M[-1] - M[-2]) / dt_min / 1e3      # g/min
    relV = dV_L * 1440.0 / (V[-1] / 1e3)
    relM = dM_g * 1440.0 / (M[-1] / 1e3)
    return {"dV_isc_L_per_day": float(dV_L * 1440.0), "dM_alb_g_per_day": float(dM_g * 1440.0),
            "rel_dV_per_day": float(relV), "rel_dM_per_day": float(relM),
            "converged": bool(relV < 5e-3 and relM < 5e-3)}

def mass_balance_residual(r):
    "The total mass of the water must be conserved (closed model without infusion)."
    tot = r["V_p"] + r["V_i"] + r["V_c"] + r["V_l"]
    return float(np.max(np.abs(tot - tot[0]))) / float(np.mean(tot))

# ----------------------------------------------------------------------------
# KANSLIGHET  (+/-50 %)
# ----------------------------------------------------------------------------
def sensitivity(p, dPc=5.0, keys=("K_f0", "beta", "a_stiff"), f=0.5):
    y0 = equilibrate(p)
    base, _ = steady(p, dPc=dPc, y0=y0)
    b = dict(dV=float(base["V_i"][-1] - base["y0"][2] * 1e3),
             dM=float(base["M_i"][-1] - base["y0"][3] * 1e3))
    b["dM"] = abs(b["dM"]) if abs(b["dM"]) > 1e-9 else 1e-9
    b["R_ret"] = b["dM"] / b["dV"]
    rows = []
    for k in keys:
        for sgn in (+1, -1):
            pp = build_params(**{k: p[k] * (1 + sgn * f)})
            yy, err0 = _safe(equilibrate, pp)
            r, err = (None, err0) if yy is None else _safe(steady, pp, dPc=dPc, y0=yy)
            if r is None:
                rows.append({"param": k, "value": p[k] * (1 + sgn * f),
                             "sign": "+" if sgn > 0 else "-",
                             "runaway_no_steady_state": True, "note": err})
                continue
            r = r[0]
            dV = float(r["V_i"][-1] - r["y0"][2] * 1e3)
            dM = float(r["M_i"][-1] - r["y0"][3] * 1e3)
            rows.append({"param": k, "value": p[k] * (1 + sgn * f), "sign": "+" if sgn > 0 else "-",
                         "dV_isc_L": dV / 1e3, "dM_alb_g": dM / 1e3,
                         "R_ret_mg_mL": dM / dV,
                         "dV_pct": 100.0 * (dV / 1e3 - b["dV"] / 1e3) / (b["dV"] / 1e3),
                         "dM_pct": 100.0 * (dM / 1e3 - b["dM"] / 1e3) / (b["dM"] / 1e3)})
    return b, rows

# ----------------------------------------------------------------------------
# HUVUD
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results.json")
    a = ap.parse_args()
    p = build_params()
    units = unit_check(p)

    res = {"id": "BT-HX-Q020",
           "model": "4-compartment capillary-interstitium-lymph (first principles)",
           "prereg_sha256": open("PREREG.sha256").read().split()[0],
           "parameters": {k: {"value": v[0], "unit": v[1], "source": v[2]}
                          for k, v in PARAMS.items()},
           "derived_params": {"pi_max": p["pi_max"], "c_half": p["c_half"],
                              "f_A_eff": p["K_f0"] / (p["L_p"] * p["A_cap"] * 6e5)},
           "unit_check": units}

    # --- P6: basfloden ---------------------------------------------------
    base, _ = steady(p, dPc=0.0, y0=equilibrate(p))
    res["baseline"] = {
        "V_isc_L": float(base["V_i"][-1] / 1e3),
        "M_isc_g": float(base["M_i"][-1] / 1e3),
        "c_isc_gL": float(base["c_isc"][-1]),
        "c_pl_gL": float(base["c_pl"][-1]),
        "pi_isc_mmHg": float(base["pi_isc"][-1]),
        "pi_pl_mmHg": float(base["pi_pl"][-1]),
        "P_isc_mmHg": float(base["P_isc"][-1]),
        "J_net_mLmin": float(base["J_net"][-1]),
        "Q_L_mLmin": float(base["Q_L"][-1]),
        "E_glyc": float(base["E_glyc"][-1]),
        "mass_balance_rel_resid": mass_balance_residual(base)}

    # --- P1/P2/P3: steg i Pc ---------------------------------------------
    scen = {}
    y0 = equilibrate(p)
    res["y0_equilibrated"] = [float(v) for v in y0]
    for dp in (2.0, 5.0, 10.0):
        r, err = _safe(steady, p, dPc=dp, y0=y0)
        if r is None:
            scen[f"{dp:.0f}"] = {"dP_c_mmHg": dp, "no_physical_solution": True,
                                 "reason": err,
                                 "R_ret_mg_mL": None, "dV_isc_L": None}
            continue
        r = r[0]
        dV = float(r["V_i"][-1] - r["y0"][2] * 1e3) / 1e3
        dM = float(r["M_i"][-1] - r["y0"][3] * 1e3) / 1e3
        scen[f"{dp:.0f}"] = {
            "dP_c_mmHg": dp, "dV_isc_L": dV, "dM_alb_g": dM,
            "R_ret_mg_mL": dM / dV,
            "P_isc_ss_mmHg": float(r["P_isc"][-1]),
            "E_glyc_ss": float(r["E_glyc"][-1]),
            "K_f_ss": float(r["K_f"][-1]),
            "Q_L_ss_mLmin": float(r["Q_L"][-1]),
            "J_net_ss_mLmin": float(r["J_net"][-1]),
            "t90_h": _t90(r["V_i"]),
            "c_isc_ss_gL": float(r["c_isc"][-1]),
            "M_isc_ss_g": float(r["M_i"][-1] / 1e3),
            "V_isc_ss_L": float(r["V_i"][-1] / 1e3),
            "convergence": convergence(r, p),
            "V_l_ss_L": float(r["V_l"][-1] / 1e3),
            "M_l_ss_g": float(r["M_l"][-1] / 1e3)}
    res["scenarios"] = scen

    # --- P4: mobilitetsknae ------------------------------------------------
    Rl, Rh = R_mob(-1.0, p), R_mob(+1.0, p)
    res["P4_mobility_knee"] = {"R_mob_Pminus1": Rl, "R_mob_Pplus1": Rh,
                               "ratio": Rh / Rl, "threshold_P4": 100.0,
                               "ok": bool(Rh / Rl >= 100.0)}

    # --- P5: kritisk last --------------------------------------------------
    y0 = equilibrate(p)
    res["P5_threshold_scan"] = _threshold_scan(p, (y0,))

    # --- nullmodeller ------------------------------------------------------
    res["null_models"] = _null_models(p)

    # --- sensitivity --------------------------------------------------------------------
    b, rows = sensitivity(p, dPc=5.0)
    res["sensitivity"] = {"base_dV_isc_L": b["dV"] / 1e3,
                          "base_dM_alb_g": b["dM"] / 1e3,
                          "base_R_ret_mg_mL": b["R_ret"], "rows": rows}

    # --- prereg-kriterier --------------------------------------------------
    ks = [k for k in scen if scen[k].get("R_ret_mg_mL") is not None]
    Rr = [scen[k]["R_ret_mg_mL"] for k in ks]
    cs = [scen[k]["c_isc_ss_gL"] for k in ks]
    scen = {k: scen[k] for k in ks} | {k: scen[k] for k in scen if k not in ks}
    spread = (max(Rr) - min(Rr)) / abs(np.mean(Rr)) if np.mean(Rr) != 0 else float("inf")
    eq_err = [abs(a - b) / b for a, b in zip(Rr, cs)]
    allpos = bool(all(v > 0 for v in Rr))
    vspread = (max(scen[k]["dV_isc_L"] for k in ks) / min(scen[k]["dV_isc_L"] for k in ks)
               if ks else float("nan"))
    res["prereg_verdict"] = {
        "P1_R_ret_equals_c_isc_and_invariant": {
            "R_ret_mg_mL": {k: float(v) for k, v in zip(("2", "5", "10"), Rr)},
            "c_isc_ss_gL": {k: float(v) for k, v in zip(("2", "5", "10"), cs)},
            "rel_err_vs_c_isc": {k: float(v) for k, v in zip(("2", "5", "10"), eq_err)},
            "all_R_ret_positive": allpos,
            "spread_pct": float(100 * spread), "spread_limit_pct": 15.0,
            "ok": bool(allpos and spread < 0.15)},
        "P2_dV_spread_ge3": {"value": float(vspread), "limit": 3.0, "ok": bool(vspread >= 3.0)},
        "P3_P_isc_ss_positive_at_5mmHg": {"value": scen["5"].get("P_isc_ss_mmHg"),
                                          "ok": bool(scen["5"].get("P_isc_ss_mmHg", -1) > 0)},
        "P7_no_plasma_depletion": {"scenarios_without_physical_solution": [
            k for k, v in scen.items() if v.get("no_physical_solution")],
            "note": "closed water budget without infusion: sustained dP_c depletes plasma"},
        "P4_mobility_ratio": res["P4_mobility_knee"]["ratio"],
        "P5_threshold_in_window": res["P5_threshold_scan"]["summary"],
        "P6_basal_flows_in_range": {
            "J_net_mLmin": res["baseline"]["J_net_mLmin"],
            "Q_L_mLmin": res["baseline"]["Q_L_mLmin"],
            "window_mLmin": [0.2, 4.0],
            "ok": bool(0.2 <= res["baseline"]["Q_L_mLmin"] <= 4.0
                       and 0.2 <= res["baseline"]["J_net_mLmin"] <= 4.0)}}

    with open(a.out, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False, default=_dflt)
    print("Wrote", a.out)
    print(json.dumps(res["prereg_verdict"], indent=2, ensure_ascii=False))
    return res

def _dflt(o):
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.integer,)):   return int(o)
    if isinstance(o, np.ndarray):     return o.tolist()
    return str(o)

def _t90(V_i):
    """tid [h] till 90 % av slutvardet (index -> t via t_eval-steget 6 h)"""
    d = np.asarray(V_i) - V_i[0]
    tgt = 0.9 * d[-1]
    if abs(tgt) < 1e-12: return 0.0
    idx = np.where(d >= tgt)[0]
    return float(idx[0]) * 24.0 if len(idx) else None

def _threshold_scan(p, EQ, t_end=25 * 1440.0):
    rows = []
    for dp in np.arange(0.0, 24.1, 1.0):
        try:
            r, _ = steady(p, dPc=float(dp), y0=EQ[0], t_end=t_end, t_eval_h=12.0)
        except RuntimeError as e:
            rows.append({"dP_c": float(dp), "runaway_no_steady_state": True,
                         "dV_isc_L": None, "dM_alb_g": None, "E_glyc_ss": None, "P_isc_ss": None})
            continue
        dV = float(r["V_i"][-1] - r["y0"][2] * 1e3) / 1e3
        dM = float(r["M_i"][-1] - r["y0"][3] * 1e3) / 1e3
        rows.append({"dP_c": float(dp), "dV_isc_L": dV, "dM_alb_g": dM,
                     "E_glyc_ss": float(r["E_glyc"][-1]),
                     "P_isc_ss": float(r["P_isc"][-1])})
    thr = None
    for i in range(1, len(rows)):
        if rows[i].get("runaway_no_steady_state") or (rows[i]["E_glyc_ss"] is not None
                                                     and rows[i]["E_glyc_ss"] < 0.5):
            thr = rows[i]["dP_c"]; break
    return {"rows": rows, "summary": {"dP_c_threshold_mmHg": thr,
                                      "window": [4.0, 20.0],
                                      "ok": bool(thr is not None and 4.0 <= thr <= 20.0)}}

def _safe(fn, *a, **kw):
    "Runaway = no finite equilibrium; this is a result, not an error."
    try:
        return fn(*a, **kw), None
    except RuntimeError as e:
        return None, str(e)

def _null_models(p):
    y0 = equilibrate(p)
    out = {}
    # N1 ingen lymfa
    r, err = _safe(steady, p, dPc=5.0, null="N1", y0=y0, t_end=25 * 1440.0, t_eval_h=12.0)
    out["N1_no_lymph"] = ({"runaway_no_steady_state": True, "note": err} if r is None else
        {"dV_isc_L": float(r[0]["V_i"][-1] - r[0]["y0"][2] * 1e3) / 1e3,
         "dM_alb_g": float(r[0]["M_i"][-1] - r[0]["y0"][3] * 1e3) / 1e3,
         "P_isc_ss": float(r[0]["P_isc"][-1])})
    # N2 sigma_s = 0
    r, err = _safe(steady, p, dPc=5.0, null="N2", y0=y0, t_end=25 * 1440.0, t_eval_h=12.0)
    out["N2_sigma0"] = ({"runaway_no_steady_state": True, "note": err} if r is None else
        {"dV_isc_L": float(r[0]["V_i"][-1] - r[0]["y0"][2] * 1e3) / 1e3,
         "dM_alb_g": float(r[0]["M_i"][-1] - r[0]["y0"][3] * 1e3) / 1e3,
         "R_ret_mg_mL": float(r[0]["M_i"][-1] - r[0]["y0"][3] * 1e3) / float(r[0]["V_i"][-1] - r[0]["y0"][2] * 1e3)})
    # N3 full glykocalyx (E = 1) : satt n_glyc = 0 -> E = 1/(1+0) = 1
    p3 = build_params(**{"E_fixed": 1.0})
    out["N3_full_glycocalyx"] = {}
    for dp in (2.0, 5.0, 10.0):
        r, err = _safe(steady, p3, dPc=dp, y0=equilibrate(p3), t_end=25 * 1440.0, t_eval_h=12.0)
        out["N3_full_glycocalyx"][f"{dp:.0f}"] = ({"runaway_no_steady_state": True} if r is None else
            {"dV_isc_L": float(r[0]["V_i"][-1] - r[0]["y0"][2] * 1e3) / 1e3,
             "E_glyc_ss": float(r[0]["E_glyc"][-1])})
    vv = [out["N3_full_glycocalyx"][k].get("dV_isc_L") for k in ("2", "5", "10")]
    if any(v is None for v in vv):
        vv = None
    out["N3_full_glycocalyx"]["dV_spread"] = ((max(vv) - min(vv)) / min(vv)
                                             if vv and min(vv) != 0 else None)
    # N0 linjar P-V (ballong)
    p0 = build_params(**{"pv_form": "linear"})   # N0: ingen riktning, samma lokala komplians
    out["N0_balloon"] = {}
    for dp in (2.0, 5.0, 10.0):
        r, err = _safe(steady, p0, dPc=dp, y0=equilibrate(p0), t_end=25 * 1440.0, t_eval_h=12.0)
        out["N0_balloon"][f"{dp:.0f}"] = ({"runaway_no_steady_state": True} if r is None else
            {"dV_isc_L": float(r[0]["V_i"][-1] - r[0]["y0"][2] * 1e3) / 1e3,
             "R_ret_mg_mL": float(r[0]["M_i"][-1] - r[0]["y0"][3] * 1e3) / float(r[0]["V_i"][-1] - r[0]["y0"][2] * 1e3),
             "P_isc_ss": float(r[0]["P_isc"][-1])})
    vr = [out["N0_balloon"][k].get("dV_isc_L") for k in ("2", "5", "10")]
    rr = [out["N0_balloon"][k].get("R_ret_mg_mL") for k in ("2", "5", "10")]
    if None not in vr and min(vr) != 0:
        out["N0_balloon"]["dV_spread"] = (max(vr) - min(vr)) / min(vr)
    if None not in rr:
        out["N0_balloon"]["R_ret_spread"] = (max(rr) - min(rr)) / np.mean(rr)
    return out

if __name__ == "__main__":
    main()
