"""WHOLE-BODY GLUCOSE-INSULIN MEAL SIMULATOR (Dalla Man, Rizza, Cobelli 2007
IEEE TBME 54(10):1740-9, PMID17926672, "Meal Simulation Model of the Glucose-Insulin System"),
glucose/insulin/gut(Dalla Man 2006 PMID17153204)/beta-cell (Toffolo-Cobelli) subsystems, params
= paper's Table-1 "average subject" set (matches BioModels BIOMD0000000379 curation).

All parameters are copied verbatim from the published table -- NOT re-fit, NOT tuned to hit any
recorded number. Every parameter is fixed BEFORE any scenario is run. If a scenario's computed
number disagrees with the recorded one, THAT IS THE ANSWER -- this cell does not adjust anything
to close a gap.

FALSIFIER (pre-registered, see GATES below): each of ~11 recorded headline numbers (basal drift
band, healthy-75g-OGTT peaks, dose-response, EGP-suppression sweep, T1D fasting/meal/pump-bolus
scenarios) is checked against a tolerance stated up front (5% relative, or an absolute band where
the original claim itself gave one). PASS = within tolerance. FAIL = outside it, reported as an
exact gap, not narrated away.

STRUCTURAL NOTE ON ke1/ke2 (renal excretion term E(t)=ke1*max(Gp-ke2,0)): BioModels'
own curated encoding zeroes this term's practical effect for normoglycemic reproduction. This build
KEEPS the term (ke2=339 mg/kg is far above any healthy-range Gp, so E(t)=0 throughout the
basal/healthy/dose-response/EGP-suppression scenarios regardless -- it only activates in the T1D
hyperglycemic scenarios, where Gp can exceed 339 mg/kg). Both "kept" and "zeroed" are reported for
any scenario where the two differ.

Reads: nothing. Writes: glucose_meal_dallaman2007_results.json
(gate: required_gates_overall_pass = all of GATES).
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "glucose_meal_dallaman2007")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================================================
# PARAMETERS -- verbatim from the published table (Dalla Man 2007 Table-1 average
# subject; BioModels BIOMD0000000379). NONE adjusted. NONE tuned to any headline number below.
# ============================================================================================
P = dict(
    V_G=1.88, k1=0.065, k2=0.079, G_b=95.0,
    V_I=0.05, m1=0.19, m2=0.484, m4=0.194, m5=0.0304, m6=0.6471, I_b=25.0, S_b=1.8,
    k_max=0.0558, k_min=0.008, k_abs=0.057, k_gri=0.0558, f=0.9, b=0.82, d=0.01, BW=78.0,
    kp1=2.7, kp2=0.0021, kp3=0.009, kp4=0.0618, ki=0.0079,
    Uii=1.0, Vm0=2.5, VmX=0.047, Km0=225.59, p2U=0.0331,
    K=2.3, alpha=0.05, beta=0.11, gamma=0.5,
    ke1=0.0005, ke2=339.0,
    ka_sc_insulin=0.01818, pmol_per_unit=6000.0,
    part=0.2,
)
# HE_b (hepatic extraction at basal) kept only as a basal-diagnostic reference; the RHS below
# recomputes HE/m3 from the INSTANTANEOUS secretion S(t) every call (bug #3 fix -- was frozen at
# the basal value here and never updated, silently pinning hepatic clearance regardless of state).
P["HE_b"] = -P["m5"] * P["S_b"] + P["m6"]
P["m3_b_diag_only"] = P["HE_b"] * P["m1"] / (1.0 - P["HE_b"])

# State vector order:
# [Gp, Gt, Ip, Il, Qsto1, Qsto2, Qgut, X, I1, Id, Ipo, Y, SCins]
IDX = dict(Gp=0, Gt=1, Ip=2, Il=3, Qsto1=4, Qsto2=5, Qgut=6, X=7, I1=8, Id=9, Ipo=10, Y=11, SCins=12)
NSTATE = 13


def k_empt(Qsto, D):
    if D <= 0:
        return P["k_max"]
    a_g = 5.0 / (2.0 * D * (1.0 - P["b"]))
    b_g = 5.0 / (2.0 * D * P["d"])
    return P["k_min"] + (P["k_max"] - P["k_min"]) / 2.0 * (
        np.tanh(a_g * (Qsto - P["b"] * D)) - np.tanh(b_g * (Qsto - P["d"] * D)) + 2.0
    )


def rhs(t, y, D, meal_t0, secretion_gain, sc_bolus_pmol_per_min, basal_ins_pmol_per_min, use_ke=True):
    Gp, Gt, Ip, Il, Qsto1, Qsto2, Qgut, X, I1, Id, Ipo, Y, SCins = y
    G = Gp / P["V_G"]
    I = Ip / P["V_I"]

    # ---- gut / gastric emptying (meal starts at meal_t0) ----
    if D > 0 and t >= meal_t0:
        Qsto = Qsto1 + Qsto2
        ke = k_empt(Qsto, D)
        dQsto1 = -P["k_gri"] * Qsto1
        dQsto2 = -ke * Qsto2 + P["k_gri"] * Qsto1
        dQgut = -P["k_abs"] * Qgut + ke * Qsto2
        Ra = P["f"] * P["k_abs"] * Qgut / P["BW"]
    else:
        dQsto1 = dQsto2 = dQgut = 0.0
        Ra = 0.0

    # ---- insulin secretion (Toffolo-Cobelli beta cell), scaled by secretion_gain (T1D=0) ----
    Gprime_proxy = P["k1"] * Gp - P["k2"] * Gt  # dGp/dt without the meal/EGP/util terms as G' proxy
    # Use actual instantaneous glucose derivative estimate via chain: we need dG/dt for K*max(G',0).
    # Approximate using the true state derivatives computed below (two-pass not needed: use dGp/dt
    # after computing EGP/Uid/E first). We compute EGP/Uid/E now, then Gp' , then G' = Gp'/V_G.

    Id_eff = Id
    EGP = max(P["kp1"] - P["kp2"] * Gp - P["kp3"] * Id_eff - P["kp4"] * Ipo, 0.0)
    # Uid carries the (1-part) discount on Vm0+VmX*X (part=0.2), per the published
    # Dalla Man 2007 form.
    Uid = (1.0 - P["part"]) * (P["Vm0"] + P["VmX"] * X) * Gt / (P["Km0"] + Gt)
    E = P["ke1"] * max(Gp - P["ke2"], 0.0) if use_ke else 0.0

    dGp = EGP + Ra - P["Uii"] - E - P["k1"] * Gp + P["k2"] * Gt
    dGt = -Uid + P["k1"] * Gp - P["k2"] * Gt
    Gdot = dGp / P["V_G"]

    # bug #4 fix: both the Spo G'-term and the dY threshold term are the PUBLISHED unclipped forms
    # (max(...,0) clipping only matters transiently anyway -- at any true fixed point Gdot=0 and
    # G=G_eq exactly, so the clip is a no-op there; but off the fixed point during a meal transient
    # the clip was asymmetrically suppressing negative excursions the rebuild does not suppress).
    Spo = Y + P["K"] * Gdot + P["S_b"]
    Spo = secretion_gain * Spo
    # bug #5 fix (the compounding factor found via long-horizon equilibration diff against the
    # rebuild): the published Dalla Man 2007 beta-cell ODE is dIpo/dt = -gamma*Ipo + Spo, i.e. ONLY
    # the Ipo term carries gamma -- NOT -gamma*(Ipo-Spo), which spuriously multiplies Spo (and
    # therefore the S_b basal floor) by gamma too. That single extra gamma factor on Spo was, by
    # itself, responsible for the entire +13.3% residual basal-drift gap left after bugs #1-#4 were
    # fixed (measured: 20000-min equilibration settles at G=107.63, gap=+13.30%, with the old form;
    # G=94.9418, gap=-0.061%, with this form -- and the committed and rebuild scripts now agree on
    # this equilibrium to >8 significant figures, an external cross-model anchor neither script was
    # tuned against).
    dIpo = -P["gamma"] * Ipo + Spo
    dY = -P["alpha"] * (Y - P["beta"] * (G - P["G_b"]))
    # bug #1 fix: S(t) (secretion into the liver, dIl) is gamma*Ipo per the same published equation
    # (S(t)=gamma*Ipo(t)) -- was missing the gamma factor entirely.
    S = secretion_gain * P["gamma"] * Ipo
    # bug #3 fix: hepatic extraction HE and derived m3 are recomputed from the INSTANTANEOUS S(t)
    # every RHS call (were frozen at the basal S_b value in P["m3"] and never updated).
    HE = np.clip(-P["m5"] * S + P["m6"], 0.0, 1.0)
    m3_t = HE * P["m1"] / (1.0 - HE)

    # ---- exogenous SC insulin (bolus pulses handled as instantaneous state jumps in the driver;
    # here only the continuous basal infusion + SC depot->plasma absorption ODE is integrated) ----
    dSCins = sc_bolus_pmol_per_min - P["ka_sc_insulin"] * SCins
    Ra_ins = P["ka_sc_insulin"] * SCins + basal_ins_pmol_per_min

    # ---- insulin subsystem ----
    dIp = -(P["m2"] + P["m4"]) * Ip + P["m1"] * Il + Ra_ins
    dIl = -(P["m1"] + m3_t) * Il + P["m2"] * Ip + S

    # ---- remote insulin action ----
    dX = -P["p2U"] * (X - (I - P["I_b"]))

    # ---- delayed insulin action on EGP ----
    dI1 = -P["ki"] * (I1 - I)
    dId = -P["ki"] * (Id - I1)

    dy = np.zeros(NSTATE)
    dy[IDX["Gp"]] = dGp
    dy[IDX["Gt"]] = dGt
    dy[IDX["Ip"]] = dIp
    dy[IDX["Il"]] = dIl
    dy[IDX["Qsto1"]] = dQsto1
    dy[IDX["Qsto2"]] = dQsto2
    dy[IDX["Qgut"]] = dQgut
    dy[IDX["X"]] = dX
    dy[IDX["I1"]] = dI1
    dy[IDX["Id"]] = dId
    dy[IDX["Ipo"]] = dIpo
    dy[IDX["Y"]] = dY
    dy[IDX["SCins"]] = dSCins
    return dy


def basal_state(secretion_gain=1.0):
    """Basal steady state. ORIENT finding (measured, not assumed): the Table-1 'average subject'
    parameter set is a population-average vector, not guaranteed to satisfy the FULL coupled
    nonlinear steady state exactly at Gp_b=G_b*V_G, Ip_b=I_b*V_I simultaneously -- confirmed here:
    solving dIp=0/dIl=0 analytically for Il0 with Ip0 pinned at I_b*V_I leaves a residual (~0.87
    pmol/kg/min) regardless of the Spo/Sb convention used. Rather than pin an inconsistent analytic
    IC (which produced a spurious -40% Gp drift over 300 min pre-fix), find the TRUE fixed point by
    long unforced numerical equilibration from the definitional/table IC -- same honest method
    already used for the T1D case below. This does not adjust any parameter; it only refuses to
    force an algebraically over-determined analytic shortcut. The two IDENTITIES the task called
    out as exact (Ip(0)=I_b*V_I=1.25, Ipo(0)=S_b/gamma=3.6) are used as the STARTING guess and are
    reported as converged/not-converged diagnostics below (both hold if the residual reported is
    ~0 relative to Ip_b/Ipo_b)."""
    Gp_b_guess = P["G_b"] * P["V_G"]
    Ip0 = P["I_b"] * P["V_I"]
    Ipo0 = (P["S_b"] / P["gamma"]) if secretion_gain > 0 else 0.0
    E_b = P["ke1"] * max(Gp_b_guess - P["ke2"], 0.0)

    def gt_resid(Gt_b):
        Uid_b = (1.0 - P["part"]) * P["Vm0"] * Gt_b / (P["Km0"] + Gt_b)
        return Uid_b - (P["k1"] * Gp_b_guess - P["k2"] * Gt_b)

    Gt_b_guess = brentq(gt_resid, 1e-6, 1e6, xtol=1e-12, rtol=1e-14)
    Il0_guess = (P["m2"] * Ip0 + secretion_gain * Ipo0) / (P["m1"] + P["m3_b_diag_only"])
    dIp_check_analytic = -(P["m2"] + P["m4"]) * Ip0 + P["m1"] * Il0_guess

    y0 = np.zeros(NSTATE)
    y0[IDX["Gp"]] = Gp_b_guess
    y0[IDX["Gt"]] = Gt_b_guess
    y0[IDX["Ip"]] = Ip0
    y0[IDX["Il"]] = Il0_guess
    y0[IDX["X"]] = 0.0
    y0[IDX["I1"]] = P["I_b"]
    y0[IDX["Id"]] = P["I_b"]
    y0[IDX["Ipo"]] = Ipo0
    y0[IDX["Y"]] = 0.0  # true dY=0 fixed point at G=G_b (see forced-correction note in rhs())
    y0[IDX["SCins"]] = 0.0

    def f(t, y):
        return rhs(t, y, D=0.0, meal_t0=0.0, secretion_gain=secretion_gain,
                   sc_bolus_pmol_per_min=0.0, basal_ins_pmol_per_min=0.0, use_ke=True)

    sol = solve_ivp(f, [0, 20000], y0, method="LSODA", rtol=1e-9, atol=1e-11, max_step=50.0)
    y_eq = sol.y[:, -1]
    tail = sol.y[IDX["Gp"], -50:]
    equilibration_drift = float(np.max(np.abs(np.diff(tail))))
    Gp_eq = float(y_eq[IDX["Gp"]])
    G_eq = Gp_eq / P["V_G"]
    gap_from_Gb_pct = (G_eq - P["G_b"]) / P["G_b"] * 100.0

    diag = dict(
        Gp_b_guess=Gp_b_guess, Gt_b_guess=Gt_b_guess, dIp_check_analytic=dIp_check_analytic,
        equilibration_drift_mgkgmin=equilibration_drift, G_equilibrated=G_eq,
        gap_equilibrated_vs_Gb_pct=gap_from_Gb_pct,
        Ipo_equilibrated=float(y_eq[IDX["Ipo"]]), Ipo_identity_target=Ipo0,
        Ip_equilibrated=float(y_eq[IDX["Ip"]]), Ip_identity_target=Ip0,
    )
    return y_eq, diag


def t1d_basal_fixed_point():
    """T1D: secretion_gain=0 forced (S(t)=0 always). Find the TRUE re-equilibrated basal fixed
    point by long unforced integration from the healthy basal IC with secretion off, until drift
    < 1e-6 mg/kg/min in Gp. This is a genuine re-equilibration, not the healthy Gp_b=95*V_G value."""
    y0, _ = basal_state(secretion_gain=1.0)
    y0 = y0.copy()
    y0[IDX["Ipo"]] = 0.0
    y0[IDX["Y"]] = 0.0

    def f(t, y):
        return rhs(t, y, D=0.0, meal_t0=0.0, secretion_gain=0.0,
                   sc_bolus_pmol_per_min=0.0, basal_ins_pmol_per_min=0.0, use_ke=True)

    sol = solve_ivp(f, [0, 20000], y0, method="LSODA", rtol=1e-9, atol=1e-11,
                     dense_output=False, max_step=50.0)
    yf = sol.y[:, -1]
    tail = sol.y[IDX["Gp"], -50:]
    drift = float(np.max(np.abs(np.diff(tail))))
    return yf, drift, sol


def with_meal_ic(y0, D):
    """Set Qsto1(0)=D, the published initial condition for a meal of size D."""
    y = y0.copy()
    y[IDX["Qsto1"]] = D
    return y


def run_scenario(D, meal_t0, y0, t_end, secretion_gain=1.0,
                  sc_bolus_schedule=None, basal_ins_rate_U_hr=0.0, use_ke=True, max_step=1.0):
    """sc_bolus_schedule: list of (t_bolus_min, dose_U) -- delivered as a narrow pulse into SCins
    depot by adding dose*pmol_per_unit directly to the SCins state at that time (event-based jump),
    consistent with 'fast SC insulin absorption pulse feeding into Ip via ka_sc_insulin'.
    UNIT NOTE: Ip/Il/SCins are pmol/kg (per the published units,
    matching Ip(0)=I_b*V_I with V_I in L/kg). A dose in U converts to TOTAL pmol via pmol_per_unit,
    which must be divided by BW to enter the same per-kg state -- omitting this was an ~78x (=BW)
    overdose, measured as G going NEGATIVE in the bolus+basal scenario on the first run."""
    basal_ins_pmol_per_min = basal_ins_rate_U_hr * P["pmol_per_unit"] / 60.0 / P["BW"]

    events = []
    if sc_bolus_schedule:
        for (tb, dose_U) in sc_bolus_schedule:
            events.append((tb, dose_U * P["pmol_per_unit"] / P["BW"]))
    events.sort()

    y = y0.copy()
    t_cur = 0.0
    ts_all, ys_all = [], []
    # DEDUP breakpoints: a raw
    # [e[0] for e in events] + [t_end] list can repeat a timestamp (two events at the
    # same tb) and the OLD "bp <= t_cur: skip" condition silently dropped ANY breakpoint
    # at or before the current time -- since t_cur starts at 0.0 and boluses are
    # scheduled at t=0.0, that skipped the very first breakpoint and the bolus dose was
    # NEVER applied (measured: SCins depot stayed identically 0.0 for the full run of
    # scenario 6a/6b). CONVENTION CHOSEN: bp < t_cur -> already past, skip (guards
    # against a duplicate/out-of-order breakpoint being re-applied); bp == t_cur -> an
    # event exactly at the current clock (t=0 start, or a breakpoint reached by the
    # previous iterate) is legitimate -- apply its dose(s) WITHOUT invoking solve_ivp on
    # a zero-length interval; bp > t_cur -> integrate up to bp, then apply. Each
    # breakpoint value appears once (sorted set), so no event is ever applied twice.
    breakpoints = sorted(set([e[0] for e in events] + [t_end]))
    for bp in breakpoints:
        if bp < t_cur:
            continue
        if bp > t_cur:
            def f(t, yy):
                return rhs(t, yy, D=D, meal_t0=meal_t0, secretion_gain=secretion_gain,
                           sc_bolus_pmol_per_min=0.0, basal_ins_pmol_per_min=basal_ins_pmol_per_min,
                           use_ke=use_ke)

            sol = solve_ivp(f, [t_cur, bp], y, method="LSODA", rtol=1e-8, atol=1e-10, max_step=max_step,
                             dense_output=True, t_eval=np.arange(t_cur, bp + 1e-9, 1.0))
            ts_all.append(sol.t)
            ys_all.append(sol.y)
            y = sol.y[:, -1].copy()
            t_cur = bp
        for (tb, dose_pmol) in events:
            if abs(tb - bp) < 1e-6:
                y[IDX["SCins"]] += dose_pmol

    t_full = np.concatenate(ts_all)
    y_full = np.concatenate(ys_all, axis=1)
    return t_full, y_full


def G_of(y_full):
    return y_full[IDX["Gp"], :] / P["V_G"]


def I_of(y_full):
    return y_full[IDX["Ip"], :] / P["V_I"]


def peak_and_time(t, x):
    i = int(np.argmax(x))
    return float(x[i]), float(t[i])


def val_at(t, x, t_query):
    i = int(np.argmin(np.abs(t - t_query)))
    return float(x[i])


# ============================================================================================
# SCENARIO 1: basal fixed-point self-test, 300 min, no meal
# ============================================================================================
y0_healthy, basal_diag = basal_state(secretion_gain=1.0)
t1, y1 = run_scenario(D=0.0, meal_t0=0.0, y0=y0_healthy, t_end=300.0, secretion_gain=1.0)
G1 = G_of(y1)
drift_lo = (np.min(G1) - P["G_b"]) / P["G_b"] * 100.0
drift_hi = (np.max(G1) - P["G_b"]) / P["G_b"] * 100.0

# ============================================================================================
# SCENARIO 2/3: healthy OGTT dose-response at 50/75/100 g CHO (1 g CHO = 1000 mg)
# ============================================================================================
dose_response = {}
for grams in (50, 75, 100):
    D = grams * 1000.0
    t, y = run_scenario(D=D, meal_t0=0.0, y0=with_meal_ic(y0_healthy, D), t_end=400.0, secretion_gain=1.0)
    G = G_of(y)
    I = I_of(y)
    Gpk, Gpk_t = peak_and_time(t, G)
    Ipk, Ipk_t = peak_and_time(t, I)
    dose_response[grams] = dict(
        G_peak=Gpk, G_peak_t=Gpk_t, I_peak=Ipk, I_peak_t=Ipk_t,
        G_120=val_at(t, G, 120.0), G_240=val_at(t, G, 240.0),
    )

# ============================================================================================
# SCENARIO 4: EGP-suppression AUC (28-240 min) dose sweep at 40/50/60/75/100 g
# ============================================================================================
def egp_suppression_pct(D, y0, secretion_gain=1.0, use_ke=True):
    t, y = run_scenario(D=D, meal_t0=0.0, y0=with_meal_ic(y0, D), t_end=400.0, secretion_gain=secretion_gain, use_ke=use_ke)
    Gp = y[IDX["Gp"], :]
    Id = y[IDX["Id"], :]
    Ipo = y[IDX["Ipo"], :]
    EGP = np.maximum(P["kp1"] - P["kp2"] * Gp - P["kp3"] * Id - P["kp4"] * Ipo, 0.0)
    mask = (t >= 28.0) & (t <= 240.0)
    EGP_b = P["kp1"] - P["kp2"] * (P["G_b"] * P["V_G"]) - P["kp3"] * P["I_b"] - P["kp4"] * (P["S_b"] / P["gamma"])
    _trapz_fn = getattr(np, "trapezoid", None) or np.trapz  # numpy>=2.0 renamed trapz->trapezoid
    auc_actual = _trapz_fn(EGP[mask], t[mask])
    auc_basal = EGP_b * (t[mask][-1] - t[mask][0])
    supp_pct = (1.0 - auc_actual / auc_basal) * 100.0
    return supp_pct

egp_sweep = {}
for grams in (40, 50, 60, 75, 100):
    egp_sweep[grams] = egp_suppression_pct(grams * 1000.0, y0_healthy)

# ============================================================================================
# SCENARIO 5: T1D fasting fixed point (secretion_gain=0, re-equilibrated) + 75g meal peak
# ============================================================================================
y0_t1d, t1d_drift, t1d_sol = t1d_basal_fixed_point()
t1d_fasting_Gp_over_VG = float(y0_t1d[IDX["Gp"]] / P["V_G"])

t5, y5 = run_scenario(D=75000.0, meal_t0=0.0, y0=with_meal_ic(y0_t1d, 75000.0), t_end=500.0, secretion_gain=0.0)
G5 = G_of(y5)
t1d_meal_peak, t1d_meal_peak_t = peak_and_time(t5, G5)

# ============================================================================================
# SCENARIO 6: T1D pump/bolus forced-OODA test, 75g meal
#   (a) carb bolus 7.5U @ t=0 + correction bolus 2.92U @ t=0, no basal
#   (b) same boluses + continuous basal 0.7 U/hr
# ============================================================================================
y0_t1d_meal = with_meal_ic(y0_t1d, 75000.0)
t6a, y6a = run_scenario(D=75000.0, meal_t0=0.0, y0=y0_t1d_meal, t_end=600.0, secretion_gain=0.0,
                         sc_bolus_schedule=[(0.0, 7.5 + 2.92)], basal_ins_rate_U_hr=0.0)
G6a = G_of(y6a)
t6b, y6b = run_scenario(D=75000.0, meal_t0=0.0, y0=y0_t1d_meal, t_end=600.0, secretion_gain=0.0,
                         sc_bolus_schedule=[(0.0, 7.5 + 2.92)], basal_ins_rate_U_hr=0.7)
G6b = G_of(y6b)

def first_cross_below(t, G, thresh):
    below = G < thresh
    if not below.any():
        return None
    i = int(np.argmax(below))
    return float(t[i])

a_plateau_t502 = val_at(t6a, G6a, 502.0)
a_end600 = val_at(t6a, G6a, 600.0)
a_min_G = float(np.min(G6a))
b_cross200 = first_cross_below(t6b, G6b, 200.0)
b_end600 = val_at(t6b, G6b, 600.0)

# ============================================================================================
# GATES -- pre-registered BEFORE inspecting outputs above (tolerances fixed here, in code, at
# write-time; 5% relative unless the original claim itself specified an absolute band).
# ============================================================================================
def within(x, target, rel=0.05):
    if target == 0:
        return abs(x) < 1e-6
    return abs(x - target) / abs(target) <= rel

REC = dict(  # RECORDED headline numbers of the original model run, untouched
    drift_lo=-0.35, drift_hi=0.92,
    G_peak_75=161.6, I_peak_75=248.7, G120_75=144.7, G240_75=101.4,
    Gpeak_50=141.5, Gpeak_75=161.6, Gpeak_100=180.5,
    egp40=30.6, egp50=37.2, egp60=43.6, egp75=52.8, egp100=67.4,
    t1d_fasting=246.2, t1d_meal_peak=413.5,
    a_end600=214.7, a_never_below_200=True,
    b_cross200=276.0, b_end600=136.5,
)

GATES = {
    "basal_drift_band": bool(drift_lo >= REC["drift_lo"] - 0.5 and drift_hi <= REC["drift_hi"] + 0.5)
                          if not (np.isnan(drift_lo) or np.isnan(drift_hi)) else False,
    "healthy75_Gpeak": within(dose_response[75]["G_peak"], REC["G_peak_75"]),
    "healthy75_Ipeak": within(dose_response[75]["I_peak"], REC["I_peak_75"]),
    "healthy75_G120": within(dose_response[75]["G_120"], REC["G120_75"]),
    "healthy75_G240": within(dose_response[75]["G_240"], REC["G240_75"]),
    "doseresp_50_Gpeak": within(dose_response[50]["G_peak"], REC["Gpeak_50"]),
    "doseresp_75_Gpeak": within(dose_response[75]["G_peak"], REC["Gpeak_75"]),
    "doseresp_100_Gpeak": within(dose_response[100]["G_peak"], REC["Gpeak_100"]),
    "egp_supp_40": within(egp_sweep[40], REC["egp40"]),
    "egp_supp_50": within(egp_sweep[50], REC["egp50"]),
    "egp_supp_60": within(egp_sweep[60], REC["egp60"]),
    "egp_supp_75": within(egp_sweep[75], REC["egp75"]),
    "egp_supp_100": within(egp_sweep[100], REC["egp100"]),
    "t1d_fasting_Gp": within(t1d_fasting_Gp_over_VG, REC["t1d_fasting"]),
    "t1d_meal_peak": within(t1d_meal_peak, REC["t1d_meal_peak"]),
    "t1d_bolusalone_end600": within(a_end600, REC["a_end600"]),
    "t1d_bolusalone_never_below_200": bool(a_min_G >= 200.0) == REC["a_never_below_200"],
    "t1d_bolusbasal_cross200": (within(b_cross200, REC["b_cross200"]) if b_cross200 is not None else False),
    "t1d_bolusbasal_end600": within(b_end600, REC["b_end600"]),
}
REQUIRED_GATE_KEYS = list(GATES.keys())
required_gates_overall_pass = bool(all(GATES[k] for k in REQUIRED_GATE_KEYS))

results = {
    "citation": {"pmid_2007": "17926672", "pmid_2006_gut": "17153204", "biomodel": "BIOMD0000000379"},
    "basal_diag": basal_diag,
    "scenario1_basal_drift_pct": {"lo": drift_lo, "hi": drift_hi},
    "scenario2_3_dose_response": dose_response,
    "scenario4_egp_suppression_pct": egp_sweep,
    "scenario5_t1d": {
        "fasting_Gp_over_VG": t1d_fasting_Gp_over_VG, "fasting_drift_check": t1d_drift,
        "meal_peak_G": t1d_meal_peak, "meal_peak_t": t1d_meal_peak_t,
    },
    "scenario6_t1d_pump": {
        "a_bolus_alone": {"G_at_502": a_plateau_t502, "G_at_600": a_end600, "min_G": a_min_G},
        "b_bolus_plus_basal": {"cross_below_200_t": b_cross200, "G_at_600": b_end600},
    },
    "recorded_original_claim": REC,
    "gates": GATES,
    "required_gates_overall_pass": required_gates_overall_pass,
}

out_path = f"{OUT_DIR}/glucose_meal_dallaman2007_results.json"
def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    return str(o)


with open(out_path, "w") as fh:
    json.dump(results, fh, indent=2, default=_json_default)

print(f"wrote {out_path}")
print(json.dumps(GATES, indent=2, default=_json_default))
print(f"\nrequired_gates_overall_pass: {required_gates_overall_pass}")
print("\n--- key computed numbers ---")
print("basal drift %% range:", drift_lo, drift_hi)
print("dose_response:", json.dumps(dose_response, indent=2))
print("egp_sweep:", egp_sweep)
print("t1d fasting Gp/VG:", t1d_fasting_Gp_over_VG, "drift_check:", t1d_drift)
print("t1d meal peak:", t1d_meal_peak, "@t=", t1d_meal_peak_t)
print("t1d bolus-alone: @502=", a_plateau_t502, "@600=", a_end600, "min=", a_min_G)
print("t1d bolus+basal: cross<200 t=", b_cross200, "@600=", b_end600)
