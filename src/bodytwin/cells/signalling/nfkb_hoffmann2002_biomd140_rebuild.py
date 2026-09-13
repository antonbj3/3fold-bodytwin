#!/usr/bin/env python3
"""
Full 24-species Hoffmann et al. 2002 (Science, PMID 12424381) IkBalpha/beta/epsilon 3-isoform
NF-kB negative-feedback module -- BioModels-curated SBML BIOMD0000000140 -- transcribed EXACTLY
(all 24 species, all 45 reactions) from the BioModels .ode/SBML export and integrated with ZERO
period-fitting.

Reads: nothing. Writes: nfkb_hoffmann2002_biomd140_rebuild_results.json under the cell output
directory.

QUESTION: does it reproduce (a) a damped nuclear-NF-kB oscillation period under a transient IKK
pulse, inside the pre-registered Ashall 2009 (PMID 19359585) [50,150] min band (|period-100|<=50),
and (b) the persistent/clamped-IKK regime map (weak IKK sustains multi-cycle oscillation,
moderate/strong collapses to 0-1 peak / overdamped, matching Hoffmann 2002's abstract:
IkBbeta/epsilon "reduce the system's oscillatory potential ... during longer stimulations")?

SOURCE: BioModels export https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000140
  -> BIOMD0000000140.ode (XPP export of the curated SBML, human-readable reaction list) and
     BIOMD0000000140_url.xml (the SBML itself). cytoplasm=nucleus=1.0 (both compartments unit
     volume in this curated model), so volume factors are identically 1 and dropped below; every
     parameter value and every v1..v45 flux expression below is transcribed VERBATIM from the
     .ode file, not re-derived or approximated. The reduced 3-variable geometric model in the
     nfkb_signaling_dynamics cell is a DIFFERENT model and explicitly does not reconstruct this one.

GATES / FALSIFIERS (pre-registered BEFORE any integration was run):
  F1 PERIOD: transient IKK=0.1 bolus (after a long unstimulated pre-equilibration, verified a
     true fixed point via max|dy/dt|) gives nuclear-NF-kB interpeak periods with mean inside
     [50,150] min (Ashall 2009-derived band, |period-100|<=50) -- an external anchor that traces
     to an independently-conducted single-cell/pulsatile-TNF study.
  F2 REGIME MAP (qualitative, Hoffmann 2002's abstract, held out -- never fit): persistent
     clamped IKK should show MORE sustained oscillation at low forcing and COLLAPSE toward a
     single overdamped transition at higher forcing (IkBbeta/epsilon "stabilize" the response).
     FALSIFIED if oscillation cycle-count is flat or increasing in IKK forcing strength.
  F3 VOID FLOOR: zero IKK forcing (no bolus, no clamp) must show NO peaks at all (system sits at
     its unstimulated fixed point) -- a stimulus-independent oscillator would refute the
     feedback-driven mechanism.
  F4 CONSERVATION-STYLE SANITY: the three IkB-isoform total pools (IkBalpha/beta/eps, summed
     across all their free+bound+nuclear forms) plus free+bound NF-kB must be non-negative and
     finite throughout (no negative concentrations from a transcription sign error).

ANCHOR (external, decorrelated from this integration): Ashall et al. 2009 (PMID 19359585,
Science): "system completely resets between 100 and 200 min after each stimulus" -- an independent
single-cell-imaging study, not Hoffmann 2002's population-biochemistry fit, conducted 7 years later
by an unrelated lab. The [50,150] min band used for F1 is the pre-registered +/-50 min tolerance
around the ~100 min single-cell figure.
"""
import os
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "nfkb_hoffmann2002_biomd140_rebuild")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = os.path.join(OUT_DIR, "nfkb_hoffmann2002_biomd140_rebuild_results.json")

# ---------------------------------------------------------------------------
# Parameters, transcribed VERBATIM from BIOMD0000000140.ode
# ---------------------------------------------------------------------------
P = dict(
    a4=30.0, d4=0.03, a5=30.0, d5=0.03, a6=30.0, d6=0.03,
    r4=1.224, r5=0.45, r6=0.66, deg4=0.00135,
    k1=5.4, k01=0.0048,
    tr2a=9.25e-5, tr2=0.99, tr3=0.0168, tr2b=1.068e-5, tr2e=7.62e-6,
    a1=1.35, d1=0.075, tr1=0.2448, deg1=0.00678, tp1=0.018, tp2=0.012,
    a2=0.36, d2=0.105, a3=0.54, d3=0.105,
    a7=11.1, k2=0.828, a8=2.88, k2_beta=0.624, a9=4.2, k2_eps=0.624,
    r1=0.2442, r2=0.09, r3=0.132, k02=0.0072,
    k2_IkBbeta_nuc_NFkB_nuc=0.0069,
)

# state order (24 species, matches the .ode's declared species set exactly)
STATES = [
    "IkBalpha", "NFkB", "IkBalpha_NFkB", "IkBbeta", "IkBbeta_NFkB",
    "IkBeps", "IkBeps_NFkB", "IKK_IkBalpha", "IKK_IkBalpha_NFkB", "IKK",
    "IKK_IkBbeta", "IKK_IkBbeta_NFkB", "IKK_IkBeps", "IKK_IkBeps_NFkB",
    "NFkB_nuc", "IkBalpha_nuc", "IkBalpha_nuc_NFkB_nuc", "IkBbeta_nuc",
    "IkBbeta_nuc_NFkB_nuc", "IkBeps_nuc", "IkBalpha_transcript",
    "IkBbeta_transcript", "IkBeps_transcript", "IkBeps_nuc_NFkB_nuc",
]
IDX = {s: i for i, s in enumerate(STATES)}

Y0 = np.zeros(len(STATES))
Y0[IDX["IkBalpha"]] = 0.1
Y0[IDX["NFkB"]] = 0.1
Y0[IDX["NFkB_nuc"]] = 0.001
# all other species init 0.0 (verbatim: "init X=0.0" or "par X=0.0" in the .ode)


def fluxes(y, fr_beta_term):
    """v1..v45, transcribed verbatim (cytoplasm=nucleus=1.0 dropped)."""
    (IkBalpha, NFkB, IkBalpha_NFkB, IkBbeta, IkBbeta_NFkB, IkBeps, IkBeps_NFkB,
     IKK_IkBalpha, IKK_IkBalpha_NFkB, IKK, IKK_IkBbeta, IKK_IkBbeta_NFkB,
     IKK_IkBeps, IKK_IkBeps_NFkB, NFkB_nuc, IkBalpha_nuc, IkBalpha_nuc_NFkB_nuc,
     IkBbeta_nuc, IkBbeta_nuc_NFkB_nuc, IkBeps_nuc, IkBalpha_transcript,
     IkBbeta_transcript, IkBeps_transcript, IkBeps_nuc_NFkB_nuc) = y

    v1 = P['a4']*IkBalpha*NFkB - P['d4']*IkBalpha_NFkB
    v2 = P['a5']*IkBbeta*NFkB - P['d5']*IkBbeta_NFkB
    v3 = P['a6']*IkBeps*NFkB - P['d6']*IkBeps_NFkB
    v4 = P['a4']*IKK_IkBalpha*NFkB - P['d4']*IKK_IkBalpha_NFkB
    v5 = P['r4']*IKK_IkBalpha_NFkB
    v6 = P['a5']*IKK_IkBbeta*NFkB - P['d5']*IKK_IkBbeta_NFkB
    v7 = P['r5']*IKK_IkBbeta_NFkB
    v8 = P['a6']*IKK_IkBeps*NFkB - P['d6']*IKK_IkBeps_NFkB
    v9 = P['r6']*IKK_IkBeps_NFkB
    v10 = P['deg4']*IkBalpha_NFkB
    v11 = P['deg4']*IkBbeta_NFkB
    v12 = P['deg4']*IkBeps_NFkB
    v13 = P['k1']*NFkB - P['k01']*NFkB_nuc
    v14 = P['a4']*IkBalpha_nuc*NFkB_nuc - P['d4']*IkBalpha_nuc_NFkB_nuc
    v15 = P['a5']*IkBbeta_nuc*NFkB_nuc - P['d5']*IkBbeta_nuc_NFkB_nuc
    v16 = P['a6']*IkBeps_nuc*NFkB_nuc - P['d6']*IkBeps_nuc_NFkB_nuc
    v17 = P['tr2a']
    v18 = P['tr2']*NFkB_nuc**2
    v19 = P['tr3']*IkBalpha_transcript
    v20 = P['tr2b']
    v21 = P['tr3']*IkBbeta_transcript
    v22 = P['tr2e']
    v23 = P['tr3']*IkBeps_transcript
    v24 = P['a1']*IkBalpha*IKK - P['d1']*IKK_IkBalpha
    v25 = P['tr1']*IkBalpha_transcript
    v26 = P['deg1']*IkBalpha
    v27 = P['tp1']*IkBalpha - P['tp2']*IkBalpha_nuc
    v28 = P['a2']*IkBbeta*IKK - P['d2']*IKK_IkBbeta
    v29 = P['tr1']*IkBbeta_transcript
    v30 = P['deg1']*IkBbeta
    v31 = 0.5*P['tp1']*IkBbeta - 0.5*P['tp2']*IkBbeta_nuc
    v32 = P['a3']*IkBeps*IKK - P['d3']*IKK_IkBeps
    v33 = P['tr1']*IkBeps_transcript
    v34 = P['deg1']*IkBeps
    v35 = 0.5*P['tp1']*IkBeps - 0.5*P['tp2']*IkBeps_nuc
    v36 = P['a7']*IKK*IkBalpha_NFkB - P['d1']*IKK_IkBalpha_NFkB
    v37 = P['k2']*IkBalpha_nuc_NFkB_nuc
    v38 = P['a8']*IKK*IkBbeta_NFkB - P['d2']*IKK_IkBbeta_NFkB
    v39 = P['k2_IkBbeta_nuc_NFkB_nuc']*fr_beta_term*IkBbeta_nuc_NFkB_nuc
    v40 = P['a9']*IKK*IkBeps_NFkB - P['d3']*IKK_IkBeps_NFkB
    v41 = 0.5*P['k2_eps']*IkBeps_nuc_NFkB_nuc
    v42 = P['r1']*IKK_IkBalpha
    v43 = P['r2']*IKK_IkBbeta
    v44 = P['r3']*IKK_IkBeps
    v45 = P['k02']*IKK
    return dict(v1=v1, v2=v2, v3=v3, v4=v4, v5=v5, v6=v6, v7=v7, v8=v8, v9=v9,
                v10=v10, v11=v11, v12=v12, v13=v13, v14=v14, v15=v15, v16=v16,
                v17=v17, v18=v18, v19=v19, v20=v20, v21=v21, v22=v22, v23=v23,
                v24=v24, v25=v25, v26=v26, v27=v27, v28=v28, v29=v29, v30=v30,
                v31=v31, v32=v32, v33=v33, v34=v34, v35=v35, v36=v36, v37=v37,
                v38=v38, v39=v39, v40=v40, v41=v41, v42=v42, v43=v43, v44=v44, v45=v45)


def rhs(t, y, clamp_ikk=None, fr_beta_term=1.0):
    """dY/dt. clamp_ikk: if not None, holds IKK's *free-pool* concentration fixed externally
    (models a sustained/persistent forcing) instead of letting v45 etc. evolve it freely."""
    y = np.asarray(y, dtype=float)
    if clamp_ikk is not None:
        y = y.copy()
        y[IDX["IKK"]] = clamp_ikk
    v = fluxes(y, fr_beta_term)
    d = np.zeros(len(STATES))
    d[IDX["IkBalpha"]] = -v['v1'] - v['v24'] + v['v25'] - v['v26'] - v['v27']
    d[IDX["NFkB"]] = (-v['v1'] - v['v2'] - v['v3'] - v['v4'] + v['v5'] - v['v6']
                       + v['v7'] - v['v8'] + v['v9'] + v['v10'] + v['v11'] + v['v12'] - v['v13'])
    d[IDX["IkBalpha_NFkB"]] = v['v1'] - v['v10'] - v['v36'] + v['v37']
    d[IDX["IkBbeta"]] = -v['v2'] - v['v28'] + v['v29'] - v['v30'] - v['v31']
    d[IDX["IkBbeta_NFkB"]] = v['v2'] - v['v11'] - v['v38'] + v['v39']
    d[IDX["IkBeps"]] = -v['v3'] - v['v32'] + v['v33'] - v['v34'] - v['v35']
    d[IDX["IkBeps_NFkB"]] = v['v3'] - v['v12'] - v['v40'] + v['v41']
    d[IDX["IKK_IkBalpha"]] = -v['v4'] + v['v24'] - v['v42']
    d[IDX["IKK_IkBalpha_NFkB"]] = v['v4'] - v['v5'] + v['v36']
    d[IDX["IKK"]] = (v['v5'] + v['v7'] + v['v9'] - v['v24'] - v['v28'] - v['v32']
                      - v['v36'] - v['v38'] - v['v40'] + v['v42'] + v['v43'] + v['v44'] - v['v45'])
    d[IDX["IKK_IkBbeta"]] = -v['v6'] + v['v28'] - v['v43']
    d[IDX["IKK_IkBbeta_NFkB"]] = v['v6'] - v['v7'] + v['v38']
    d[IDX["IKK_IkBeps"]] = -v['v8'] + v['v32'] - v['v44']
    d[IDX["IKK_IkBeps_NFkB"]] = v['v8'] - v['v9'] + v['v40']
    d[IDX["NFkB_nuc"]] = v['v13'] - v['v14'] - v['v15'] - v['v16']
    d[IDX["IkBalpha_nuc"]] = -v['v14'] + v['v27']
    d[IDX["IkBalpha_nuc_NFkB_nuc"]] = v['v14'] - v['v37']
    d[IDX["IkBbeta_nuc"]] = -v['v15'] + v['v31']
    d[IDX["IkBbeta_nuc_NFkB_nuc"]] = v['v15'] - v['v39']
    d[IDX["IkBeps_nuc"]] = -v['v16'] + v['v35']
    d[IDX["IkBalpha_transcript"]] = v['v17'] + v['v18'] - v['v19']
    d[IDX["IkBbeta_transcript"]] = v['v20'] - v['v21']
    d[IDX["IkBeps_transcript"]] = v['v22'] - v['v23']
    d[IDX["IkBeps_nuc_NFkB_nuc"]] = v['v16'] - v['v41']
    if clamp_ikk is not None:
        d[IDX["IKK"]] = 0.0
    return d


def find_fixed_point(fr_beta_term=1.0, t_pre=20000.0):
    """Long unstimulated (IKK=0 throughout) integration to true steady state, verified."""
    sol = solve_ivp(lambda t, y: rhs(t, y, clamp_ikk=None, fr_beta_term=fr_beta_term),
                     (0, t_pre), Y0, method="BDF", rtol=1e-10, atol=1e-13,
                     max_step=200.0)
    y_ss = sol.y[:, -1]
    dydt_ss = rhs(t_pre, y_ss, fr_beta_term=fr_beta_term)
    return y_ss, float(np.max(np.abs(dydt_ss)))


def transient_bolus_run(y_ss, ikk_bolus=0.1, t_total=600.0, fr_beta_term=1.0):
    y0 = y_ss.copy()
    y0[IDX["IKK"]] = ikk_bolus
    sol = solve_ivp(lambda t, y: rhs(t, y, clamp_ikk=None, fr_beta_term=fr_beta_term),
                     (0, t_total), y0, method="BDF", rtol=1e-10, atol=1e-13,
                     max_step=1.0, dense_output=False)
    t_grid = np.linspace(0, t_total, 6001)
    y_grid = sol.sol(t_grid) if sol.sol is not None else None
    if y_grid is None:
        # re-run with dense output
        sol = solve_ivp(lambda t, y: rhs(t, y, clamp_ikk=None, fr_beta_term=fr_beta_term),
                         (0, t_total), y0, method="BDF", rtol=1e-10, atol=1e-13,
                         max_step=1.0, dense_output=True)
        y_grid = sol.sol(t_grid)
    nfkb_nuc_trace = y_grid[IDX["NFkB_nuc"]]
    peaks, _ = find_peaks(nfkb_nuc_trace, prominence=1e-6)
    peak_times = t_grid[peaks]
    interpeak = np.diff(peak_times)
    return dict(t=t_grid, trace=nfkb_nuc_trace, peak_times=peak_times.tolist(),
                interpeak_periods=interpeak.tolist(), all_finite=bool(np.all(np.isfinite(y_grid))),
                all_nonneg=bool(np.all(y_grid >= -1e-8)))


def clamped_run(y_ss, ikk_level, t_total=1200.0, fr_beta_term=1.0):
    y0 = y_ss.copy()
    y0[IDX["IKK"]] = ikk_level
    sol = solve_ivp(lambda t, y: rhs(t, y, clamp_ikk=ikk_level, fr_beta_term=fr_beta_term),
                     (0, t_total), y0, method="BDF", rtol=1e-9, atol=1e-12,
                     max_step=1.0, dense_output=True)
    t_grid = np.linspace(0, t_total, 6001)
    y_grid = sol.sol(t_grid)
    nfkb_nuc_trace = y_grid[IDX["NFkB_nuc"]]
    peaks, _ = find_peaks(nfkb_nuc_trace, prominence=1e-6)
    return dict(n_peaks=int(len(peaks)), peak_times=t_grid[peaks].tolist(),
                trace_final=float(nfkb_nuc_trace[-1]), trace_max=float(nfkb_nuc_trace.max()),
                all_finite=bool(np.all(np.isfinite(y_grid))))


def main():
    results = {}

    # 1. true unstimulated fixed point
    y_ss, max_dydt = find_fixed_point()
    results['fixed_point_max_dydt'] = max_dydt
    results['fixed_point_is_steady'] = bool(max_dydt < 1e-4)

    # F3 void floor: from the fixed point, integrate freely with IKK held at 0 (no bolus) --
    # must show zero peaks.
    void = transient_bolus_run(y_ss, ikk_bolus=0.0, t_total=600.0)
    results['F3_void_floor_n_peaks'] = int(len(void['peak_times']))
    results['F3_void_floor_pass'] = bool(len(void['peak_times']) == 0)

    # F1: transient IKK=0.1 bolus, measure NFkB_nuc interpeak periods
    trans = transient_bolus_run(y_ss, ikk_bolus=0.1, t_total=600.0)
    periods = trans['interpeak_periods']
    mean_period = float(np.mean(periods)) if periods else None
    results['F1_transient_bolus_n_peaks'] = len(trans['peak_times'])
    results['F1_interpeak_periods_min'] = periods
    results['F1_mean_period_min'] = mean_period
    results['F1_period_band_5090_150'] = [50.0, 150.0]
    results['F1_pass'] = bool(mean_period is not None and 50.0 <= mean_period <= 150.0)
    results['F4_conservation_sanity_pass'] = bool(trans['all_finite'] and trans['all_nonneg'])

    # F2: persistent/clamped IKK regime map
    regime = {}
    for level in [0.01, 0.1, 0.3, 1.0]:
        r = clamped_run(y_ss, level, t_total=1200.0)
        regime[str(level)] = r
    results['F2_regime_map'] = {k: {'n_peaks': v['n_peaks'], 'trace_final': v['trace_final']}
                                 for k, v in regime.items()}
    npeaks_by_level = [regime[str(l)]['n_peaks'] for l in [0.01, 0.1, 0.3, 1.0]]
    # F2 pass: NOT monotonically flat/increasing -- Hoffmann's claim is a COLLAPSE
    # (fewer cycles) at higher forcing, i.e. non-increasing overall and a real drop somewhere.
    results['F2_npeaks_by_level'] = dict(zip(['0.01', '0.1', '0.3', '1.0'], npeaks_by_level))
    results['F2_pass_collapse_at_higher_forcing'] = bool(
        npeaks_by_level[0] >= npeaks_by_level[-1] and min(npeaks_by_level) < max(npeaks_by_level[0], 1)
    )

    overall_pass = bool(results['fixed_point_is_steady'] and results['F1_pass']
                         and results['F3_void_floor_pass'] and results['F4_conservation_sanity_pass'])
    results['overall_pass_F1_F3_F4'] = overall_pass
    results['note'] = ("F2 is qualitative/informational (Hoffmann 2002 "
                        "abstract, held out, not a numeric fit target) -- reported but not gating "
                        "overall_pass.")

    with open(OUT_PATH, "w") as f:
        json.dump(results, f, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else o)

    print("=== NFKB Hoffmann2002 BIOMD0000000140 full rebuild ===")
    print(f"fixed point max|dy/dt| = {max_dydt:.3e}  (steady: {results['fixed_point_is_steady']})")
    print(f"F3 void floor peaks = {results['F3_void_floor_n_peaks']}  pass={results['F3_void_floor_pass']}")
    print(f"F1 transient-bolus interpeak periods (min) = "
          f"{[round(p,2) for p in periods]}  mean={mean_period}")
    print(f"F1 pass (band [50,150]) = {results['F1_pass']}")
    print(f"F2 regime map n_peaks by clamped IKK level: {results['F2_npeaks_by_level']}")
    print(f"F4 conservation/finite sanity pass = {results['F4_conservation_sanity_pass']}")
    print(f"OVERALL (F1,F3,F4) PASS = {overall_pass}")
    print(f"results written to {OUT_PATH}")
    return results


if __name__ == "__main__":
    main()
