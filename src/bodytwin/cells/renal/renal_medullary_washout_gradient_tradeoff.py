"""Medullary washout vs corticopapillary gradient trade-off probe.

Sweeps the vasa-recta washout term of the countercurrent_multiplier cell's validated simulator
(washout used as a medullary-blood-flow proxy) and reports the papillary osmolality at each value,
the local d(papilla)/d(washout) slope around the human-calibrated point W_CAL=0.03, and the effect
of doubling washout (0.03 -> 0.06).

Reads: the countercurrent_multiplier module (imported, not read as data). Writes: nothing --
this is a print-only probe; the printed slope and the 2x-washout comparison are the result.
"""
import os, json
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))
_CCM_PATH = os.path.join(os.path.dirname(_HERE), "organ_systems", "countercurrent_multiplier.py")
spec = importlib.util.spec_from_file_location("ccm", _CCM_PATH)
ccm = importlib.util.module_from_spec(spec); spec.loader.exec_module(ccm)

C0=ccm.C0; DC=ccm.DC_SINGLE_EFFECT; N=ccm.N_CAL; NI=ccm.N_ITER
# washout as a medullary-blood-flow proxy (vasa recta washout term already in the validated model)
washouts = [0.005,0.01,0.02,0.03,0.04,0.06,0.08,0.10,0.15,0.20]
rows=[]
for w in washouts:
    r = ccm.run_ccm(N, DC, C0, w, NI, topology="counter")
    rows.append((w, round(r["papilla"],1), r["converged"]))
for w,p,c in rows:
    print(w,p,c)

# local slope d(papilla)/d(washout) around the calibrated human point W_CAL=0.03
import numpy as np
w_arr=np.array([r[0] for r in rows]); p_arr=np.array([r[1] for r in rows])
# finite-difference slope bracketing 0.03
idx = list(w_arr).index(0.03)
slope_up = (p_arr[idx+1]-p_arr[idx])/(w_arr[idx+1]-w_arr[idx])
slope_down = (p_arr[idx]-p_arr[idx-1])/(w_arr[idx]-w_arr[idx-1])
print("slope_around_calibrated(mOsm per unit washout, up/down):", slope_up, slope_down)
# doubling washout from calibrated 0.03->0.06 (proxy: medullary flow doubles)
r_double = ccm.run_ccm(N, DC, C0, 0.06, NI, topology="counter")
print("papilla at 2x washout (0.06):", r_double["papilla"], "vs calibrated 0.03:", ccm.run_ccm(N,DC,C0,0.03,NI,topology="counter")["papilla"])
