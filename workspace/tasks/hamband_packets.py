"""BT-HB: one The swarm packet per L1 gait trial (113). Two tests per packet on real data:
(1) an independent implementation of the ε-band (contact under stress-2 ≤ (1+ε)·min) — coverage/width against measured implant force;
(2) does the measured force's excess over the LP minimum (meas−lo) follow MEASURED hamstring EMG over time (U447's structural claim)?
Real files only (no symlinks). The L1 selection contains no jw_lungef1."""
import json, shutil
from pathlib import Path
W = Path(''); R = W / 'results'
loc = json.load(open(R / 'CX-INVERSEOC/locations.json'))['trials']
q = []; prof = 'ABC'
for k, (key, t) in enumerate(sorted(loc.items())):
    npz = R / 'L1/prep' / f'{key}.npz'
    if not npz.is_file() or 'lungef' in key: continue
    j = 'BT-HB-' + key.replace('__', '-').replace('_', '-')
    d = R / j; inp = d / 'inputs'
    if (d / 'RESULTS.md').exists(): continue
    inp.mkdir(parents=True, exist_ok=True)
    shutil.copy2(npz, inp / 'trial.npz')
    json.dump({'key': key, **{x: t[x] for x in ('person', 'activity', 'group', 'n_frames')}, 'rows': t['rows']},
              open(inp / 'lp_bounds_inverseoc.json', 'w'))
    shutil.copy2(R / 'CX-INVERSEOC/epsilon.py', inp / 'ref_inverseoc_epsilon.py')
    shutil.copy2(W / 'tasks/NIGHT_PREAMBLE.md', inp / 'NIGHT_PREAMBLE.md')
    (d / 'BRIEF.md').write_text(f"""# {j} — ε-band + hamstring co-contraction on {t['person']} {t['activity']} (L1 trial {key})
Read inputs/NIGHT_PREAMBLE.md. Work only in this directory, 1 thread, < 1 GB, total ≤ 20 min. All data is in inputs/ (real files; read nothing outside).

## Data (inputs/trial.npz, per frame i, n={t['n_frames']})
- A[i] (6×166) f = b[i]: equilibrium; 0 ≤ f ≤ F0 (166 muscle elements).
- total knee contact = C0[i] + cj[i]·f.
- medial = 0.5·C0 − mx/d + (0.5·cj − Mkx/d)·f ≥ 0; lateral = 0.5·C0 + mx/d + (0.5·cj + Mkx/d)·f ≥ 0 (d scalar).
- meas_tot[i] = measured implant force (N); bwN = body weight (N); inr = frames the source analysis accepted.
- emg[i] (15 channels, normalised envelope, order: semimem, bifem, vasmed, vaslat, rf, medgas, latgas, tfl, tibant, peronl, soleus, addmagnus, gmax, gmed, sartorius).
- PCSA_i = F0_i/60 N/cm²; stress-2 E(f) = Σ(f_i/PCSA_i)².
- inputs/lp_bounds_inverseoc.json holds an earlier analysis's lo_BW/hi_BW/meas_BW per frame, for CROSS-CHECKING your own LP only.
- inputs/ref_inverseoc_epsilon.py shows how the earlier analysis did this. Write your OWN code; you may read the reference.

## Frozen criteria (write PREREG.md + PREREG.sha256 BEFORE computing)
1. **The LP:** per frame, lo/hi = min/max of total contact under equilibrium, box and med/lat ≥ 0 (scipy linprog highs). Cross-check against lp_bounds_inverseoc.json (report the maximum deviation in BW).
2. **The ε-band:** E_min = min E(f) over the same set (a QP; SLSQP or trust-constr is fine; report solver status). For ε ∈ {{0.05, 0.10, 0.20}}: [lo_ε, hi_ε] = min/max total contact under the same constraints + E(f) ≤ (1+ε)·E_min. Report per ε:
   - (a) the share of frames where meas lies in [lo_ε, hi_ε], over ALL LP-feasible frames and over those with meas ≥ lo;
   - (b) median width in BW;
   - (c) the band midpoint's RMSE in BW against meas, and against N1g = 2.1153·‖grf‖ (grf in N? check the magnitude against bwN and state the unit).
   Criterion (frozen): at ε=0.10, coverage ≥ 90 % of frames with meas ≥ lo, AND median width ≤ 0.8 BW.
3. **Hamstrings:** r = the Pearson correlation between (meas − lo) and the hamstring EMG (mean of semimem+bifem), over LP-feasible frames. Compare against:
   - quadriceps (vasmed+vaslat+rf), gastroc (medgas+latgas), and GRF magnitude;
   - a lag ±100 ms, reporting the lag with the maximum r;
   - a control: hamstring EMG circularly shifted by half the trial (r must fall).
   Criterion (frozen): the hamstring r > the quadriceps r AND > the gastroc r, AND r > r_shift + 0.2.
4. Solve at most 3 ε values × 2 per frame; if time is short, take every other frame and say so.

Every outcome is a data point: report the numbers and which node the deviation points to (C0/geometry if meas < lo, EMG normalisation, phase). No verdict words.
Deliver results.json (all numbers above) and RESULTS.md starting with the line "{j}".
""")
    q.append(f'{prof[len(q) % 3]} swarm {j}')
open(W / 'tasks/lanes/bt_queue.txt', 'a').write('\n'.join(q) + '\n')
print(len(q), 'BT-HB packets queued')
