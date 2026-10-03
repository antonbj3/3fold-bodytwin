"""Test of breakthrough criterion 1 (coordinator 24/9 ~22:00): does the parameter-free knee law hold for OTHER persons and activities
than it was developed on (A1391/A1398)? One packet per trial with raw implant force + GRF + markers + the person's static trial + the ID code.
Real files only (no symlinks). Excludes jw_lungef1 (F-8)."""
import csv, os, shutil, glob
W = '.'; R = W + '/results'
rows = list(csv.DictReader(open(R + '/DATAMATRIX_TABLES/P-CONTACT-COMP.csv')))
def dm(ds, unit, proto):
    d = f"{R}/BT-DM-{ds}-{unit.replace('_', '-')}-{proto}/inputs/raw.csv"
    return d if os.path.isfile(d) else None
statics = {}
for d in glob.glob(R + '/BT-DM-*-P-MARKER-GEOM'):
    b = os.path.basename(d)
    if 'static' in b.lower():
        parts = b.split('-'); ds, person = parts[2], parts[3].upper()
        statics.setdefault((ds, person), d + '/inputs/raw.csv')
q = []; prof = 'ABC'; k = 0
for r in rows:
    ds, u, p, a = r['dataset'], r['unit'], r['person'], r['activity']
    if 'lungef' in u or p == 'JW' and 'lunge' in a: continue
    g, m, c = dm(ds, u, 'P-GRF-PEAK'), dm(ds, u, 'P-IMU-VIRTUAL'), dm(ds, u, 'P-CONTACT-COMP')
    if not (g and m and c): continue
    j = f"BT-LG-{ds}-{u.replace('_', '-')}"; d = f'{R}/{j}'
    if os.path.exists(d + '/RESULTS.md'): continue
    os.makedirs(d + '/inputs', exist_ok=True)
    shutil.copy2(g, d + '/inputs/grf.csv'); shutil.copy2(m, d + '/inputs/markers.csv'); shutil.copy2(c, d + '/inputs/implant_etibia_raw.csv')
    st = statics.get((ds, p.upper()))
    if st: shutil.copy2(st, d + '/inputs/static_markers.csv')
    for f in ('build.py', 'README_FOR_FIELD.md'):
        s = f'{R}/CX-LUNGEID/{f}'
        if os.path.isfile(s): shutil.copy2(s, d + '/inputs/lungeid_' + f)
    for f in ('RESULTS.md',):
        shutil.copy2(f'{R}/CX-LAW2/{f}', d + '/inputs/LAW2_RESULTS.md')
        shutil.copy2(f'{R}/CX-CONTACTDEF/{f}', d + '/inputs/CONTACTDEF_RESULTS.md')
    shutil.copy2(W + '/tasks/NIGHT_PREAMBLE.md', d + '/inputs/NIGHT_PREAMBLE.md')
    open(d + '/BRIEF.md', 'w').write(f'# {j} — keeps the parameter-free knee law for {p} i aktiviteten {a}? (genombrottskriterium 1)\nRead inputs/NIGHT_PREAMBLE.md. Only work in this directory; ≤ 60 s per run, 1 thread, < 1 GB. All files are in inputs/ (real files; read nothing outside).\nThe Law (A1398, CX-LAW2; read inputs/LAW2_RESULTS.md): F_law = ||GRF|| + |M_knee|/r_q + 0,29·|M_ankle|/0,045, r_q = 0,045 m (literature tape 0,030–0,060). NO adaptation to the implant force.\nData: inputs/grf.csv (GRF in N; select plate under implant bone and justify), inputs/markers.csv (marker trajectories), inputs/static_markers.csv (if available; static calibration), inputs/implant_etibia_raw.csv (eTibia in POUND-FORCE: multiply by 4,4482216152605 → N; total force = norm(Fx,Fy,Fz) or Fz according to inputs/CONTACTDEF_RESULTS.md).\nSteps: (1) calculate net sagittal moment at knee and ankle from markers+GRF (quasi-static or Newton–Euler; feel free to reuse the ideas in inputs/lungeid_build.py; report segment definitions and characters). (2) F_law per ruta. (3) Compare to measured implant force: top (BW and %), time series-RMSE, korrelation, k_law = median(F_law/||GRF||) versus measured k = median(F_meas/||GRF||) there GRF ≥ 0,2 BW. (4) Compare with walk-N1g k = 2,81 (A1390). (5) Sensitivity r_q 0,030/0,045/0,060.\nThe outcome is a data point: report the numbers and which node the deviation points to (arm curve in deep flexion, co-contraction, dynamics, error platta/bone). No judgement.\nPREREG.md + PREREG.sha256 before calculation. Deliver results.json and RESULTS.md starting with the line "{j}".\n')
    q.append(f'{prof[k % 3]} swarm {j}'); k += 1
open(W + '/tasks/lanes/bt_queue.txt', 'a').write('\n'.join(q) + '\n')
print(len(q), "layer generalization package queued")
