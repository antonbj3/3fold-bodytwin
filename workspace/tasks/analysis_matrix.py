"""Analysis matrix (coordinator 24/9 ~21:10): outcome × signal family × validation design → Space The swarm analysis packets.
Deterministic, no LLM per packet. Every packet: a real held-out prediction question plus a permutation null, on the CORRECTED tables
(CX-CONTACTDEF A1390: eTibia lbf→N, GRF per single plate). Consumers: sensor models (A321), N1g/activity prior, CX-BAYESOED priors.
F-8's jw_lungef1 is not in the tables (CX-DMCOMPUTE)."""
import os, shutil, itertools
W = '.'; T = W + '/results/DATAMATRIX_TABLES'; R = W + '/results'
OUTCOMES = {'kontakttopp': 'contact_peak_N', 'kontaktimpuls': 'contact_impulse_Ns', 'kontaktmedel': 'contact_mean_N'}
FAMILIES = {'GRF-topp': ['P-GRF-PEAK.csv'], 'GRF-impuls': ['P-GRF-IMPULSE.csv'], 'GRF-balans': ['P-GRF-BALANCE.csv'],
            'EMG-timing': ['P-EMG-TIMING.csv'], 'EMG-amplitud': ['P-EMG-AMPLITUDE.csv'],
            "marker_geometry": ['P-MARKER-GEOM.csv'], 'virtuell-IMU': ['P-IMU-VIRTUAL.csv']}
DESIGNS = {'LOPO': "leave one person out at a time (unit = person; gait cycles/trials from the same person are NOT independent)",
           'LOAO': "leave one activity out at a time (generalisation to an unseen activity)",
           'LODO': "leave one dataset (GC competition) out at a time (generalisation between protocols)"}
queue = []; prof = 'ABC'; k = 0
for (on, oc), (fn, ff), (dn, dd) in itertools.product(OUTCOMES.items(), FAMILIES.items(), DESIGNS.items()):
    j = f'BT-AM-{on}-{fn}-{dn}'.replace('ö', 'o').replace('å', 'a').replace('ä', 'a')
    d = f'{R}/{j}'
    if os.path.exists(f'{d}/RESULTS.md'): continue
    os.makedirs(d + '/inputs', exist_ok=True)
    for f in set(['P-CONTACT-COMP.csv', 'P-GRF-PEAK.csv'] + ff):
        shutil.copy2(f'{T}/{f}', f'{d}/inputs/{f}')
    shutil.copy2(W + '/tasks/NIGHT_PREAMBLE.md', d + '/inputs/NIGHT_PREAMBLE.md')
    open(d + '/BRIEF.md', 'w').write(f'''# {j} — add {fn} information on measured knee contact ({oc}) other than: |GRF|?\nRead inputs/NIGHT_PREAMBLE.md. Work only in this directory; ≤ 60 s per run, 1 thread, < 1 GB. Everything you need is in inputs/ (corrected tables, CX-CONTACTDEF: contact forces in N, GRF = largest single plate).\nUtfall: `{oc}` i P-CONTACT-COMP.csv (measured implant force, Grand Challenge). Join tables on `packet`-prefix/`dataset`+`unit` (same trial); report the number of matched attempts per person/aktivitet.\nBaselines (same validation): (1) training median, (2) k·peak_resultant_N through the origin (N1-like), (3) k per activity (N1g-like).\nQuestion: improves the addition of {fn}- The fields. ({', '.join(ff)}) prediction beyond baseline 3? Validering: {dd}.\nMotprov: permutera {fn}-the fields inside the person (eller inom aktivitet) 200 administration → win must disappear. Report full-out MAE/RMSE (N and %BW if body weight is available), vinst mot baslinje 3 with uncertainty (bootstrap on PERSONER), permutations-p.\nThe outcome is a data point that expands a node: which signals carry load information, under which generalisation. No verdict words; numbers + which node (sensor model/N1g prior/BAYESOED prior) is affected.\nPREREG.md + PREREG.sha256 before calculation. Deliver results.json and RESULTS.md starting with the line "{j}".\n''')
    queue.append(f'{prof[k % 3]} swarm {j}'); k += 1
open(W + '/tasks/lanes/bt_queue.txt', 'a').write('\n'.join(queue) + '\n')
print(len(queue), "analysis packets queued")
