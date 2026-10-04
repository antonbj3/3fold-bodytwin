from dental_release.paths import expand as _release_expand
import os
import json, hashlib, time
from pathlib import Path
import numpy as np
from calibration import fit_affine, apply, bone_branch_b
P = Path(__file__).resolve().parent
OUT = Path(os.environ.get('X24_OUTPUT_DIR', str(P)))
OUT.mkdir(parents=True, exist_ok=True)
pre = json.loads((P / 'PREREG_R1.json').read_text())
assert hashlib.sha256((P / 'PREREG_R1.json').read_bytes()).hexdigest() == (P / 'PREREG_R1.sha256').read_text().strip()
R = json.loads((P / 'PHANTOM_BASELINE.json').read_text())
old = json.loads(Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/bodytwin/artifacts/cbct_multidevice_sensitometry.json')).read_text())
ref = next((r for r in R if r['device'] == 'Varian'))
V = {r['material']: r['gray_median'] for r in ref['inserts']}
rows = []
replay = []
prop = []
for dev in R:
    g = {r['material']: r['gray_median'] for r in dev['inserts']}
    anc = pre['anchors']
    x = np.array([g[m] for m in anc])
    y = np.array([V[m] for m in anc])
    c = fit_affine(x, y)
    ctrl = np.linalg.lstsq(np.column_stack([np.ones(2), x]), y, rcond=None)[0]
    for m in pre['held_out']:
        h = float(apply(c, g[m]))
        hc = float(apply(ctrl, g[m]))
        row = {'device': dev['device'], 'material': m, 'g': g[m], 'H_reference': V[m], 'H_calibrated': h, 'error_practice_HU': g[m] - V[m], 'error_calibrated_HU': h - V[m], 'control_disagreement_HU': abs(h - hc), 'pass': abs(h - V[m]) <= 40, 'resolution': 'PER_SURFACE_REGION'}
        if m == 'Delrin':
            row['injected_100gray_error_HU'] = float(apply(c, g[m] + 100) - V[m])
            row['injected_fail'] = abs(row['injected_100gray_error_HU']) > 40
        rows.append(row)
    legacy_rho = {'AIR': 0.001, 'PMP': 0.853, 'LDPE': 0.945, 'Polystyrene': 1.017}
    xx = np.array(list(legacy_rho.values()))
    yy = np.array([g[m] for m in legacy_rho])
    slope = float(np.linalg.lstsq(np.column_stack([xx, np.ones(4)]), yy, rcond=None)[0][0])
    o = old['devices'][dev['device']]
    replay.append({'device': dev['device'], 'air': g['AIR'], 'legacy_tissue_slope': slope, 'air_replay_error': g['AIR'] - o['air'], 'slope_replay_error': slope - o['slope_tissue'], 'a_inverse_Href': c.tolist(), 'gain_g_per_Href': 1 / c[1], 'water_g_affine_extrapolation': float(-c[0] / c[1]), 'water_status': 'UNMEASURED_AFFINE_EXTRAPOLATION', 'resolution': 'PER_SURFACE_REGION'})
    for h in pre['propagation']['evaluation_HU']:
        gs = float((h - c[0]) / c[1])
        e = bone_branch_b(h)
        bad = bone_branch_b(gs)
        prop.append({'device': dev['device'], 'H_reference_scenario': h, 'gray_affine_scenario': gs, 'E_reference_MPa': float(e['E_MPa']), 'E_naive_MPa': float(bad['E_MPa']), 'E_naive_ratio': float(bad['E_MPa'] / e['E_MPa']), 'density_naive_valid': bool(bad['density_valid']), 'resolution': 'PHENOMENOLOGICAL', 'status': 'conditional scenario from observed affine scanner transform, not measured bone', 'replacement_measurement': 'paired HA phantom plus mandibular density/fabric/compression'})
ele = [r for r in rows if r['device'] != 'Varian']
improvement = float(np.median([abs(r['error_practice_HU']) for r in ele]) / np.median([abs(r['error_calibrated_HU']) for r in ele]))
calmax = max((abs(r['error_calibrated_HU']) for r in ele))
checkmax = max((r['control_disagreement_HU'] for r in rows))
gates = {'holdout_40HU': calmax <= 40, 'improvement_4x': improvement >= 4, 'equal_information_check': checkmax <= 1e-08, 'legacy_measurement_replay': max((max(abs(r['air_replay_error']), abs(r['slope_replay_error'])) for r in replay)) <= 1, 'injected_error_rejected': all((r['injected_fail'] for r in ele if r['material'] == 'Delrin'))}
out = {'round': 'R1', 'claim_type': 'information_link', 'outcome': 'PASS_REFERENCE_TRANSFER' if all(gates.values()) else 'NEGATIVE_TWO_ANCHOR_SUFFICIENCY', 'prereg_sha256': hashlib.sha256((P / 'PREREG_R1.json').read_bytes()).hexdigest(), 'external_referent': pre['external_referent'], 'gates': gates, 'max_heldout_error_HU': calmax, 'median_error_improvement_factor': improvement, 'equal_information_max_disagreement_HU': checkmax, 'rows': rows, 'scanner_calibration': replay, 'conditional_E_propagation': prop, 'dropout': {'devices_eligible': 4, 'devices_rejected': 0, 'material_ROIs_eligible': 24, 'material_ROIs_rejected': 0, 'reason': 'none; 12 held-out Elekta/material comparisons, 4 Varian identity checks excluded from efficacy denominator'}, 'limitations': ['reference Varian not absolute density truth', 'no water measured', 'polymer response not HA/bone response', 'ROI medians not pointwise calibration', 'Keller and NewTom laws are separate closures', 'current K37 internally normalizes gray, so raw-HU error is counterfactual practice branch, not demonstrated current-code bug']}
(OUT / 'results_R1.json').write_text(json.dumps(out, indent=2) + '\n')
(OUT / 'HANDOFF_R1.md').write_text(f"# R1 : two references tested against external phantom\n\nUtfall {out['outcome']}. Max holding error {calmax:.3f} HU_ref, median improvement {improvement:.3f}×. Gates {gates}.\n\nNext construction: move the second anchor to the leg's attaching area; freeze R2 before the next numerical. Luft/vatten empiriskt UNKNOWN; materiallagen separat.\n")
s = json.loads((P / 'CURRENT_WORK_STATE.json').read_text())
s.update(phase='R1_COMPLETE', latest_gate=gates, next_operation='freeze R2 bone-range calibration with independent protocol/retest holdout')
(OUT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(s, indent=2) + '\n')
print(out['outcome'], gates, 'max', calmax, 'improvement', improvement)
print(replay)
