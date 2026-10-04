import hashlib, json, time
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1]

def write(p, x):
    if Path(p).name == 'CURRENT_WORK_STATE.json':
        x = dict(x, lane='X82-height-to-force', updated_at_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat())
        with (R / 'raw/MILESTONES.jsonl').open('a') as stream:
            stream.write(json.dumps(x) + '\n')
    dest = R / p
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.name == 'CURRENT_WORK_STATE.json':
        temp = dest.with_suffix('.tmp')
        temp.write_text(json.dumps(x, indent=2, allow_nan=False))
        temp.replace(dest)
    else:
        dest.write_text(json.dumps(x, indent=2, allow_nan=False))

def main():
    t = time.perf_counter()
    p = R / 'PREREG_R2.json'
    assert hashlib.sha256(p.read_bytes()).hexdigest() == (R / 'PREREG_R2.json.sha256').read_text().split()[0]
    s = json.loads((R / 'sources/CHO2014_NUMBERS.json').read_text())
    groups = []
    rows = []
    for food in s['measurements']:
        for (j, h) in enumerate(s['height_mm']):
            means = np.array(food['means_N'][j])
            total = float(means.sum())
            T = food['applied_N']
            admissible = bool(total <= T + 1)
            groups.append({'food': food['food'], 'height_mm': h, 'table': food['table'], 'published_signal_sum_N': total, 'applied_load_N': T, 'excess_above_applied_N': total - T, 'necessary_partial_axial_sum_pass': admissible, 'absolute_force_port_status': 'UNKNOWN_OBSERVATION_LAW_AND_LOAD_PATH'})
            for (tooth, mean, sd) in zip(s['tooth_order'], means, food['sd_N'][j]):
                rows.append({'food': food['food'], 'height_mm': h, 'fdi': tooth, 'signal_N': float(mean), 'sd_N': sd, 'n': 20, 'resolution': 'PER_TOOTH', 'table': food['table'], 'unit': 'N equivalent strain signal', 'locator': s['doi'], 'absolute_reaction_eligible': False, 'reason': 'Observation map/load-path not matched; SD is not instrument uncertainty'})
    delta = []
    for f in s['measurements']:
        x = np.array(f['means_N'])
        delta.append({'food': f['food'], 'height_change_mm': 0.023, 'baseline_height_mm': 0.0, 'FDI': s['tooth_order'], 'measured_signal_change_N': (x[3] - x[2]).tolist(), 'uncertainty': 'UNKNOWN paired difference uncertainty; correlation/raw strokes unavailable', 'quantity': 'Change in strain-calibrated signal, not certified reaction force'})
    M1 = np.eye(2)
    M2 = np.diag([0.5, 1.5])
    f1 = np.array([10.0, 10.0])
    f2 = np.array([15.0, 5.0])
    b2 = f1 - M2 @ f2
    y1 = M1 @ f1
    y2 = M2 @ f2 + b2
    suff = {'same_sensor_signal_N': [y1.tolist(), y2.tolist()], 'same_total_force_N': [float(f1.sum()), float(f2.sum())], 'signal_identity_error_N': float(np.max(np.abs(y1 - y2))), 'total_identity_error_N': float(abs(f1.sum() - f2.sum())), 'bit_identical': bool(y1.tobytes() == y2.tobytes()), 'downstream_tooth_difference_N': float(abs(f1[0] - f2[0])), 'minimum_extension': 'Known per-channel observation gain/offset (or direct reaction load cells) at the query contact direction and time.', 'worlds': [{'M': M1.tolist(), 'b': [0, 0], 'true_force_N': f1.tolist()}, {'M': M2.tolist(), 'b': b2.tolist(), 'true_force_N': f2.tolist()}], 'kind': 'our_own_fixture; exact algebraic ambiguity, no empirical calibration'}
    loads = np.array([0.0, 5.0, 10.0, 20.0])
    M = np.array([[1.2, 0.1], [0.05, 0.8]])
    b = np.array([0.25, -0.5])
    F = np.array([[u, v] for u in loads for v in loads])
    Y = F @ M.T + b
    fit = np.linalg.lstsq(np.c_[F, np.ones(len(F))], Y, rcond=None)[0]
    held = np.array([12.0, 7.0])
    y = M @ held + b
    recovered = np.linalg.solve(fit[:2].T, y - fit[2])
    err = float(np.max(np.abs(recovered - held)))
    corrupt = fit[:2].T.copy()
    corrupt[0] *= 0.8
    bad = np.linalg.solve(corrupt, y - fit[2])
    bad_error = float(abs(bad.sum() - held.sum()))
    calibration = {'kind': 'our_own_fixture', 'fit_measurements': len(F), 'held_out_force_N': held.tolist(), 'recovered_force_N': recovered.tolist(), 'max_error_N': err, 'injected_gain_bias': 0.2, 'injected_total_balance_error_N': bad_error, 'fault_detected': bool(bad_error > 1), 'empirical_sensor_calibration': 'NOT_RUN'}
    out = {'round': 'R2', 'claim_type': 'capability', 'external_referent': json.loads(p.read_text())['external_referent'], 'prereg_sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'source_sha256': hashlib.sha256((R / 'sources/CHO2014_NUMBERS.json').read_bytes()).hexdigest(), 'external_groups': len(groups), 'necessary_partial_sum_failures': sum((not x['necessary_partial_axial_sum_pass'] for x in groups)), 'groups': groups, 'signal_deltas': delta, 'rows_considered': len(rows), 'rows_numeric_retained': len(rows), 'rows_admitted_as_absolute_reaction': 0, 'attrition_absolute_reaction_fraction': 1.0, 'attrition_reasons': {'unknown_observation_operator_and_load_path': len(rows), 'missing_unit_or_locator': 0}, 'observation_sufficiency': suff, 'calibration_fixture': calibration, 'gates': {'external_absolute_reaction_transfer': False, 'synthetic_calibration': err <= 1e-07, 'gain_fault_detected': calibration['fault_detected']}, 'empirical_height_to_force': 'UNKNOWN; external experiment is relevant but reaction port incompatible', 'cost': {'wall_s': time.perf_counter() - t, 'physical_acquisition': 'NOT_RUN', 'external_fit': 'NONE; no patient/support fit performed', 'source_screening': 'See sources/SEARCH_LEDGER.json'}}
    write('raw/R2_OBSERVATION_ROWS.json', rows)
    write('rounds/R2/results.json', out)
    write('CURRENT_WORK_STATE.json', {'milestone': 'R2_EXTERNAL_PORT_REFUTED', 'latest_gate': out['gates'], 'next_operation': 'Minimal extra load-response acquisition for the R1 new-wrench counterexample'})
    print(json.dumps({k: out[k] for k in ['external_groups', 'necessary_partial_sum_failures', 'rows_admitted_as_absolute_reaction', 'gates']}, indent=2))
if __name__ == '__main__':
    main()
