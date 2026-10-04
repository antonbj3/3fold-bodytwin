import datetime, hashlib, json, time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from crown_corridor import predict, inverse
from run_r0_r1 import conventional
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
    pr = R / 'PREREG_R4.json'
    sha = hashlib.sha256(pr.read_bytes()).hexdigest()
    assert sha == (R / 'PREREG_R4.json.sha256').read_text().split()[0]
    cfg = json.loads(pr.read_text())
    m = cfg['metrics']
    s = m['probe_height_mm']
    eps = m['raw_per_channel_error_N']
    rows = []
    inverses = []
    calibrations = []
    states = []
    fullcount = 0
    for (ix, p) in enumerate(sorted((R / 'inputs').glob('*_R4_state.json'))):
        old = json.loads(p.read_text())
        f0 = np.array(old['baseline_force_N'])
        A = np.array(old['A'])
        N = np.array(old['basis_N'])
        H = np.array(old['H_N_per_mm'])
        i = old['predicted_fdi'].index(16) if 16 in old['predicted_fdi'] else 0
        n = len(f0)
        j = N @ H @ N[i]
        P = N @ N.T
        observations = [f0 + eps * 0.6 * np.sin(np.arange(n) + ix), f0 - j * s + eps * 0.8 * np.sin(np.arange(n) + ix + 1), f0 + j * s + eps * 0.7 * np.cos(np.arange(n) + ix + 2)]
        raw0 = observations[0]
        b = raw0 - A @ np.linalg.solve(A.T @ A, A.T @ raw0 - np.array(old['w_N']))
        jj = P @ (observations[2] - observations[1]) / (2 * s)
        state = {'schema': 'one-crown-positive-corridor-v1', 'case': old['case'], 'geometry_sha256': old['geometry_sha256'], 'A': old['A'], 'w_N': old['w_N'], 'predicted_fdi': old['predicted_fdi'], 'action_fdi': old['predicted_fdi'][i], 'baseline_force_N': b.tolist(), 'height_column_N_per_mm': jj.tolist(), 'probe_height_mm': s, 'raw_channel_error_N': eps, 'linear_full_contact_closure': True, 'closure_provenance': 'Declared passive linear simulation; not externally validated tissue', 'raw_probe_force_N': [x.tolist() for x in observations], 'provenance_kind': 'SIMULATED_THREE_HEIGHT_ACQUISITION'}
        write('exports/' + old['case'] + '_corridor_state.json', state)
        states.append((old, state))
        fullcount += 1 + 2 * (n - 3)
        cal = predict(state, 0)
        calibrations.append({'case': old['case'], 'valid': cal['status'] == 'CONDITIONAL_CORRIDOR_BOUND', 'full_energy_vectors': 1 + 2 * (n - 3), 'targeted_vectors': 3})
        for h in m['query_heights_mm']:
            ans = predict(state, h)
            hh = np.zeros(n)
            hh[i] = h
            rows.append({'case': old['case'], 'height_mm': h, 'answer': ans, 'h_vector': hh.tolist()})
        inv = inverse(state, 0.6 * f0[i], 1)
        inverses.append({'case': old['case'], 'answer': inv})
    frozen = {'frozen_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha, 'kind': 'Synthetic predictions before independent QP validation; physical NOT_RUN', 'predictions': rows, 'inverse': inverses}
    pf = R / 'FROZEN_PREDICTIONS_R4.json'
    if not pf.exists():
        write(pf.name, frozen)
        (R / (pf.name + '.sha256')).write_text(hashlib.sha256(pf.read_bytes()).hexdigest() + '  ' + pf.name + '\n')
    else:
        write('raw/R4_REPLAY_PREDICTIONS.json', frozen)
    lookup = {o['case']: (o, s) for (o, s) in states}
    records = []
    for row in rows:
        (old, ss) = lookup[row['case']]
        ans = row['answer']
        if ans.get('force_interval_N') is None:
            records.append({'case': row['case'], 'covered': False, 'reason': ans['reason']})
            continue
        fc = conventional(old, np.array(row['h_vector']))
        bound = np.array(ans['force_interval_N'])
        bad = fc.copy()
        bad[0] += 3
        records.append({'case': row['case'], 'height_mm': row['height_mm'], 'true_force_N': fc.tolist(), 'covered': bool(np.all(fc >= bound[:, 0] - 1e-06) and np.all(fc <= bound[:, 1] + 1e-06)), 'error_l2_N': float(np.linalg.norm(fc - ans['force_N'])), 'radius_N': ans['joint_l2_error_radius_N'], 'fault_detected': bool(not (np.all(bad >= bound[:, 0] - 1e-06) and np.all(bad <= bound[:, 1] + 1e-06)))})
    for row in inverses:
        (old, ss) = lookup[row['case']]
        a = row['answer']
        i = old['predicted_fdi'].index(ss['action_fdi'])
        if 'at_selected_height' in a:
            hh = np.zeros(len(old['baseline_force_N']))
            hh[i] = a['nominal_height_change_mm']
            fc = conventional(old, hh)
            row['true_target_band_pass'] = bool(abs(fc[i] - a['target_N']) <= 1)
            row['true_force_N'] = float(fc[i])
    Hs = [[[F(128), F(-64)], [F(-64), F(128)]], [[F(128), F(-64)], [F(-64), F(256)]]]
    base = [F(16), F(2), F(14)]
    col = [F(128), F(-64), F(-64)]
    probe = F(1, 64)
    hq = F(3, 64)
    summaries = [[float(v) for v in base + col + [base[k] + sign * probe * col[k] for sign in [-1, 1] for k in range(3)]] for H in Hs]
    future = []
    for H in Hs:
        det = H[0][0] * H[1][1] - H[0][1] ** 2
        z2 = F(-2)
        z1 = (hq * det + H[0][1] * z2) / H[1][1]
        future.append([float(base[0] + z1), float(base[1] + z2), float(base[2] - z1 - z2)])
    suff = {'summary_identity_error': float(np.max(np.abs(np.array(summaries[0]) - summaries[1]))), 'summary_bit_identical': np.array(summaries[0]).tobytes() == np.array(summaries[1]).tobytes(), 'summary_values': summaries, 'downstream_height_mm': float(hq), 'downstream_force_N': future, 'tooth1_difference_N': abs(future[0][0] - future[1][0]), 'minimum_extension': 'For this3-support contact branch, one extra transverse energy scalar H22, or a direct force vector at the new height. General branch changes need full H or another measured corridor.', 'status': 'ONE_COLUMN_SUMMARY_INSUFFICIENT_AFTER_CONTACT_CHANGE', 'referent': 'Closed-form active branch z2=-2; G=H^-1; z1=(h*det(H)+H12*z2)/H22'}
    ss = states[0][1]
    unknown = predict(ss, 0.046875)['status'] == 'UNKNOWN'
    bad = dict(ss, baseline_force_N=[0.001] * len(ss['baseline_force_N']))
    bad = dict(ss, height_column_N_per_mm=(np.array(ss['height_column_N_per_mm']) * 100).tolist())
    guard = predict(bad, 0)['status'] == 'UNKNOWN'
    out = {'round': 'R4', 'claim_type': 'capability', 'prereg_sha256': sha, 'external_referent': cfg['external_referent'], 'calibrations': calibrations, 'valid_calibrations': sum((x['valid'] for x in calibrations)), 'forward_queries': len(rows), 'contained': sum((x['covered'] for x in records)), 'max_radius_N': max((x.get('radius_N', 0) for x in records)), 'inverse_queries': len(inverses), 'inverse_guaranteed': sum((x['answer']['status'] == 'CONDITIONAL_TARGET_BAND' for x in inverses)), 'inverse_true_band_pass': sum((x.get('true_target_band_pass', False) for x in inverses)), 'inverse': inverses, 'summary_sufficiency': suff, 'gates': {'coverage': all((x['covered'] for x in records)), 'fault_injection': all((x.get('fault_detected', False) for x in records)), 'outside_corridor_unknown': unknown, 'possible_contact_change_unknown': guard}, 'cost': {'wall_s': time.perf_counter() - t, 'targeted_force_vectors': 36, 'full_H_force_vectors': fullcount, 'physical_acquisition': 'NOT_RUN', 'fit': '3 simulated full force vectors per layout; no full H used by query consumer', 'forward_query_calls': len(rows), 'inverse_queries': len(inverses), 'validation_QP_calls': len(rows) + len(inverses)}, 'resolution': 'PER_TOOTH; joint L2 radius PER_ARCH', 'physical_validity': 'UNKNOWN; three settings do not establish nonlinear tissue law'}
    write('raw/R4_CONTROLS.json', records)
    write('raw/R4_SUFFICIENCY.json', suff)
    write('rounds/R4/results.json', out)
    write('CURRENT_WORK_STATE.json', {'milestone': 'R4_TARGETED_ACQUISITION_COMPLETE', 'latest_gate': out['gates'], 'next_operation': 'Package standalone demo and lab protocol; independent review pending'})
    print(json.dumps({k: out[k] for k in ['valid_calibrations', 'forward_queries', 'contained', 'max_radius_N', 'inverse_guaranteed', 'inverse_true_band_pass', 'gates', 'cost']}, indent=2))
if __name__ == '__main__':
    main()
