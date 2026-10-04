import copy, datetime, hashlib, json, resource, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from height_force import predict, make_edit, inverse, at_load
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

def full_control(A, C, g, w, h, fstart):
    T = w[0]
    sol = minimize(lambda x: 0.5 * T * x @ C @ x + (g - h) @ x, np.array(fstart) / T, jac=lambda x: T * C @ x + g - h, constraints={'type': 'eq', 'fun': lambda x: A.T @ x - w / T, 'jac': lambda x: A.T}, bounds=[(0, None)] * len(fstart), method='SLSQP', options={'ftol': 1e-14, 'maxiter': 1000})
    if not sol.success:
        raise ArithmeticError(sol.message)
    f = sol.x * T
    if np.max(np.abs(A.T @ f - w)) > 1e-07:
        raise ArithmeticError('Full-control balance failure')
    return f

def main():
    t = time.perf_counter()
    cpu = time.process_time()
    p = R / 'PREREG_R3.json'
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    assert sha == (R / 'PREREG_R3.json.sha256').read_text().split()[0]
    cfg = json.loads(p.read_text())
    predictions = []
    worlds = {}
    inverse_rows = []
    states = []
    old_guards = []
    for (ix, path) in enumerate(sorted((R / 'inputs').glob('*_R4_state.json'))):
        s = json.loads(path.read_text())
        A = np.array(s['A'])
        N = np.array(s['basis_N'])
        H = np.array(s['H_N_per_mm'])
        f0 = np.array(s['baseline_force_N'])
        w0 = np.array(s['w_N'])
        n = len(f0)
        G = np.linalg.inv(H)
        P = np.eye(n) - N @ N.T
        v = N[:, 0]
        a = np.ones(n) / np.sqrt(n)
        C = N @ G @ N.T + P / 750.0 + 0.0005 * (np.outer(v, a) + np.outer(a, v))
        if np.linalg.eigvalsh(C)[0] <= 0:
            raise ArithmeticError('Non-passive hidden instrument')
        g = -C @ f0
        ciA = np.linalg.solve(C, A)
        K = A.T @ ciA
        qtrue = ciA @ np.linalg.solve(K, w0 / w0[0])
        fplus = f0 + qtrue * 10
        fminus = f0 - qtrue * 10
        if min(fplus.min(), fminus.min()) <= 0:
            raise ArithmeticError('Central load probes leave full contact')
        direction = N[:, 1]
        eta = cfg['metrics']['load_response_error_l2_N_per_N']
        observed_plus = fplus + direction * eta * 10 * 0.5
        observed_minus = fminus - direction * eta * 10 * 0.5
        qhat = (observed_plus - observed_minus) / 20
        ss = dict(s, load_response_N_per_N=qhat.tolist(), load_response_error_l2_N_per_N=eta, load_change_limit_N=15.0, provenance_kind='SIMULATED_LOAD_AND_HEIGHT_ACQUISITION', measurement_contract=s['measurement_contract'] + '; extra load-ray column from T0+/-10N, L2 error<=0.0005N/N, same rate/history')
        write('exports/' + s['case'] + '_load_state.json', ss)
        states.append(ss)
        worlds[s['case']] = {'A': A.tolist(), 'C': C.tolist(), 'g': g.tolist(), 'f0': f0.tolist(), 'w0': w0.tolist(), 'load_probe_plus_N': observed_plus.tolist(), 'load_probe_minus_N': observed_minus.tolist(), 'true_response_N_per_N': qtrue.tolist(), 'measurement_error_l2_N_per_N': float(np.linalg.norm(qhat - qtrue)), 'scope': 'our_own_fixture; hidden support is synthetic and used only by independent validation'}
        i = s['predicted_fdi'].index(16) if 16 in s['predicted_fdi'] else 0
        for dT in cfg['metrics']['total_load_changes_N']:
            T = w0[0] + dT
            for h in cfg['metrics']['height_edits_mm']:
                e = make_edit(s, i, h)
                e['load_N'] = float(T)
                ans = predict(ss, e)
                predictions.append({'case': s['case'], 'fdi': s['predicted_fdi'][i], 'load_change_N': dT, 'edit': e, 'answer': ans})
                if dT:
                    old_guards.append(predict(s, e))
            target = 0.6 * (f0 + qhat * dT)[i]
            inv = inverse(ss, s['predicted_fdi'][i], target, load_N=float(T))
            inverse_rows.append({'case': s['case'], 'load_N': float(T), 'answer': inv})
    frozen = {'frozen_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha, 'kind': 'Synthetic predictions frozen before full-compliance validation; no physical measurement', 'predictions': predictions, 'inverse': inverse_rows}
    pf = R / 'FROZEN_PREDICTIONS_R3.json'
    if not pf.exists():
        write(pf.name, frozen)
        (R / (pf.name + '.sha256')).write_text(hashlib.sha256(pf.read_bytes()).hexdigest() + '  ' + pf.name + '\n')
    else:
        write('raw/R3_REPLAY_PREDICTIONS.json', frozen)
    records = []
    lookup = {s['case']: s for s in states}
    for row in predictions:
        s = lookup[row['case']]
        world = worlds[row['case']]
        A = np.array(world['A'])
        C = np.array(world['C'])
        g = np.array(world['g'])
        w = np.array(world['w0']) * row['edit']['load_N'] / world['w0'][0]
        h = np.array(row['edit']['height_change_mm'])
        fc = full_control(A, C, g, w, h, np.array(world['f0']) * row['edit']['load_N'] / world['w0'][0])
        ans = row['answer']
        bounds = np.array(ans['force_interval_N'])
        equal = conventional(at_load(s, row['edit']['load_N']), h)
        records.append({'case': row['case'], 'load_change_N': row['load_change_N'], 'height_mm': float(h[np.argmax(np.abs(h))]), 'true_control_force_N': fc.tolist(), 'true_error_l2_N': float(np.linalg.norm(fc - ans['force_N'])), 'covered': bool(np.all(fc >= bounds[:, 0] - 1e-06) and np.all(fc <= bounds[:, 1] + 1e-06)), 'joint_radius_pass': bool(np.linalg.norm(fc - ans['force_N']) <= ans['joint_l2_error_radius_N'] + 1e-06), 'same_information_difference_N': float(np.max(np.abs(equal - ans['force_N'])))})
    for row in inverse_rows:
        ans = row['answer']
        world = worlds[row['case']]
        A = np.array(world['A'])
        C = np.array(world['C'])
        g = np.array(world['g'])
        T = row['load_N']
        w = np.array(world['w0']) * T / world['w0'][0]
        if 'at_selected_height' in ans:
            h = np.array(ans['at_selected_height']['height_change_mm'])
            fc = full_control(A, C, g, w, h, np.array(world['f0']) * T / world['w0'][0])
            i = lookup[row['case']]['predicted_fdi'].index(ans['fdi'])
            row['true_force_N'] = float(fc[i])
            row['true_target_band_pass'] = bool(ans['target_band_N'][0] <= fc[i] <= ans['target_band_N'][1])
    inj = copy.deepcopy(states[0])
    inj['load_response_N_per_N'][0] *= 1.2
    try:
        predict(inj, {**make_edit(inj, 0, 0.01), 'load_N': 315})
        fault = False
    except ValueError:
        fault = True
    out = {'round': 'R3', 'claim_type': 'capability', 'prereg_sha256': sha, 'external_referent': cfg['external_referent'], 'forward_queries': len(records), 'contained': sum((x['covered'] for x in records)), 'joint_radius_passes': sum((x['joint_radius_pass'] for x in records)), 'max_radius_N': max((x['answer']['joint_l2_error_radius_N'] for x in predictions)), 'max_true_error_l2_N': max((x['true_error_l2_N'] for x in records)), 'max_same_information_difference_N': max((x['same_information_difference_N'] for x in records)), 'inverse_queries': len(inverse_rows), 'inverse_guaranteed': sum((x['answer']['status'] == 'CONDITIONAL_TARGET_BAND' for x in inverse_rows)), 'inverse_true_band_pass': sum((x.get('true_target_band_pass', False) for x in inverse_rows)), 'old_consumer_unknown_count': sum((x['status'] == 'UNKNOWN' for x in old_guards)), 'gates': {'coverage': all((x['covered'] for x in records)), 'joint_radius': all((x['joint_radius_pass'] for x in records)), 'same_information_control': max((x['same_information_difference_N'] for x in records)) <= 0.001, 'gain_fault_detected': fault, 'unmeasured_load_guard': all((x['status'] == 'UNKNOWN' for x in old_guards))}, 'inverse': inverse_rows, 'cost': {'wall_s': time.perf_counter() - t, 'cpu_s': time.process_time() - cpu, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'new_simulated_force_vectors': 24, 'physical_force_vectors': 'NOT_RUN; +2 per layout relative to prior full height-response calibration', 'inverse_query_calls': sum((x['answer']['query_calls'] for x in inverse_rows)), 'controls': 2 * len(records) + len(inverse_rows)}, 'physical_validity': 'UNKNOWN; declared coupled passive linear instrument only'}
    write('raw/R3_SYNTHETIC_INSTRUMENTS.json', worlds)
    write('raw/R3_CONTROLS.json', records)
    write('rounds/R3/results.json', out)
    write('CURRENT_WORK_STATE.json', {'milestone': 'R3_LOAD_EXTENSION_COMPLETE', 'latest_gate': out['gates'], 'next_operation': 'Perturbed-response enclosure tests + lab protocol + consumer export'})
    print(json.dumps({k: out[k] for k in ['forward_queries', 'contained', 'joint_radius_passes', 'max_radius_N', 'max_true_error_l2_N', 'inverse_guaranteed', 'inverse_true_band_pass', 'gates']}, indent=2))
if __name__ == '__main__':
    main()
