from dental_release.paths import expand as _release_expand
import ast, copy, hashlib, importlib, json, sys, time, resource
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, LinearConstraint
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
from height_force import predict, inverse, make_edit, project_reaction

def write(p, x):
    if Path(p).name == 'CURRENT_WORK_STATE.json':
        x = dict(x, lane='X82-height-to-force', updated_at_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat())
        with (ROOT / 'raw/MILESTONES.jsonl').open('a') as stream:
            stream.write(json.dumps(x) + '\n')
    dest = ROOT / p
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.name == 'CURRENT_WORK_STATE.json':
        temp = dest.with_suffix('.tmp')
        temp.write_text(json.dumps(x, indent=2, allow_nan=False))
        temp.replace(dest)
    else:
        dest.write_text(json.dumps(x, indent=2, allow_nan=False))

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def check(tag):
    p = ROOT / f'PREREG_{tag}.json'
    assert sha(p) == (ROOT / (p.name + '.sha256')).read_text().split()[0]
    return sha(p)

def conventional(s, h):
    N = np.array(s['basis_N'])
    G = np.linalg.inv(s['H_N_per_mm'])
    f0 = np.array(s['baseline_force_N'])
    r = N.T @ h
    sol = minimize(lambda z: 0.5 * z @ G @ z - r @ z, np.zeros(N.shape[1]), jac=lambda z: G @ z - r, constraints=[LinearConstraint(N, -f0, np.full(len(f0), np.inf))], method='SLSQP', options={'ftol': 1e-13, 'maxiter': 1000})
    f = f0 + N @ sol.x
    if not sol.success or f.min() < -1e-06:
        raise ArithmeticError(sol.message)
    return f

def first_round():
    t = time.perf_counter()
    check('R0')
    s = json.loads((ROOT / 'inputs/F5367_R4_state.json').read_text())
    e = make_edit(s, s['predicted_fdi'].index(16), -0.04)
    from active_contact import solve
    source_a = (ROOT / 'inputs/predict_edit.py').read_text()
    node_a = next((x for x in ast.parse(source_a).body if isinstance(x, ast.FunctionDef) and x.name == 'predict'))
    ns_a = {'np': np, 'solve': solve}
    exec(ast.get_source_segment(source_a, node_a), ns_a)
    from types import SimpleNamespace
    legacy = SimpleNamespace(predict=ns_a['predict'])
    bad = copy.deepcopy(s)
    bad['baseline_error_N'] = -1
    old = legacy.predict(bad, e)
    write('raw/R0_LEGACY_NEGATIVE_UNCERTAINTY.json', old)
    source = (ROOT / 'inputs/run_r6.py').read_text()
    node = next((x for x in ast.parse(source).body if isinstance(x, ast.FunctionDef) and x.name == 'project'))
    ns = {}
    exec(ast.get_source_segment(source, node), ns)
    flags = {'matched_velocity': True, 'matched_history': True, 'position_closure': True}
    oldb = ns['project'](1, 2, 0.005, 0.01, 0, -0.1, **flags)
    write('raw/R0_LEGACY_REVERSED_INTERVAL.json', oldb)
    injections = []
    for name in ['spectral_error_N_per_mm', 'baseline_error_N']:
        for (value, label) in [(-1, 'negative'), (float('nan'), 'NaN'), (float('inf'), 'infinity'), (True, 'boolean')]:
            b = copy.deepcopy(s)
            b[name] = value
            try:
                predict(b, e)
                reject = False
            except ValueError:
                reject = True
            injections.append({'consumer': 'balanced energy', 'field': name, 'value': label, 'rejected': reject})
    for ix in ['de', 'fe']:
        for (value, label) in [(-0.1, 'negative'), (float('nan'), 'NaN'), (float('inf'), 'infinity'), (True, 'boolean')]:
            args = [1, 2, 0.005, 0.01, 0, 0.01]
            args[4 if ix == 'de' else 5] = value
            try:
                project_reaction(*args, **flags)
                reject = False
            except ValueError:
                reject = True
            injections.append({'consumer': 'reaction', 'field': ix, 'value': label, 'rejected': reject})
    a = predict(s, e)
    aold = legacy.predict(s, e)
    valid_diff = np.max(np.abs(np.array(a['force_N']) - aold['force_N']))
    out = {'round': 'R0', 'claim_type': 'capability', 'prereg_sha256': check('R0'), 'legacy_failure_A': {'status': old['status'], 'radius_N': old['joint_l2_error_radius_N'], 'reversed_intervals': sum((l > u for (l, u) in old['force_interval_N']))}, 'legacy_failure_B': oldb, 'injections': injections, 'rejected': sum((r['rejected'] for r in injections)), 'total': len(injections), 'valid_force_difference_N': float(valid_diff), 'gate': bool(all((r['rejected'] for r in injections)) and valid_diff <= 1e-07), 'cost_wall_s': time.perf_counter() - t}
    write('rounds/R0/results.json', out)
    write('CURRENT_WORK_STATE.json', {'milestone': 'R0_REPAIR_PASS', 'latest_gate': out['gate'], 'next_operation': 'Forward/inverse fixed-wrench questions plus exact identity test'})
    return out

def second_round():
    t = time.perf_counter()
    cpu = time.process_time()
    cfg = json.loads((ROOT / 'PREREG_R1.json').read_text())
    check('R1')
    predictions = []
    inverse_rows = []
    states = []
    for p in sorted((ROOT / 'inputs').glob('*_R4_state.json')):
        s = json.loads(p.read_text())
        states.append(s)
        for (i, tooth) in enumerate(s['predicted_fdi']):
            for h in cfg['metrics']['forward_points_per_tooth']:
                e = make_edit(s, i, h)
                a = predict(s, e)
                predictions.append({'case': s['case'], 'fdi': tooth, 'height_mm': h, 'edit': e, 'answer': a})
        i = s['predicted_fdi'].index(16) if 16 in s['predicted_fdi'] else 0
        target = 0.6 * s['baseline_force_N'][i]
        ans = inverse(s, s['predicted_fdi'][i], target)
        inverse_rows.append({'case': s['case'], 'answer': ans})
    frozen = {'at_utc': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), 'kind': 'Synthetic predictions frozen before numerical controls; not before external publication or a physical measurement', 'prereg_sha256': check('R1'), 'state_hashes': {p.name: sha(p) for p in sorted((ROOT / 'inputs').glob('*_R4_state.json'))}, 'predictions': predictions, 'inverse': inverse_rows}
    p = ROOT / 'FROZEN_PREDICTIONS.json'
    if not p.exists():
        write('FROZEN_PREDICTIONS.json', frozen)
        (ROOT / 'FROZEN_PREDICTIONS.json.sha256').write_text(sha(p) + '  FROZEN_PREDICTIONS.json\n')
    else:
        write('raw/R1_REPLAY_PREDICTIONS.json', frozen)
    controls = []
    lookup = {s['case']: s for s in states}
    for row in predictions:
        s = lookup[row['case']]
        fc = conventional(s, np.array(row['edit']['height_change_mm']))
        ans = row['answer']
        intervals = np.array(ans['force_interval_N'])
        err = float(np.max(np.abs(fc - ans['force_N'])))
        covered = bool(np.all(fc >= intervals[:, 0] - 1e-06) and np.all(fc <= intervals[:, 1] + 1e-06))
        injected = fc.copy()
        injected[0] += 3
        detects = bool(not (np.all(injected >= intervals[:, 0] - 1e-06) and np.all(injected <= intervals[:, 1] + 1e-06)))
        controls.append({'case': s['case'], 'fdi': row['fdi'], 'height_mm': row['height_mm'], 'control_force_N': fc.tolist(), 'difference_N': err, 'covered': covered, 'injected_3N_detected': detects, 'contact_released': sum((x < 1e-07 for x in ans['force_N']))})
    for row in inverse_rows:
        a = row['answer']
        i = lookup[row['case']]['predicted_fdi'].index(a.get('fdi', 16))
        if 'at_selected_height' in a:
            fc = conventional(lookup[row['case']], np.array(a['at_selected_height']['height_change_mm']))
            row['independent_force_N'] = float(fc[i])
            row['independent_target_band_pass'] = bool(a['target_band_N'][0] <= fc[i] <= a['target_band_N'][1])
    from fractions import Fraction as F
    C1 = [[F(1, 64), F(0)], [F(0), F(1, 64)]]
    C2 = [[F(3, 128), F(0)], [F(0), F(1, 128)]]
    f0 = [F(16), F(16)]
    T0 = F(32)
    T1 = F(40)
    G = [C[0][0] + C[1][1] - C[0][1] - C[1][0] for C in [C1, C2]]
    gs = [[-sum((C[i][j] * f0[j] for j in range(2))) for i in range(2)] for C in [C1, C2]]
    forces = [(C[1][1] * T1 + g[1] - g[0]) / (C[0][0] + C[1][1]) for (C, g) in zip([C1, C2], gs)]
    sufficient = {'resolution': 'PER_TOOTH', 'summary': {'projected_compliance_mm_per_N': [float(v) for v in G], 'reference_force_N': [16, 16], 'reference_total_N': 32}, 'identity_error': float(abs(G[0] - G[1])), 'bit_identical': np.array([float(G[0])]).tobytes() == np.array([float(G[1])]).tobytes(), 'new_total_force_N': 40, 'new_tooth_1_force_N': [float(x) for x in forces], 'downstream_difference_N': float(abs(forces[0] - forces[1])), 'minimum_extension': 'One per-tooth reaction response to total-load change at matched height/history, or enforce identical wrench.', 'status': 'FIXED_WRENCH_SUMMARY_INSUFFICIENT_FOR_NEW_WRENCH', 'external_referent': {'kind': 'closed_form', 'locator': 'code/run_r0_r1.py dyadic two-support equal-displacement formula', 'compared_quantity': 'Exact equilibrium reactions in declared linear support model', 'refutes_us': True}}
    write('raw/SUFFICIENCY_NEW_WRENCH.json', sufficient)
    one = predictions[0]
    newedit = copy.deepcopy(one['edit'])
    newedit['w_N'] = (np.array(lookup[one['case']]['w_N']) * 1.05).tolist()
    guard = predict(lookup[one['case']], newedit)
    out = {'round': 'R1', 'claim_type': 'capability', 'prereg_sha256': check('R1'), 'forward_queries': len(predictions), 'contained': sum((x['covered'] for x in controls)), 'force_radius_le1N': sum((x['answer']['joint_l2_error_radius_N'] <= 1 for x in predictions)), 'max_radius_N': max((x['answer']['joint_l2_error_radius_N'] for x in predictions)), 'max_control_difference_N': max((x['difference_N'] for x in controls)), 'contact_release_queries': sum((x['contact_released'] > 0 for x in controls)), 'inverse_queries': len(inverse_rows), 'inverse_guaranteed': sum((x['answer']['status'] == 'CONDITIONAL_TARGET_BAND' for x in inverse_rows)), 'inverse_independent_band_pass': sum((x.get('independent_target_band_pass', False) for x in inverse_rows)), 'fault_injection_detected': sum((x['injected_3N_detected'] for x in controls)), 'fault_injection_total': len(controls), 'gates': {'coverage': all((x['covered'] for x in controls)), 'same_information_control': max((x['difference_N'] for x in controls)) <= 0.001, 'fault_injection': all((x['injected_3N_detected'] for x in controls)), 'new_wrench_guard': guard['status'] == 'UNKNOWN'}, 'summary_sufficiency': sufficient, 'inverse': inverse_rows, 'controls_file': 'raw/R1_CONTROLS.json', 'cost': {'wall_s': time.perf_counter() - t, 'cpu_s': time.process_time() - cpu, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'forward_calls': len(predictions), 'inverse_forward_calls': sum((x['answer']['query_calls'] for x in inverse_rows)), 'QP_control_calls': len(controls) + len(inverse_rows), 'physical_acquisition': 'NOT_RUN; reused earlier prospective simulation states; prior 1+2(n-3) readings per layout remain charged prerequisites'}}
    write('raw/R1_CONTROLS.json', controls)
    write('rounds/R1/results.json', out)
    write('exports/F5367_inverse.json', next((x['answer'] for x in inverse_rows if x['case'] == _release_expand('@DENTAL_CASE_ID@'))))
    write('CURRENT_WORK_STATE.json', {'milestone': 'R1_COMPLETE', 'latest_gate': out['gates'], 'next_operation': 'External observable transfer + query-specific minimal measurement protocol'})
    return {k: out[k] for k in ['forward_queries', 'contained', 'max_radius_N', 'max_control_difference_N', 'inverse_guaranteed', 'gates']}
if __name__ == '__main__':
    print(json.dumps(first_round(), indent=2))
    print(json.dumps(second_round(), indent=2))
