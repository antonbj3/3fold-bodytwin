import csv, json, time, hashlib, datetime, resource, copy, zipfile
import numpy as np
from scipy.stats import norm, spearmanr
from geometry import ROOT, DATA, apply, load_cloud
from wrench import build, solve_load, gate, PARENT
start = time.monotonic()
rows = json.loads((ROOT / 'raw/R1_POSES.json').read_text())
ix = {(r['patient'], r['jaw'], r['source'], r['week'], r['fdi']): r for r in rows}
prereg = json.loads((ROOT / 'PREREG_R3.json').read_text())
out = []
pred = []
failures = []
regions = ['Buccogingival', 'Buccal', 'Incisal/occlusal', 'Palatal', 'Palatogingival']

def upper(f):
    return f if f // 10 in [1, 2] else (2 if f // 10 == 3 else 1) * 10 + f % 10
parent_hashes = [{'file': str(PARENT / n), 'sha256': hashlib.sha256((PARENT / n).read_bytes()).hexdigest()} for n in ['whole_arch.py', 'whole_arch_refined.py', 'contact_model.py']]
for p in ['3485', '6457']:
    zf = zipfile.ZipFile(DATA / f'{p}_Clinical_Trial.zip')
    for j in ['OK', 'UK']:
        ids = sorted({r['fdi'] for r in rows if r['patient'] == p and r['jaw'] == j})
        for week in range(1, 10):
            prev = {f: np.array(ix[p, j, 'Bottmedical', week - 1, f]['H']['all_crowns']) for f in ids}
            next_ = {f: np.array(ix[p, j, 'Bottmedical', week, f]['H']['all_crowns']) for f in ids}
            centers = {f: apply(prev[f], np.array(ix[p, j, 'Bottmedical', week - 1, f]['centroid_base_mm'])) for f in ids}
            mean = np.mean(list(centers.values()), axis=0)
            arch = {'frame': {'occlusal_unit': [0.0, 1.0, 0.0]}, 'teeth': {}}
            back = {upper(f): f for f in ids}
            for f in ids:
                cc = centers[f]
                e = cc - mean
                e[1] = 0
                e /= np.linalg.norm(e)
                cloud = load_cloud(zf, p, 'Bottmedical', week - 1, f)
                yy = cloud['points'][:, 1]
                height = float(np.quantile(yy, 0.95) - np.quantile(yy, 0.05))
                arch['teeth'][str(upper(f))] = {'centroid_mm': cc.tolist(), 'crown_height_mm': height, 'buccal_unit': e.tolist()}
            for sc in prereg['scenarios']:
                title = f"{p}_{j}_W{week}_{sc['name']}"
                records = [{'arm': 'TS', 'quantity': q, 'region': rg, 'median_mm': sc['thickness_mm'] if q == 'finished_thickness' else sc['gap_mm']} for q in ['finished_thickness', 'passive_gap'] for rg in regions]
                try:
                    s = build(arch, records)
                    b = []
                    for patch in s['patch']:
                        f = back[patch['fdi']]
                        D = next_[f] @ np.linalg.inv(prev[f])
                        point = patch['point']
                        activation = point - apply(D, point)
                        b.append(patch['normal'] @ activation - patch['gap'])
                    b = np.array(b)
                    a = solve_load(s, b, 'dual')
                    c = solve_load(s, b, 'primal')
                    g = gate(a, c)
                    inj = copy.deepcopy(a)
                    inj['force_balance_N'] += 0.01
                    inj_m = copy.deepcopy(a)
                    inj_m['moment_balance_Nmm'] += 0.01
                    inj_gap = copy.deepcopy(a)
                    inj_gap['minimum_gap_mm'] = -0.01
                    injected = {'bad_force_rejected': not gate(inj, c)['pass'], 'bad_moment_rejected': not gate(inj_m, c)['pass'], 'bad_gap_rejected': not gate(inj_gap, c)['pass']}
                    stage = {'key': title, 'patient': p, 'jaw': j, 'week': week, 'scenario': sc['name'], 'dual': a, 'primal': c, 'gate': g, 'injections': injected, 'moment_origin': 'planned previous-stage crown material centroid, canonical all-crown frame', 'physical_force_status': 'UNKNOWN_UNMEASURED_NATURALIGNER_PARAMETERS_AND_SEATED_GEOMETRY', 'actual_load_history': 'UNKNOWN', 'input_kind': 'planned consecutive stage shapes; not actual lagged seated shape', 'source_patient_match': True, 'time_scale': 'HANDOVER'}
                    out.append(stage)
                    if not g['pass']:
                        failures.append({'key': title, 'reason': 'numerical_gate', 'gate': g})
                    for (key, w) in a['wrenches'].items():
                        f = back[int(key)]
                        pred.append({'patient': p, 'jaw': j, 'fdi': f, 'week': week, 'scenario': sc['name'], 'F_N': w[:3], 'M_Nmm': w[3:], 'force_norm_N': float(np.linalg.norm(w[:3])), 'moment_norm_Nmm': float(np.linalg.norm(w[3:])), 'origin_mm': centers[f].tolist(), 'numerical_admissible': g['pass'], 'resolution': 'PER_TOOTH', 'physical_status': 'UNKNOWN', 'time_scale_to_movement': 'HANDOVER', 'moment_origin': 'planned previous-stage crown material centroid'})
                except Exception as e:
                    out.append({'key': title, 'exception': type(e).__name__ + ': ' + str(e), 'numerical_admissible': False})
                    failures.append({'key': title, 'reason': 'solver_exception', 'error': str(e)})
                print(title, 'last_failure', bool(failures and failures[-1]['key'] == title), 'elapsed_s', round(time.monotonic() - start, 1), flush=True)
for (n, d) in [('R3V2_STAGES.json', out), ('R3V2_WRENCHES.json', pred), ('R3V2_PARENT_HASHES.json', parent_hashes)]:
    (ROOT / 'raw' / n).write_text(json.dumps(d, indent=2) + '\n')
frozen = {'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': hashlib.sha256((ROOT / 'PREREG_R3.json').read_bytes()).hexdigest(), 'prediction_file': 'raw/R3V2_WRENCHES.json', 'sha256': hashlib.sha256((ROOT / 'raw/R3V2_WRENCHES.json').read_bytes()).hexdigest(), 'before_descriptive_comparison': True, 'before_original_measurement': False, 'blinding': 'Motion already measured in R1. This freeze checks no fit during comparison, not prospective validation.'}
freeze = ROOT / 'FROZEN_PREDICTIONS_R3V2.json'
if not freeze.exists():
    freeze.write_text(json.dumps(frozen, indent=2) + '\n')
else:
    old = json.loads(freeze.read_text())
    if old['sha256'] != frozen['sha256']:
        raise RuntimeError('Prediction drift against original freeze')
moves = json.loads((ROOT / 'raw/R1_WEEKLY_MEASUREMENTS.json').read_text())
mi = {(r['patient'], r['jaw'], r['week'], r['fdi']): r for r in moves}
joined = []
for row in pred:
    m = mi.get((row['patient'], row['jaw'], row['week'], row['fdi']))
    if m is not None:
        d = {**row, 'observed_weekly_mm': m['observed_weekly_mm'], 'planned_weekly_mm': m['planned_weekly_mm'], 'weekly_sensitivity_budget_mm': m['weekly_sensitivity_budget_mm'], 'absolute_motion_status': m['absolute_movement_status'], 'force_explains_motion': 'UNKNOWN'}
        joined.append(d)
diag = []
for p in ['3485', '6457']:
    for sc in prereg['scenarios']:
        r = [x for x in joined if x['patient'] == p and x['scenario'] == sc['name'] and x['numerical_admissible']]
        f = [x['force_norm_N'] for x in r]
        o = [x['observed_weekly_mm'] for x in r]
        pl = [x['planned_weekly_mm'] for x in r]
        diag.append({'patient': p, 'scenario': sc['name'], 'rows': len(r), 'independent_patients': 1, 'force_observed_spearman': float(spearmanr(f, o).statistic) if len(r) > 2 and np.ptp(f) > 0 else None, 'plan_observed_spearman': float(spearmanr(pl, o).statistic) if len(r) > 2 else None, 'p_value': None, 'interpretation': 'Descriptive scenario only; dependencies and reference error prevent force explanatory inference.'})
zsum = norm.ppf(0.975) + norm.ppf(0.8)
power = [{'partial_rho': rho, 'k_adjustment_covariates': 2, 'independent_patients_required': int(np.ceil(5 + zsum ** 2 / np.arctanh(rho) ** 2)), 'alpha': 0.05, 'power': 0.8, 'resolution': 'POPULATION', 'claim': 'conditional asymptotic patient-summary correlation power, not pilot-estimated n'} for rho in [0.3, 0.5, 0.7]]

def port_accept(r):
    return r['patient'] == '3485' and r['fdi'] == 11 and (r['week'] == 1) and (r.get('force_unit') == 'N') and (r.get('moment_unit') == 'Nmm') and (r.get('moment_origin_mm') is not None) and (r.get('response_interval') == [0, 1])
base = {'patient': '3485', 'fdi': 11, 'week': 1, 'force_unit': 'N', 'moment_unit': 'Nmm', 'moment_origin_mm': [0, 0, 0], 'response_interval': [0, 1]}
ports = {'matched_structural_row_accepts': port_accept(base), 'wrong_patient_rejected': not port_accept({**base, 'patient': 'different'}), 'wrong_force_unit_rejected': not port_accept({**base, 'force_unit': 'mN'}), 'wrong_stage_rejected': not port_accept({**base, 'week': 2}), 'missing_origin_rejected': not port_accept({**base, 'moment_origin_mm': None}), 'wrong_time_rejected': not port_accept({**base, 'response_interval': [1, 2]})}
result = {'round': 'R3', 'implementation_attempt': 'V2_FLOAT_FRAME_INPUT', 'claim_type': 'information_link', 'external_referent': prereg['external_referent'], 'stage_cases': len(out), 'numerically_admissible_stages': sum((x.get('gate', {}).get('pass', False) for x in out)), 'wrench_rows': len(pred), 'joined_rows': len(joined), 'physical_force_validation': 'UNKNOWN_NO_MATCHED_WRENCH', 'force_explanation': 'UNKNOWN_NOT_IDENTIFIABLE_FROM_TWO_PATIENTS', 'failures': failures, 'diagnostics': diag, 'sample_size_scenarios': power, 'sample_size_actual': 'UNKNOWN effect size/correlation/noise and load history not estimated', 'port_fault_checks': ports, 'elapsed_s': time.monotonic() - start, 'peak_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'parent_physical_gate': 'X27 R8 failed external neighbour forces; not promoted by reusing energy'}
for (n, d) in [('R3V2_RESULTS.json', result), ('R3V2_JOINED.json', joined)]:
    (ROOT / 'raw' / n).write_text(json.dumps(d, indent=2, allow_nan=False) + '\n')
print(json.dumps(result, indent=2))
state = json.loads((ROOT / 'CURRENT_WORK_STATE.json').read_text())
state.update(status='R3_COMPLETE', latest_gate='Conditional force interface executed; physical force explanation UNKNOWN', next_operation='R4 test available non-dental stable reference ROI, then package one-command demo', updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2) + '\n')
