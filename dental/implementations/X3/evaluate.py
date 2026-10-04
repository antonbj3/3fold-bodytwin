"""Only this stage reads hidden original tooth surfaces, after prediction freeze."""
from crownbench import *
import argparse, resource, collections
METRICS = ('surface_rms_mm', 'surface_p95_mm', 'occlusal_p95_mm', 'proximal_gap_error_mm', 'proximal_patch_p95_mm', 'contact_location_mm')
TYPES = {2: 'incisor', 3: 'canine', 5: 'premolar', 6: 'molar'}

def metric(pred, ref, neighbors, n, key, count=2048):
    (pv, pf) = pred
    (rv, rf) = ref
    pp = sampled(pv, pf, count, seed(key + '/surface'))
    rp = sampled(rv, rf, count, seed(key + '/surface'))
    pd = distances(pp, rv, rf)
    rd = distances(rp, pv, pf)
    po = pp @ n >= np.quantile(pp @ n, 0.65)
    ro = rp @ n >= np.quantile(rp @ n, 0.65)
    ans = {'surface_rms_mm': float(np.sqrt((np.mean(pd ** 2) + np.mean(rd ** 2)) / 2)), 'surface_p95_mm': float(max(np.quantile(pd, 0.95), np.quantile(rd, 0.95))), 'occlusal_p95_mm': float(max(np.quantile(pd[po], 0.95), np.quantile(rd[ro], 0.95)))}
    gaps = []
    patches = []
    locations = []
    raw = []
    for (idx, (nv, nf)) in enumerate(neighbors):
        ns = sampled(nv, nf, 1024, seed(key + '/neighbor/' + str(idx)))
        dt = distances(ns, rv, rf)
        dp = distances(ns, pv, pf)
        patch = dt <= np.quantile(dt, 0.05)
        gaps.append(float(abs(np.quantile(dp, 0.05) - np.quantile(dt, 0.05))))
        patches.append(float(np.quantile(abs(dp[patch] - dt[patch]), 0.95)))
        pc = distances(pp, nv, nf)
        rc = distances(rp, nv, nf)
        cp = pp[pc <= np.quantile(pc, 0.02)].mean(0)
        cr = rp[rc <= np.quantile(rc, 0.02)].mean(0)
        locations.append(float(np.linalg.norm(cp - cr)))
        raw.append({'neighbor': idx, 'reference_gap_q05_mm': float(np.quantile(dt, 0.05)), 'predicted_gap_q05_mm': float(np.quantile(dp, 0.05)), 'reference_patch_centroid': cr.tolist(), 'predicted_patch_centroid': cp.tolist(), 'patch_samples': int(patch.sum())})
    ans.update(proximal_gap_error_mm=float(np.mean(gaps)), proximal_patch_p95_mm=float(max(patches)), contact_location_mm=float(np.mean(locations)), proximal_details=raw, antagonist_contact=None, antagonist_status='UNKNOWN_NO_VERIFIED_BITE_TRANSFORM')
    return ans

def balanced(rows, method, metric_name):
    vals = collections.defaultdict(list)
    for r in rows:
        if method in r['metrics']:
            vals[r['case']].append(r['metrics'][method][metric_name])
    return {c: float(np.mean(v)) for (c, v) in vals.items()}

def gate_summary(rows, method, prereg, scheduled, roundtrip):
    base = 'mirror'
    b = balanced(rows, base, 'occlusal_p95_mm')
    a = balanced(rows, method, 'occlusal_p95_mm')
    cases = sorted(set(a) & set(b))
    if not cases:
        return {'decision': 'UNKNOWN_NO_OUTPUTS'}
    av = np.array([a[c] for c in cases])
    bv = np.array([b[c] for c in cases])
    diff = av - bv
    rng = np.random.default_rng(6103)
    ii = rng.integers(0, len(cases), (1000, len(cases)))
    ci = np.percentile(diff[ii].mean(1), [2.5, 97.5])
    proxA = balanced(rows, method, 'proximal_patch_p95_mm')
    proxB = balanced(rows, base, 'proximal_patch_p95_mm')
    prox_gain = np.mean([proxA[c] - proxB[c] for c in sorted(set(proxA) & set(proxB))])
    complete = sum((method in r['metrics'] for r in rows))
    types = {TYPES[r['fdi'] % 10] for r in rows if method in r['metrics']}
    minimum = 8 if prereg['round'].startswith('R1') else 6
    data = len(cases) >= minimum and len(types) >= 3 and (complete / scheduled >= 0.8)
    primary = av.mean() <= 0.9 * bv.mean() and ci[1] < 0
    absolute = av.mean() <= 0.5
    prox = prox_gain <= 0.05
    num = roundtrip['passes']
    control = 'neighbor_mesh' if method == 'neighbor_sdf' else 'ring_mesh'
    ca = balanced(rows, control, 'occlusal_p95_mm')
    cc = sorted(set(ca) & set(a))
    cd = float(np.mean([a[c] - ca[c] for c in cc])) if cc else None
    return {'decision': 'POSITIVE_ANATOMICAL_RECONSTRUCTION' if all([data, primary, absolute, prox, num]) else 'NEGATIVE', 'gates': {'data': 'PASS' if data else 'FAIL', 'primary_relative_and_bootstrap': 'PASS' if primary else 'FAIL', 'absolute_0.50mm': 'PASS' if absolute else 'FAIL', 'proximal_noninferiority_0.05mm': 'PASS' if prox else 'FAIL', 'SDF_train_roundtrip': 'PASS' if num else 'FAIL'}, 'method': method, 'scheduled': scheduled, 'completed': complete, 'eligible_fraction': complete / scheduled, 'test_patients': len(cases), 'mean_occlusal_p95_mm': float(av.mean()), 'mirror_mean_occlusal_p95_mm': float(bv.mean()), 'ratio_to_mirror': float(av.mean() / bv.mean()), 'paired_patient_difference_mm': float(diff.mean()), 'paired_patient_bootstrap95ci_mm': ci.tolist(), 'paired_proximal_patch_difference_mm': float(prox_gain), 'strongest_equal_information_control': control, 'candidate_minus_control_occlusal_mm': cd, 'algorithmic_comparison': 'TIE_WITHIN_0.05mm' if cd is not None and abs(cd) <= 0.05 else 'NO_NOVELTY_CLAIM', 'clinical_acceptability': 'UNKNOWN', 'antagonist_contact': 'UNKNOWN_NO_VERIFIED_BITE_TRANSFORM'}

def controls(ref, neighbors, n, F, key):
    (rv, rf) = ref
    identity = metric(ref, ref, neighbors, n, key + '/identity')
    shifted = metric((rv + 3 * F[:, 1], rf), ref, neighbors, n, key + '/translation')
    cusp = rv.copy()
    z = rv @ n
    top = z >= np.quantile(z, 0.65)
    cusp[top] += 0.8 * n
    raised = metric((cusp, rf), ref, neighbors, n, key + '/raised_cusp')
    (nv, nf) = neighbors[0]
    dv = distances(rv, nv, nf)
    ind = dv <= np.quantile(dv, 0.15)
    direction = unit(rv[ind].mean(0) - nv.mean(0))
    cv = rv.copy()
    cv[ind] += direction
    proximal = metric((cv, rf), ref, neighbors, n, key + '/proximal_patch')
    checks = {'identity_surface': identity['surface_p95_mm'] <= 1e-07, 'translation3mm_fails_surface0.50mm': shifted['surface_p95_mm'] > 0.5, 'raised0.8mm_fails_occlusal0.50mm': raised['occlusal_p95_mm'] > 0.5, 'local_proximal1mm_fails_patch0.25mm': proximal['proximal_patch_p95_mm'] > 0.25}
    injected_values = {'coverage': {'value': 0.0, 'gate_pass': 0.0 >= 0.8}, 'primary_ratio': {'value': 2.0, 'gate_pass': 2.0 <= 0.9}, 'bootstrap_upper': {'value': 1.0, 'gate_pass': 1.0 < 0}, 'absolute': {'value': 2.0, 'gate_pass': 2.0 <= 0.5}, 'proximal': {'value': 1.0, 'gate_pass': 1.0 <= 0.05}, 'numerical': {'value': 1.0, 'gate_pass': 1.0 <= 0.3}, 'same_information_tie': {'value': 1.0, 'gate_pass': abs(1.0) <= 0.05}}
    checks['all_injected_threshold_values_rejected'] = all((not v['gate_pass'] for v in injected_values.values()))
    return {'reference_key': key, 'kind': 'our_own_fixture', 'external_reference': False, 'checks': checks, 'all_pass': all(checks.values()), 'identity': identity, 'shift3mm': shifted, 'raised_topband0.8mm': raised, 'local_proximal1mm': proximal, 'injected_threshold_values': injected_values}

def evaluate(round_name):
    start = time.perf_counter()
    group = 'test_R1' if round_name == 'R1' else 'reserve_R2'
    fr = json.load(open(P / f'FROZEN_PREDICTIONS_{round_name}.json'))
    prereg = json.load(open(P / f'PREREG_{round_name}.json'))
    assert fr['all_predictions_complete'] and fr['hidden_evaluation_started'] == False
    assert sha(fr['predictions_index']) == fr['index_sha256'], 'PREDICTION_INDEX_HASH_DRIFT'
    for r in fr['prediction_files']:
        assert sha(r['path']) == r['sha256'], 'PREDICTION_HASH_DRIFT'
    hidden = {r['key']: r for r in json.load(open(DATA / 'hidden' / f'{group}_INDEX.json'))}
    pred = json.load(open(fr['predictions_index']))
    rows = []
    probe = None
    state(round_name + '_EVALUATING', 'FROZEN_HASHES_VERIFIED', 'measure withheld real surfaces and local neighbor patches')
    raw = P / f'RAW_METRICS_{round_name}.jsonl'
    with open(raw, 'w', buffering=1) as out:
        for rec in pred:
            row = {'key': rec['key'], 'case': rec['case'], 'jaw': rec['jaw'], 'fdi': rec['fdi'], 'tooth_type': TYPES[rec['fdi'] % 10], 'eligible': rec['eligible'], 'metrics': {}}
            if not rec['methods']:
                row['failure'] = rec.get('generation_failure', rec.get('reason'))
                rows.append(row)
                out.write(json.dumps(row) + '\n')
                continue
            hr = hidden[rec['key']]
            assert sha(hr['reference']) == hr['reference_sha256'], 'EXTERNAL_REFERENCE_HASH_DRIFT'
            ref = read_surface(hr['reference'])
            with np.load(hr['context']) as d:
                cv = d['vertices']
                cf = d['faces']
                cl = d['labels']
            k = rec['fdi']
            assert k not in np.unique(cl)
            neighbors = [select(cv, cf, cl, k - 1), select(cv, cf, cl, k + 1)]
            F = np.array(rec['frame'])
            n = F[:, 2]
            for (name, det) in rec['methods'].items():
                ts = time.perf_counter()
                try:
                    row['metrics'][name] = dict(metric(read_surface(det['path']), ref, neighbors, n, rec['key']), evaluation_seconds=time.perf_counter() - ts)
                except Exception as e:
                    row.setdefault('evaluation_failures', {})[name] = type(e).__name__ + ': ' + str(e)
            if probe is None and 'neighbor_sdf' in row['metrics']:
                probe = (ref, neighbors, n, F, rec['key'], read_surface(rec['methods']['neighbor_sdf']['path']), row['metrics']['neighbor_sdf'])
            rows.append(row)
            out.write(json.dumps(row, allow_nan=False) + '\n')
            if len(rows) % 4 == 0:
                print('EVALUATE', round_name, len(rows), 'of', len(pred), 'seconds', round(time.perf_counter() - start, 1), flush=True)
    methods = sorted({m for r in rows for m in r['metrics']})
    tables = []
    for jaw in ['all', 'upper', 'lower']:
        for tooth_type in ['all', 'incisor', 'canine', 'premolar', 'molar']:
            sub = [r for r in rows if (jaw == 'all' or r['jaw'] == jaw) and (tooth_type == 'all' or r['tooth_type'] == tooth_type)]
            for m in methods:
                cases = balanced(sub, m, 'occlusal_p95_mm')
                table = {'jaw': jaw, 'tooth_type': tooth_type, 'method': m, 'teeth': sum((m in r['metrics'] for r in sub)), 'patients': len(cases)}
                if cases:
                    table.update({q: float(np.mean(list(balanced(sub, m, q).values()))) for q in METRICS})
                    tables.append(table)
    candidate = 'neighbor_sdf' if round_name == 'R1' else 'ring_sdf'
    train = json.load(open(P / 'TRAIN_FIT.json'))
    decision = gate_summary(rows, candidate, prereg, len(pred), train['roundtrip'])
    control = None
    sampling = None
    if probe:
        (ref, neighbors, n, F, key, predprobe, original) = probe
        control = controls(ref, neighbors, n, F, key)
        sampling2 = metric(predprobe, ref, neighbors, n, key, count=8192)
        sampling = {'reference_key': key, 'counts': [2048, 8192], 'metric_differences_mm': {q: sampling2[q] - original[q] for q in METRICS}, 'meaning': 'conditional sensitivity for one frozen output, not a general numerical guarantee'}
        if not control['all_pass']:
            decision['decision'] = 'NEGATIVE_METRIC_SENSITIVITY_FAILURE'
    put(P / f'METRIC_CONTROLS_{round_name}.json', control)
    put(P / f'SAMPLING_CHECK_{round_name}.json', sampling)
    res = {'schema': 'X3-crown-results-v1', 'round': round_name, 'evaluated_utc': now(), 'external_referent': prereg['external_referent'], 'prereg_sha256': sha(P / f'PREREG_{round_name}.json'), 'split_sha256': sha(P / 'SPLIT.json'), 'frozen_predictions_sha256': sha(P / f'FROZEN_PREDICTIONS_{round_name}.json'), 'raw_metrics_sha256': sha(raw), 'decision': decision, 'metric_controls_all_pass': control['all_pass'] if control else False, 'tables': tables, 'scheduled': len(pred), 'eligible': sum((r['eligible'] for r in pred)), 'generation_failures': [{'key': r['key'], 'reason': r['failure']} for r in rows if r.get('failure')], 'evaluation_seconds': time.perf_counter() - start, 'peak_rss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'preparation': json.load(open(P / f'PREPARATION_{group}.json')), 'generation': {k: fr[k] for k in ['fit_and_generation_seconds', 'peak_rss_kb', 'failure_count', 'data_bytes']}, 'train_seconds': train['seconds'], 'source_measurement_uncertainty': 'UNKNOWN scanner/annotation uncertainty; observed errors conditional on supplied geometry', 'scope': 'Natural exposed tooth completion, provisional mm. No clinical crowns, force contact, registered antagonist, prep margin or manufacturing validation. No direct CrownGen/MADCrowner performance comparison: those models not executed.', 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    put(P / f'RESULTS_{round_name}.json', res)
    state(round_name + '_COMPLETE', decision, 'preserve HANDOFF then change operation toward same parent capability')
    print(json.dumps({'round': round_name, 'decision': decision, 'controls': control['checks'] if control else None, 'seconds': res['evaluation_seconds']}), flush=True)
    return res
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--round', default='R1', choices=['R1', 'R2'])
    a = ap.parse_args()
    evaluate(a.round)
