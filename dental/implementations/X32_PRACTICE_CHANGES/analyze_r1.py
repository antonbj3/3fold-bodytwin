from dental_release.paths import expand as _release_expand
import csv, json, hashlib, time, collections
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent
PKG = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts/DEMO48_PACKAGE'))

def read(demo, name):
    p = ROOT / 'raw/frozen_decision_sources' / demo / name
    record = next((x for x in json.loads((ROOT / 'raw/DECISION_SOURCE_MANIFEST.json').read_text()) if x['demo'] == demo and Path(x['snapshot_path']).name == name))
    assert hashlib.sha256(p.read_bytes()).hexdigest() == record['sha256'], 'Frozen decision input drift'
    return (list(csv.DictReader(p.open())), str(p), hashlib.sha256(p.read_bytes()).hexdigest())

def bootstrap(rows, cluster_key='patient', replicates=10000):
    g = collections.defaultdict(list)
    for r in rows:
        g[r[cluster_key]].append(r)
    a = np.array([[len(v), sum((x['changed'] for x in v)), sum((x['more_restrictive'] for x in v)), sum((x['less_restrictive'] for x in v)), int(any((x['changed'] for x in v)))] for v in g.values()], dtype=float)
    rng = np.random.default_rng(3202)
    draws = rng.integers(0, len(a), (replicates, len(a)))
    s = a[draws].sum(axis=1)
    vals = 100 * s[:, 1:4] / s[:, 0, None]
    return {'unit_rate_ci95': np.quantile(vals, [0.025, 0.975], axis=0).T.tolist(), 'any_patient_change_ci95': np.quantile(100 * s[:, 4] / len(a), [0.025, 0.975]).tolist(), 'clusters': len(a), 'replicates': replicates, 'seed': 3202}

def summarize(rows, bootstrap_allowed=True):
    n = len(rows)
    p = collections.defaultdict(list)
    for r in rows:
        p[r['patient']].append(r)
    out = {'n': n, 'changed': sum((r['changed'] for r in rows)), 'more_restrictive': sum((r['more_restrictive'] for r in rows)), 'less_restrictive': sum((r['less_restrictive'] for r in rows)), 'clinical_direction': 'UNKNOWN_NO_MATCHED_CLINICAL_OUTCOME', 'unit': 'decisions per 100 eligible units', 'resolution_level': rows[0]['resolution_level']}
    out.update({k + '_per100': 100 * out[k] / n for k in ['changed', 'more_restrictive', 'less_restrictive']})
    out['patients_with_any_change'] = sum((any((r['changed'] for r in v)) for v in p.values()))
    out['patients'] = len(p)
    out['patients_any_change_per100'] = 100 * out['patients_with_any_change'] / len(p)
    if bootstrap_allowed:
        out['intervals'] = bootstrap(rows)
    out['cluster_type'] = rows[0].get('cluster_type', 'patient')
    if out['cluster_type'] != 'patient':
        for k in list(out):
            if 'patient' in k:
                out[k.replace('patients', 'source_families').replace('patient', 'source_family')] = out.pop(k)
        if 'intervals' in out:
            z = out['intervals']
            z['any_source_family_change_ci95'] = z.pop('any_patient_change_ci95')
    return out

def decision(patient, unit, q0, q1, threshold, kind, level, evidence):
    d0 = q0 >= threshold if kind == 'clearance' else q0 <= threshold
    d1 = q1 >= threshold if kind == 'clearance' else q1 <= threshold
    return dict(patient=patient, cluster_type='patient', object=unit, q0=q0, q1=q1, threshold=threshold, kind=kind, practice_accept=d0, demo_accept=d1, changed=d0 != d1, more_restrictive=d0 and (not d1), less_restrictive=not d0 and d1, resolution_level=level, evidence=evidence)

def require_pair(quantity0, quantity1, unit0, unit1, locator):
    assert quantity0 == quantity1, 'quantity substitution'
    assert unit0 == unit1, 'unit substitution'
    assert locator, 'missing source locator'
    return True

def main():
    start = time.perf_counter()
    ledger = []
    bindings = []
    comparisons = []
    (raw, path, h) = read('X8', 'PER_SITE_FULL_GEOMETRY_RISK.csv')
    bindings.append(dict(path=path, sha256=h))
    require_pair('canal_clearance', 'canal_clearance', 'mm', 'mm', 'sources/iti2018.pdf; TF2')
    seen = {}
    rows = []
    for r in raw:
        k = (r['case'], r['fdi'])
        values = tuple((float(r[x]) for x in ['old_K3_gap_mm', 'whole_cylinder_label_gap_lower_mm', 'whole_cylinder_label_gap_upper_mm']))
        if k in seen:
            assert np.allclose(values, seen[k], rtol=0, atol=1e-08)
            continue
        seen[k] = values
        (q0, lo, hi) = values
        assert (lo >= 2) == (hi >= 2), 'interval crosses threshold'
        a = decision(r['case'], r['fdi'], q0, lo, 2, 'clearance', 'PER_TOOTH', 'ANNOTATION_VIRTUAL_GEOMETRY')
        a['demo'] = 'X8'
        a['quantity'] = 'implant-body to annotated canal clearance'
        a['unit'] = 'mm'
        rows.append(a)
    ledger += rows
    x = summarize(rows)
    x.update(demo='X8', comparison='Old nominal 2 mm vs whole-body 2 mm annotated clearance', locator='sources/iti2018.pdf recommendations; TF2 labels', n_source_rows=len(raw), eligible_independent_units=len(rows), n_repeated_guide_rows=len(raw) - len(rows), excluded_original_sites=125, original_sites=2706, status='GEOMETRIC_DECISION_DISAGREEMENT; CLINICAL_BENEFIT_UNKNOWN')
    comparisons.append(x)
    require_pair('canal_clearance', 'canal_clearance', 'mm', 'mm', 'sources/iti2018.pdf; declared scenario')
    (raw, path, h) = read('X5', 'DECISION_MARGINS.csv')
    bindings.append(dict(path=path, sha256=h))
    rows = []
    for r in raw:
        if r['selection'] != 'ARCHIVED_RB' or r['region'] != 'canal_body':
            continue
        q = float(r['distance_mm'])
        b = float(r['scenario_budget_mm'])
        a = decision(r['case'], r['site'], q, q - b, float(r['threshold_mm']), 'clearance', 'PER_TOOTH', 'UNCALIBRATED_ERROR_BUDGET_SCENARIO')
        a.update(demo='X5', quantity='canal_body', unit='mm')
        rows.append(a)
    ledger += rows
    x = summarize(rows)
    x.update(demo='X5', comparison='Nominal 2 mm vs scenario-budget robust 2 mm', locator='sources/iti2018.pdf; X5 DECISION_MARGINS scenario_budget_mm', status='HYPOTHETICAL_ROBUSTNESS_CHANGE; PHYSICAL_BENEFIT_UNKNOWN', n_source_rows=len(raw), rejected_rows=len(raw) - len(rows), rejection_reason='other guard quantities or alternative selection')
    comparisons.append(x)
    require_pair('marginal_gap_after_explicit_proxy_mapping', 'marginal_gap_after_explicit_proxy_mapping', 'um', 'um', 'X13 FACIT.csv; historical rule UNVERIFIED_PRIMARY')
    (raw, path, h) = read('X13', 'measurements.csv')
    bindings.append(dict(path=path, sha256=h))
    rows = []
    reasons = collections.Counter()
    for r in raw:
        if r['region'] != 'marginal':
            reasons['different_quantity'] += 1
            continue
        if not r['marginal_spacer_um']:
            reasons['missing_explicit_marginal_setting'] += 1
            continue
        a = decision(r['source_family'], r['row_id'], float(r['marginal_spacer_um']), float(r['measured_mean_um']), 120, 'gap', 'PER_SURFACE_REGION', 'PUBLISHED_LAB_GROUP_MEAN')
        a.update(demo='X13', quantity='marginal_gap_group_mean', unit='um', cluster_type='source_family')
        rows.append(a)
    ledger += rows
    x = summarize(rows)
    x.update(demo='X13', comparison='Explicit CAD marginal setting-as-gap proxy vs measured marginal group mean <=120 um', locator='X13 FACIT.csv row DOI/table/cell; 120 um historical convention primary origin UNVERIFIED', status='COUNTERFACTUAL_PROXY_ONLY_NOT_VERIFIED_CURRENT_PRACTICE', n_source_rows=len(raw), rejected_rows=sum(reasons.values()), rejection_reasons=dict(reasons), cluster_meaning='source family, not patient', inference_scope='purposive laboratory groups; no specimen-tail acceptance or patient interval')
    comparisons.append(x)
    check = sum(((r['q0'] >= r['threshold']) != (r['q1'] >= r['threshold']) if r['kind'] == 'clearance' else (r['q0'] <= r['threshold']) != (r['q1'] <= r['threshold']) for r in ledger))
    assert check == sum((r['changed'] for r in ledger))
    faults = {}
    try:
        assert check + 1 == sum((r['changed'] for r in ledger))
    except AssertionError:
        faults['wrong_count_rejected'] = True
    try:
        require_pair('internal_gap', 'marginal_gap', 'um', 'um', 'source/table/cell')
    except AssertionError:
        faults['internal_as_marginal_rejected'] = True
    try:
        require_pair('clearance', 'clearance', 'mm', 'um', 'source/table/cell')
    except AssertionError:
        faults['unit_substitution_rejected'] = True
    try:
        require_pair('clearance', 'clearance', 'mm', 'mm', '')
    except AssertionError:
        faults['missing_locator_rejected'] = True
    faults['duplicate_patient_interval_invariant'] = bootstrap(ledger[:2581]) == bootstrap(ledger[:2581] * 2)
    faults['uniform_duplicate_cluster_invariant'] = abs(summarize(ledger[:2581], False)['changed_per100'] - summarize(ledger[:2581] * 2, False)['changed_per100']) < 1e-12
    assert all(faults.values())
    with (ROOT / 'raw/decision_ledger.csv').open('w') as f:
        public_rows = [dict(r, cluster_id=r['patient'], patient=r['patient'] if r['cluster_type'] == 'patient' else None) for r in ledger]
        w = csv.DictWriter(f, fieldnames=list(public_rows[0]))
        w.writeheader()
        w.writerows(public_rows)
    out = {'round': 'R1', 'claim_type': 'information_link', 'comparisons': comparisons, 'same_information_control': 'Exact agreement; independent scalar thresholds executed', 'fault_injections': faults, 'input_bindings': bindings, 'clinical_failure_coverage': 'UNKNOWN_FROM_OVERLAPPING_REVIEW_MARGINALS', 'cost': {'analysis_wall_seconds': time.perf_counter() - start, 'fit': 'none', 'inherited_extraction_cost': 'UNKNOWN; original demos not rerun'}}
    (ROOT / 'rounds/RESULTS_R1.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps([{k: v for (k, v) in x.items() if k in ['demo', 'n', 'changed', 'changed_per100', 'patients', 'patients_with_any_change', 'intervals']} for x in comparisons], indent=2))
if __name__ == '__main__':
    main()
