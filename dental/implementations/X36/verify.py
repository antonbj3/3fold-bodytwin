from pathlib import Path
import json, hashlib, time
import numpy as np
from meta_engine import fit
from analyze import read, dump, P, validate_unit
from contrasts import build_contrasts

def main():
    start = time.perf_counter()
    checks = []

    def check(name, passed, wrong_rejected, detail):
        checks.append(dict(name=name, pass_gate=bool(passed), injected_wrong_value_rejected=bool(wrong_rejected), detail=detail))
    for manifest in read('INPUT_MANIFEST.json'):
        content = (P / manifest['snapshot']).read_bytes()
        h = hashlib.sha256(content).hexdigest()
        check('input_sha_' + manifest['snapshot'], h == manifest['sha256'], hashlib.sha256(content + b'INJECTED').hexdigest() != manifest['sha256'], h)
    for name in ['PREREG_R1', 'PREREG_R2', 'PREREG_R3']:
        raw = (P / (name + '.json')).read_bytes()
        h = hashlib.sha256(raw).hexdigest()
        frozen = (P / (name + '.sha256')).read_text().split()[0]
        check('frozen_' + name, h == frozen, hashlib.sha256(raw + b'FAULT').hexdigest() != frozen, h)
    for path in sorted(P.glob('MODEL_*.json')):
        m = json.loads(path.read_text())
        c = m.get('dense_control')
        if not c:
            continue
        rel = c['variance_relative_max_error']
        obj = c['objective_at_same_theta_abs_error']
        b = np.array(m['beta'])
        berr = c['beta_max_abs_error'] / (1 + np.max(abs(b)))
        check('actual_table_dense_REML_' + path.stem, max(rel, obj, berr) < 1e-05, 1 + obj > 1e-05, c)
    data = read('NORMALIZED_DATA.json')
    for kind in ['cement', 'crown']:
        row = data[kind][0]
        fault = dict(row, unit='WRONG_UNIT')
        rejected = False
        try:
            validate_unit(fault)
        except ValueError:
            rejected = True
        check('declared_unit_' + kind, validate_unit(row), rejected, dict(expected_unit=row['unit'], injected_unit=fault['unit']))
    a = [r for r in data['cement'] if r['eligible']][:12]
    X = np.ones((len(a), 2))
    y = [r['outcome'] for r in a]
    S = np.diag([r['variance'] for r in a])
    bad = fit(y, X, S, [r['cluster'] for r in a], ['intercept', 'duplicate'], profile=False)
    check('rank_gate', bad['status'] == 'UNKNOWN_RANK_DEFICIENT', bad['status'] == 'UNKNOWN_RANK_DEFICIENT', bad['support'])
    r = read('RESULTS_R2.json')['controls']
    original = build_contrasts(data['cement'], data['crown'])
    faulty = [dict(x) for x in data['cement']]
    endpoint = next((x for x in original if x['dataset'] == 'cement'))['high_row_id']
    next((x for x in faulty if x['row_id'] == endpoint))['outcome'] += 100
    badcontrast = build_contrasts(faulty, data['crown'])
    wrong_cancel = max((abs(a['estimate'] - b['estimate']) for (a, b) in zip(original, badcontrast)))
    check('offset_scale_cancellation', r['offset_scale_cancellation']['pass_gate'], wrong_cancel > 1e-09, dict(correct=r['offset_scale_cancellation'], injected_nonuniform_offset_error=wrong_cancel))
    repair = read('SOURCE_CORRECTION_R2.json')
    exact = repair['corrected_mean_N']
    wrong = repair['old_mean_N']
    check('correct_source_cell_accepted', abs(exact - 2119.94) <= 0.011, abs(wrong - 2119.94) > 0.011, repair)
    check('wrong_prefixed_cell_rejected', r['wrong_prefixed_cell_rejected']['pass_gate'], abs(wrong - exact) > 0.011, r['wrong_prefixed_cell_rejected'])
    for c in r['endpoint_control']:
        check('contrast_reconstruction_' + c['stratum'], c['correct_pass'], c['injected_wrong_estimate_rejected'], c)
    fl = read('PRISMA_FLOW.json')
    for (kind, v) in fl.items():
        ok = v['unique_records'] == v['inherited_extracted_reports'] + v['not_extracted_or_eligibility_not_assessed']
        wrong = v['unique_records'] + 1 == v['inherited_extracted_reports'] + v['not_extracted_or_eligibility_not_assessed']
        check('flow_accounting_' + kind, ok, not wrong, v)
    for (kind, v) in read('FUNNEL_DATA.json').items():
        valid = v['egger']['formal_gate_pass']
        actual_k = len(v['study_summaries'])
        fake_k = sum((x['n_rows'] for x in v['study_summaries']))
        check('Egger_study_count_' + kind, valid == (actual_k >= 10), fake_k != actual_k, dict(correct=v['egger'], injected_row_count_as_study_count=fake_k))
    for (kind, v) in read('LOSO_R1.json').items():
        rows = [r for r in v['rows'] if r['status'] == 'SCORED']
        clusters = {r['held_cluster'] for r in rows}

        def rmse(shift):
            return np.sqrt(np.mean([np.mean([(r['prediction'] + shift - r['observed']) ** 2 for r in rows if r['held_cluster'] == s]) for s in clusters]))
        original = rmse(0)
        wrong = rmse(1000 if kind == 'cement' else 10)
        check('heldout_error_' + kind, abs(original - v['summary']['study_balanced_RMSE']) < 1e-07, wrong > 0.8 * v['summary']['practice_study_balanced_RMSE'], dict(original=float(original), injected=float(wrong)))
    ext = read('STATSMODELS_CONTROL.json')
    for c in ext['checks']:
        check('external_GLS_' + c['model'], c['pass_gate'], c['wrong_intercept_rejected'], c)
    fx = read('RESULTS_R3.json')
    check('physical_fixture_alias', fx['control']['direct_alias_pass'], fx['control']['injected_indenter_0_25mm_rejected'], fx['control'])
    out = dict(all_pass=all((c['pass_gate'] and c['injected_wrong_value_rejected'] for c in checks)), n_checks=len(checks), checks=checks, elapsed_seconds=time.perf_counter() - start, scope='Source integrity, flow arithmetic and numerical implementation. Scientific support gates are separate and fail. No independent scientific review.', known_limitations=['LINEAGE_DRIFT.json: exact first parent result bytes not retained before replay', 'Original source-corrupted row is intentionally rejected, not admitted', 'Source extractor and this verifier written by same producer; independent dual extraction pending'])
    dump('VERIFICATION.json', out)
    print('Verification', out['n_checks'], 'checks', out['all_pass'])
    assert out['all_pass']
if __name__ == '__main__':
    main()
