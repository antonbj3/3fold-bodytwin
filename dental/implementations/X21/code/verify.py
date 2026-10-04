"""Source-bound scientific reproducibility and fault rejection; no clinical admission."""
from dental_release.paths import expand as _release_expand
import json, zipfile, time, copy, tempfile
from pathlib import Path
import numpy as np
from contact import P, D, sha, write, load_case, prepare, boxes, buckets, extrema

def run():
    st = time.perf_counter()
    checks = []

    def check(name, passed, **kw):
        checks.append(dict(name=name, passed=bool(passed), **kw))
        if not passed:
            raise AssertionError(name)
    for k in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']:
        p = P / f'PREREG_{k}.json'
        check('prereg_' + k, sha(p) == Path(str(p) + '.sha256').read_text().strip())
    m = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    rows = {c: json.load(open(P / 'raw/cases' / (c + '.json'))) for c in m['cases']}
    check('993_nonempty_complete', len(rows) == 993 and all((r.get('error') is None for r in rows.values())))
    check('explicit_empty_source_refusal', any((r['case_id'] == _release_expand('@DENTAL_CASE_ID@') and (not r['complete']) for r in m['excluded'])))
    for name in ['FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R2.json']:
        freeze = json.load(open(P / name))
        check(name + '_probabilities', sha(Path(freeze['predictions_path'])) == freeze['predictions_sha256'])
        reference = P / 'raw/failed_checks/report_capability_R1_frozen.py' if name == 'FROZEN_PREDICTIONS.json' else P / 'code/report_capability.py'
        check(name + '_code', sha(reference) == freeze['code_sha256'])
        check(name + '_case_inputs', all((sha(P / 'raw/cases' / (r['case_id'] + '.json')) == r['case_file_sha256'] and rows[r['case_id']]['map_sha256'] == r['map_sha256'] for r in freeze['map_manifest'])))
    freeze = json.load(open(P / 'FROZEN_PREDICTIONS_R5.json'))
    check('R5_frozen_prediction', sha(P / 'raw/PREDICTIONS_R5.json') == freeze['prediction_sha256'])
    check('R5_frozen_source', sha(P / 'raw/CORRECTED_REPORTS_R5.json') == freeze['corrected_source_reports_sha256'])
    check('R5_frozen_code', sha(P / 'code/repair_reports.py') == freeze['code_sha256'])
    manifest = json.load(open(P / 'raw/DATA_ARTIFACT_MANIFEST.json'))
    check('data_files_hashes', all((Path(r['path']).stat().st_size == r['bytes'] and sha(Path(r['path'])) == r['sha256'] for r in manifest['files'])))
    check('3GB_budget', sum((p.stat().st_size for p in D.rglob('*') if p.is_file())) < 3000000000)
    with zipfile.ZipFile(m['zip']) as z:
        import hashlib
        central = hashlib.sha256(json.dumps([(i.filename, i.file_size, i.CRC) for i in z.infolist()]).encode()).hexdigest()
        check('source_zip_central_directory', central == m['zip_central_directory_sha256'])
        c = m['numerical_panel'][0]
        (meta, R, center, arches) = load_case(z, c)
        (U, ul, uc, uf, us) = arches[0]
        (L, ll, lc, lf, ls) = arches[1]
        ub = boxes(U)
        lb = boxes(L)
        origin = np.minimum(ub[:, :2].min(0), lb[:, :2].min(0)) - 0.8
        hi = np.maximum(ub[:, 3:5].max(0), lb[:, 3:5].max(0))
        shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
        (ptr, ids) = buckets(lb, origin, 0.8, shape)
        (best, wi, xy, n, nt) = extrema(U, L, ul, ll, uc, lc, ub, lb, origin, 0.8, shape, ptr, ids, True)
        saved = np.load(rows[c]['map_path'])['minimum_gap']
        finite = np.isfinite(saved)
        err = float(np.max(np.abs(saved[finite] - best[finite])))
        check('live_full_case_minima', np.array_equal(finite, np.isfinite(best)) and err <= 1e-07, max_error_mm=err, case_id=c)
        bad = copy.deepcopy(meta['arches']['upper'])
        bad['input_sha256'] = '0' * 64
        rejected = False
        try:
            prepare(z.read(c + '/ios/ios_upper.stl'), bad, 'upper', R, center)
        except ValueError as e:
            rejected = str(e) == 'SOURCE_HASH_MISMATCH'
        check('wrong_source_hash_rejected', rejected)
    for k in ['R3', 'R4']:
        rr = json.load(open(P / 'raw' / f'RESULTS_{k}.json'))
        check(k + '_wrong_bound_controls', rr['all_wrong_lower_bound_controls_rejected'] if k == 'R3' else rr['all_wrong_bounds_rejected'])
    r6 = json.load(open(P / 'raw/RESULTS_R6.json'))
    check('125_negative_scope_injections_rejected', r6['literal_no_lateral_controls'] == 125 and r6['all_wrong_region_labels_rejected'])
    from compare_measurement import compare
    pf = json.load(open(P / 'raw/MEASUREMENT_PREDICTIONS.json'))
    p = pf['records'][0]
    bad = copy.deepcopy(p)
    bad['xyz_mm'][0] += 1
    probe = P / 'raw/INJECTED_POINT_ERROR_1MM.json'
    write(probe, dict(records=[bad]))
    out = P / 'raw/INJECTED_POINT_COMPARISON.json'
    compare(probe, out)
    r = json.load(open(out))
    check('future_point_reader_1mm_error_rejected', r['gate'] == 'FAIL' and abs(r['records'][0]['error_mm'] - 1) < 1e-12, external_truth='our_own_fixture; detector check only, no anatomy validation')
    rawr1 = json.load(open(P / 'raw/failed_checks/RESULTS_R1_INITIAL.json'))
    current = json.load(open(P / 'raw/RESULTS_R1.json'))
    check('failed_R1_gate_preserved', current['gates'] == rawr1['gates'] and current['primary_gate'] == 'FAIL')
    check('failed_R5_global_scope_gate_preserved', not json.load(open(P / 'raw/RESULTS_R5.json'))['gates']['explicit_negation_injection'])
    write(P / 'raw/VERIFICATION.json', dict(passed=True, checks=checks, wall_s=time.perf_counter() - st, independent_review=False, clinical_validation='UNKNOWN'))
    print('verified', len(checks), 'checks, live case error', err, 's', round(time.perf_counter() - st, 2))
if __name__ == '__main__':
    run()
