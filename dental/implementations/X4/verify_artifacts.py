"""Recompute frozen primary scores from raw vectors, allowing only proven rounding."""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '2'
import datetime, json, pathlib
import numpy as np
from freeze import ROOT, sha, write, now

def run():
    checks = []
    summaries = 0
    worst = 0.0
    casecount = 0
    for rid in ('R1', 'R2', 'R3', 'R4'):
        reg = ROOT / f'PREREG_{rid}.json'
        expected = (ROOT / f'PREREG_{rid}.sha256').read_text().split()[0]
        assert sha(reg) == expected
        result = json.loads((ROOT / f'RESULTS_{rid}.json').read_text())
        assert result['prereg_sha256'] == expected
        frozenfile = ROOT / f'FROZEN_PREDICTIONS_{rid}.json'
        frozen = json.loads(frozenfile.read_text())
        assert sha(frozenfile) == result['prediction_manifest_sha256']
        archive = ROOT / 'sources' / f'planner_{rid}.py' if rid != 'R3' else ROOT / 'sources/R3_code/implant_variant.py'
        assert sha(archive) == frozen['predictor_code_sha256']
        ts = datetime.datetime.fromisoformat(frozen['frozen_utc']).timestamp()
        rows = {r['case']: r for r in result['cases']}
        for rr in frozen['predictions']:
            assert sha(rr['path']) == rr['sha256']
            prediction = json.loads(pathlib.Path(rr['path']).read_text())
            casecount += 1
            if rid in ('R1', 'R4'):
                assert prediction['target_files_opened'] == ['Pre']
            if 'arrays' in prediction:
                assert sha(prediction['arrays']['path']) == prediction['arrays']['sha256']
            row = rows[rr['case']]
            if 'raw' not in row:
                continue
            assert sha(row['raw']['path']) == row['raw']['sha256']
            assert pathlib.Path(row['raw']['path']).stat().st_mtime >= ts
            with np.load(row['raw']['path']) as raw:
                for (name, method) in row['methods'].items():
                    fd = raw[name + '_forward']
                    rd = raw[name + '_reverse']
                    if method['primary_p95_mm'] is None:
                        assert not len(fd) or not len(rd)
                        continue
                    value = max(float(np.quantile(fd.astype(float), 0.95)), float(np.quantile(rd.astype(float), 0.95)))
                    spacing = max(float(np.spacing(fd).max()) if len(fd) else 0, float(np.spacing(rd).max()) if len(rd) else 0)
                    delta = abs(value - method['primary_p95_mm'])
                    bound = spacing / 2 + np.finfo(float).eps * max(1, value) * 8
                    assert delta <= bound, (rid, rr['case'], name, delta, bound)
                    worst = max(worst, delta)
                    summaries += 1
        checks.append(dict(round=rid, prereg_hash_matches=True, prediction_hashes_match=True, prediction_code_archive_matches=True, raw_hashes_match=True, frozen_before_raw_measurement=True, n_cases=len(rows), n_available=sum((r['numerically_available'] for r in rows.values()))))
    report = json.loads((ROOT / 'results.json').read_text())
    feedback = json.loads((ROOT / 'GRAPH_FEEDBACK.json').read_text())
    assert feedback['sha256'] == sha(ROOT / 'results.json')
    assert feedback['review_state'] == 'PENDING_INDEPENDENT_REVIEW'
    assert not report['scientific_admission']
    assert json.loads((ROOT / 'CONTROL_FALSIFICATION.json').read_text())['passed']
    write(ROOT / 'FINAL_VERIFICATION.json', dict(passed=True, verified_utc=now(), rounds=checks, n_case_runs=casecount, n_raw_method_summaries_recomputed=summaries, max_quantized_raw_vs_report_difference_mm=worst, tolerance='derived float32 half-ULP plus double arithmetic bound, not a retuned empirical gate', source_roles='virtual Post and defective Original retained', graph_result_hash_matches=True, review_state='PENDING_INDEPENDENT_REVIEW', scientific_admission=False))
    print('VERIFIED', casecount, 'case runs', summaries, 'raw summaries; worst roundtrip', worst)
if __name__ == '__main__':
    run()
