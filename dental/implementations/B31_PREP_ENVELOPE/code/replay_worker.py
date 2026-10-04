from common import *
import importlib.metadata as metadata
from run_r2 import solve2
from run_r3 import main as run_paths
from run_r4 import main as run_exact_surfaces
from controls import run_controls
from report import main as report

def verify_exterior(rows):
    cohort = {r['key']: r for r in read(ROOT / 'FROZEN_COHORT.json')['records']}
    checks = []
    for row in rows:
        if 'candidate_path' not in row:
            continue
        d = load_np(row['candidate_path'])
        t = load_np(cohort[row['key']]['private_path'])['target']
        e = d['vertices'][d['faces'][d['roles'] == 0]]
        canonical = lambda a: [sorted((tuple((float(x) for x in p)) for p in tr)) for tr in a]
        eq = canonical(t) == canonical(e)
        tampered = e.copy()
        tampered[0, 0, 0] += 0.2
        reject = canonical(t) != canonical(tampered)
        checks.append({'key': row['key'], 'original_triangles_identical': eq, 'injected_0.2mm_exterior_shift_rejected': reject, 'identity_error_mm': 0.0 if eq else 'MISMATCH', 'resolution': 'PER_POINT'})
    assert all((c['original_triangles_identical'] and c['injected_0.2mm_exterior_shift_rejected'] for c in checks))
    dump(ROOT / 'raw/EXTERIOR_IDENTITY_CONTROLS.json', checks)

def verify_frozen_detail():
    frozen = read(ROOT / 'FROZEN_BASELINE_DETAIL.json')
    checks = []
    for (rd, want) in frozen['rounds'].items():
        actual = {r['key']: r for r in read(ROOT / 'raw' / (rd + '_RESULTS.json'))}
        assert set(actual) == {r['key'] for r in want}
        for w in want:
            r = actual[w['key']]
            values = {k: r.get('canonical_intaglio', {}).get('exact_volume_mm3') if k == 'exact_volume_mm3' else r.get(k) for k in w if k != 'key'}
            expected = {k: v for (k, v) in w.items() if k != 'key'}
            checks.append({'round': rd, 'key': w['key'], 'all_values_identical': values == expected, 'differing_fields': [k for k in expected if values[k] != expected[k]]})
    ok = all((r['all_values_identical'] for r in checks))
    dump(ROOT / 'raw/FROZEN_DETAIL_COMPARISON.json', {'all_pass': ok, 'rows': checks, 'scope': frozen['scope']})
    if not ok:
        raise RuntimeError('Frozen case-level artifacts changed; preserve failed replay')
    return ok

def main(full_source=False):
    start = time.perf_counter()
    lock = read(ROOT / 'RELEASE_LOCK.json')
    for (p, expect) in lock['code_and_contracts'].items():
        if sha(ROOT / p) != expect:
            raise RuntimeError('Executed-copy code/contract hash mismatch ' + p)
    cache = read(ROOT / 'FROZEN_CACHE.json')
    bykey = {r['key']: r for r in cache['rows']}
    expected = read(ROOT / 'BASELINE_COMPARISON.json')
    rows = []
    for r in read(ROOT / 'FROZEN_COHORT.json')['records']:
        assert sha(r['private_path']) == r['private_sha256'] and sha(r['public_path']) == r['public_sha256']
        try:
            c = bykey.get(r['key'])
            rr = solve2(r, None if full_source else c)
        except Exception as e:
            rr = {'key': r['key'], 'family': r['family'], 'status': 'UNKNOWN_CONSTRUCTION', 'reason': repr(e), 'complete_digital_status': 'UNKNOWN', 'resolution': 'PER_TOOTH'}
            cp = ROOT / 'raw/r2' / (r['key'] + '_FIELD_CHECKPOINT.json')
            if cp.exists():
                rr.update(read(cp))
        rows.append(rr)
        dump(ROOT / 'raw/R2_RESULTS.json', rows)
        print('replay', rr['key'], rr['status'], flush=True)
    dump(ROOT / 'raw/R2_COST.json', {'seconds': time.perf_counter() - start, **usage(), 'source_preprocessing': 'FULL_RECOMPUTED' if full_source else 'frozen full-box sign/distance fields; source hashes and all dictionaries regenerated; failed-family fields also consumed', 'legacy_acquisition': 'UNKNOWN'})
    verify_exterior(rows)
    run_controls()
    run_paths()
    run_exact_surfaces()
    detail_ok = verify_frozen_detail()
    out = report()
    observed = {k: out['summaries'][k] for k in ['R2', 'R3', 'R4']}
    ok = observed == expected['summary']
    dump(ROOT / 'raw/REPLAY_COMPARISON.json', {'matches_frozen_summary': ok, 'expected': expected['summary'], 'observed': observed, 'all18_reconstructed': len(rows) == 18, 'all54_case_artifacts_and_values_identical': detail_ok, 'scope': 'R2 all18 plus R3 full paths, R4 exact positive surface-pair retest and falsifiers; R1 is preserved historical failed branch, not freshly regenerated', 'full_source': full_source, 'source_fields_sha_checked': True, 'copied_code_hashes_checked': True, 'seconds': time.perf_counter() - start, **usage()})
    if not ok:
        raise RuntimeError('Frozen summary changed; preserve failed replay')
    state('REPLAY_COMPLETE', observed, 'Independent review; next actual marginal-volume / exact shoulder Boolean construction')
    print('REPLAY_COMPLETE', ROOT, flush=True)
if __name__ == '__main__':
    import sys
    main('--full-source' in sys.argv)
