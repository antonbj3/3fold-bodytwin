"""Hash current ZIP/cache bytes and recompute the headline from current tooth poses."""
import hashlib, json, math, os
from pathlib import Path
P = Path(__file__).resolve().parent

def read(name):
    return json.loads((P / name).read_text())

def hashes(path):
    s = hashlib.sha256()
    m = hashlib.md5()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda : f.read(1048576), b''):
            s.update(chunk)
            m.update(chunk)
    return (s.hexdigest(), m.hexdigest())

def transformed(row):
    H = row['H']['bilateral_molars']
    c = row['centroid_base_mm']
    return [sum((H[i][j] * c[j] for j in range(3))) + H[i][3] for i in range(3)]

def distance(a, b):
    return math.sqrt(sum(((x - y) ** 2 for (x, y) in zip(a, b))))

def verify_integrity():
    checks = []
    for r in read('inputs/ZIP_INTEGRITY.json')['files']:
        (s, m) = hashes(r['file'])
        good = s == r['sha256'] and m == r['md5'] and ('md5:' + m == r['published_checksum']) and (Path(r['file']).stat().st_size == r['bytes'])
        checks.append(dict(file=r['file'], pass_gate=good, injected_wrong_hash_rejected=s != '0' * 64))
    manifest = read('DATA_MANIFEST.json')
    records = manifest if isinstance(manifest, list) else manifest['files']
    for r in records:
        path = Path(r.get('path', r.get('file', '')))
        (s, _) = hashes(path)
        checks.append(dict(file=str(path), pass_gate=s == r['sha256'] and path.stat().st_size == r['bytes'], injected_wrong_hash_rejected=s != '0' * 64))
    if not checks or not all((x['pass_gate'] and x['injected_wrong_hash_rejected'] for x in checks)):
        raise ValueError('Current ZIP/cache bytes disagree with frozen hashes')
    return checks

def main():
    checks = verify_integrity()
    poses = read('raw/R1_POSES.json')
    index = {(r['patient'], r['jaw'], r['source'], r['week'], r['fdi']): r for r in poses}
    count = 0
    maximum = 0.0
    max_error = 0.0
    by_patient = {}
    pairs = read('raw/R2_PAIRS.json')
    for pair in pairs:
        (p, j, t, a, b) = (pair[k] for k in ['patient', 'jaw', 'week', 'fdi_i', 'fdi_j'])
        req = [index[p, j, s, w, f] for s in ['Sirona', 'Bottmedical'] for w in [t, 0] for f in [a, b]]

        def d(s, w):
            return distance(transformed(index[p, j, s, w, a]), transformed(index[p, j, s, w, b]))
        tracking = abs(d('Sirona', t) - d('Sirona', 0) - (d('Bottmedical', t) - d('Bottmedical', 0)))
        budget = 0.8 + sum((r['sampling_spread_mm'] for r in req))
        good = all((r['surface_fit_pass'] for r in req + [index[p, j, 'Sirona', w, f] for w in [t, t - 1] for f in [a, b]]))
        resolved = good and tracking > budget
        count += resolved
        by_patient[p] = by_patient.get(p, 0) + resolved
        maximum = max(maximum, tracking)
        max_error = max(max_error, abs(tracking - pair['relative_tracking_error_mm']), abs(budget - pair['tracking_sensitivity_budget_mm']))
    primary = read('results.json')['primary_finding']
    main_pass = count == primary['pair_tracking_rows_above_declared_sensitivity_budget'] and len(pairs) == primary['pair_week_rows'] and (len({r['patient'] for r in poses}) == primary['independent_patients']) and (abs(maximum - primary['max_pair_tracking_error_mm']) < 1e-09) and (max_error < 1e-09)
    out = dict(status='PASS' if checks and all((x['pass_gate'] and x['injected_wrong_hash_rejected'] for x in checks)) and main_pass else 'FAIL', scope='Present bytes and derived relative-distance arithmetic; n=2 pipeline, no clinical or force validation.', integrity_checked_files=len(checks), checks=checks, primary_pass=main_pass, recomputed_pair_rows=len(pairs), recomputed_resolved=count, by_patient=by_patient, recomputed_max_mm=maximum, max_raw_error_mm=max_error, injected_count_rejected=count != primary['pair_tracking_rows_above_declared_sensitivity_budget'] + 1)
    (P / 'REVIEW_CONTROL_RESULT.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: v for (k, v) in out.items() if k != 'checks'}, indent=2))
    if out['status'] != 'PASS':
        raise SystemExit(1)
if __name__ == '__main__':
    main()
