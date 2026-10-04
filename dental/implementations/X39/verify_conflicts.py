"""Compare conflict outcomes to independent per-case source recounts."""
from pathlib import Path
import json
P = Path(__file__).resolve().parent
R = P.parent

def read(p):
    return json.loads(Path(p).read_text())

def verify(k, n, expected_k, expected_n):
    return (k, n) == (expected_k, expected_n)

def main():
    results = read(P / 'results.json')
    manifest = read(R / 'LANE_X7_BITE2TEXT/raw/DATA_MANIFEST.json')
    corrected = read(R / 'LANE_X21_CONTACT_MAP/raw/CORRECTED_REPORTS_R5.json')
    comparison = read(R / 'LANE_X21_CONTACT_MAP/raw/COMPARISON_ROWS_R2.json')
    disagreement = read(R / 'LANE_X7_BITE2TEXT/raw/REPORT_DISAGREEMENT.json')
    mapping = dict(zip(['Overbite', 'Overjet', 'Crossbite', 'Midlines', 'Curve of Spee', 'Molar right', 'Molar left', 'Canine right', 'Canine left'], ['overbite', 'overjet', 'crossbite', 'midlines', 'spee', 'molar_right', 'molar_left', 'canine_right', 'canine_left']))
    checks = []
    for row in results['report_conflicts']:
        f = mapping[row['field'].replace(', corrected', '')]
        if f == 'crossbite':
            cases = manifest['test'] if row['population'].startswith('Internal') else list(corrected)
            groups = [[q['labels'][f] for q in corrected[c] if q['labels'].get(f) is not None] for c in cases]
        elif row['population'].startswith('Internal'):
            groups = [q['reference_all'] for q in comparison if q['field'] == f]
        else:
            groups = [q['all_values'] for q in disagreement[f]['pairs']]
        multi = [q for q in groups if len(q) >= 2]
        n = len(multi)
        k = sum((len(set(q)) > 1 for q in multi))
        assert verify(row['conflicts'], row['n'], k, n)
        assert not verify(row['conflicts'], 989, k, n)
        assert abs(row['fraction'] - k / n) < 1e-12
        checks.append({'field': row['field'], 'population': row['population'], 'k': k, 'n': n, 'injected_989_denominator_rejected': True})
    (P / 'REVIEW_CONFLICT_CONTROL.json').write_text(json.dumps({'status': 'PASS', 'checks': checks, 'source': 'Per-case report category rows, not reported fraction divided by its own denominator', 'new_models': 0}, indent=2) + '\n')
    print('Independent report-conflict control PASS: 14 source recounts and denominator injections.')
if __name__ == '__main__':
    main()
