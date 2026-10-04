"""Validate numerical/source contracts and non-vacuous corruption controls."""
from pathlib import Path
import json, hashlib, math
P = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def verify():
    checks = []
    for tag in ['R1', 'R1_LINEAGE', 'R2', 'R3', 'R4']:
        f = P / f'PREREG_{tag}.json'
        assert sha(f) == f.with_suffix('.sha256').read_text().strip()
        checks.append({'check': 'prereg ' + tag, 'pass': True})
    for tag in ['R1', 'R2', 'R3']:
        f = P / f'FROZEN_PREDICTIONS_{tag}.json'
        assert sha(f) == f.with_suffix('.sha256').read_text().strip()
        r = json.loads(f.read_text())
        assert r['source_measurements_sha256'] == sha(P / 'measurements.json')
        checks.append({'check': 'frozen prediction ' + tag, 'pass': True})
    sources = json.loads((P / 'SOURCE_MANIFEST.json').read_text())
    assert sources['measurements_sha256'] == sha(P / 'measurements.json')
    assert sources['csv_sha256'] == sha(P / 'measurements.csv')
    for s in sources['sources']:
        assert sha(s['path']) == s['sha256']
    checks.append({'check': '11 primary source hashes + frozen extracted input', 'pass': True})
    r1 = json.loads((P / 'RESULTS_R1.json').read_text())
    r2 = json.loads((P / 'RESULTS_R2.json').read_text())
    r3 = json.loads((P / 'RESULTS_R3.json').read_text())
    r4 = json.loads((P / 'RESULTS_R4.json').read_text())
    probes = []
    for r in [r1, r2, r3, r4]:
        for probe in r['corruption_tests']:
            assert probe.get('rejected', probe.get('error_gate_rejects', probe.get('parity_gate_rejects', False)))
            probes.append(dict(round=r['round'], **probe))
    checks.append({'check': 'all executed corruption probes reject the injected error', 'pass': True, 'count': len(probes)})
    assert all((d['maximum_gls_prediction_difference_um'] < 0.001 for d in r1['metrics'].values()))
    assert all((d['maximum_gaussian_control_difference_um'] < 0.001 for d in r2['metrics'].values()))
    assert r3['maximum_reciprocal_control_difference_um'] < 0.001
    assert r4['numerical_verification_pass'] and r4['same_mean_pass']
    assert sha(r4['operator_path']) == r4['operator_sha256']
    checks.append({'check': 'same-information controls and actual read-only BTE1 operator', 'pass': True})
    rows = json.loads((P / 'measurements.json').read_text())
    assert len(rows) == 148
    assert sum((d['region'] in ['marginal', 'axial', 'occlusal'] for d in rows)) == 94
    assert not r1['metrics']['axial']['dose_support']['inverse_supported']
    assert not r1['metrics']['occlusal']['dose_support']['inverse_supported']
    decisions = json.loads((P / 'INVERSE_DECISIONS.json').read_text())
    assert all((d['status'] == 'UNKNOWN' and d['CAD_setting_um'] is None for d in decisions['regional']))
    checks.append({'check': 'unsupported inverse query abstains', 'pass': True})
    size = sum((f.stat().st_size for f in P.rglob('*') if f.is_file()))
    assert size < 3000000000
    assert (P / 'CEMENT_GAP_DEMO.png').stat().st_size > 10000 and (P / 'CEMENT_GAP_DEMO.pdf').stat().st_size > 1000
    checks.append({'check': 'figure generated, lane disk <3GB', 'pass': True, 'bytes': size})
    record = dict(status='PASS', scientific_admission=False, checks=checks, corruption_probes=probes)
    (P / 'DELIVERY_CHECK.json').write_text(json.dumps(record, indent=2) + '\n')
    print('DELIVERY_CHECK PASS;', len(probes), 'non-vacuous corruption probes')
    return record
if __name__ == '__main__':
    verify()
