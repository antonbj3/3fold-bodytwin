from common import *
from fractions import Fraction as F
import numpy as np

def main():
    gates = []
    for p in ROOT.glob('PREREG_*.json'):
        assert sha(p) == p.with_suffix(p.suffix + '.sha256').read_text().split()[0], p
    for p in ROOT.glob('FROZEN_PREDICTIONS*.json'):
        assert sha(p) == p.with_suffix(p.suffix + '.sha256').read_text().split()[0], p
    gates.append('prereg_and_predictions_unchanged')
    for x in read(ROOT / 'SOURCE_LOCK.json')['files']:
        assert sha(x['path']) == x['sha256'], x['path']
    gates.append('original_source_locks')
    r1 = read(ROOT / 'raw/R1_RESULTS.json')
    r2 = read(ROOT / 'raw/R2_RESULTS.json')
    r3 = read(ROOT / 'raw/R3_RESULTS.json')
    r4 = read(ROOT / 'raw/R4_RESULTS.json')
    for x in [r1['point_file'], r2['artifact']]:
        assert sha(x['path']) == x['sha256']
    gates.append('all_array_hashes')
    distances = list(r2['pulp_queries']) + [d for q in r2['canal_queries'] for d in q['distances'].values()]
    for d in distances:
        squared = F(d['squared_mm2'])
        (lo, hi) = map(F.from_float, d['distance_interval_mm'])
        assert lo * lo <= squared <= hi * hi and 0 <= lo <= hi
        if squared > 0:
            assert max(F(0), lo - F(1, 100)) ** 2 < squared
    gates.append('exact_rational_distance_enclosures_and_wrong_values')
    assert r2['control_max_error_mm'] <= 1e-08 and r2['controls_count'] == 159
    assert r3['summary_identity_error'] == 0 and r3['strongest_witness']['downstream_difference_lower_mm'] >= 0.01
    assert r3['enclosure_violations'] == 0 and r3['common_rigid_invariance_error'] == 0
    for c in r3['comparisons']:
        assert c['summary_squared_norm_mm2'] == ['9/100', '9/100']
        for k in ('negative_distance', 'positive_distance'):
            assert c['predicted_interval_mm'][0] <= c[k][0] <= c[k][1] <= c['predicted_interval_mm'][1]
    gates.append('summary_and_every_directed_state_enclosure')
    assert all((f['rejected'] for f in read(ROOT / 'raw/CONTROLS_R1.json')))
    for r in (r2, r3, r4):
        assert all((f['rejected'] for f in r['faults']))
    for q in r4['queries']:
        assert q['support_voxels'] == 27
        (lo, hi) = q['conditional_gray_interval']
        assert all((lo <= w['gray'] <= hi for w in q['witnesses']))
        assert q['physical_gray_interval'] is None
    gates.append('raw_gray_support_and_all_fault_rejections')
    assert all((q['status'] == 'UNKNOWN' and q['nominal_value'] is None for q in r1['decisions']))
    points = read(ROOT / 'raw/X77_POINT_DECISIONS.json')
    assert len(points) == 101 and all((p['patient_id'] == 'LiuHao2023Demo1' and p['status'] == 'UNKNOWN' and (p['pulp_distance_mm'] is None) for p in points))
    gates.append('no_validation_patient_substitution')
    if (ROOT / 'CODE_LOCK.json').exists():
        for x in read(ROOT / 'CODE_LOCK.json')['files']:
            assert sha(x['path']) == x['sha256']
        gates.append('final_code_lock')
    dump(ROOT / 'VERIFICATION.json', dict(status='PASS', gates=gates, review_state='PRODUCER_VERIFICATION_NOT_INDEPENDENT_REVIEW', data_bytes=diskcheck()))
    print('Verification PASS:', len(gates), 'gate families')
if __name__ == '__main__':
    main()
