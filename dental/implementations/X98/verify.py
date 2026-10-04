"""Read-only final gate; validates raw numbers, not producer PASS booleans."""
import argparse
import csv
import json
import math
from pathlib import Path
from implant_safety import SafetyModule
from implant_safety.module import sha, load
from validation import check_distance, controls as recompute_controls
ROOT = Path(__file__).resolve().parent

def verify(output):
    module = SafetyModule(ROOT)
    fp = ROOT / 'FROZEN_PREDICTIONS_R2.json'
    assert sha(fp) == fp.with_suffix('.json.sha256').read_text().strip(), 'frozen prediction hash drift'
    if (ROOT / 'CODE_MANIFEST.json').exists():
        for r in load(ROOT / 'CODE_MANIFEST.json'):
            assert sha(ROOT / r['path']) == r['sha256'], 'code hash drift: ' + r['path']
    with (output / 'TABLE.csv').open() as f:
        tab = list(csv.DictReader(f))
    frozen = load(fp)['rows']
    assert len(tab) == len(frozen) == 176
    for (actual, expected) in zip(tab, frozen):
        for (k, v) in expected.items():
            assert actual[k] == str(v), 'frozen table field changed: ' + k
    rows = [json.loads(l) for l in (output / 'QUERIES.jsonl').read_text().splitlines() if l.strip()]
    assert len(rows) == 176
    by = {}
    for (r, t) in zip(rows, tab):
        q = r['query']
        site = r['site_id']
        g = r['terms']['guide_combination']
        c = r['combined_conditional']
        geo = r['digital_geometry']
        expected_query = module.query(site_id=site, **q)
        assert r == expected_query, 'raw query differs from source-bound recomputation: ' + site
        assert site == t['site_id'] and q['guide_type'] == t['guide'] and (q['implant_system'] == t['system'])
        assert r['physical_safety']['status'] == 'UNKNOWN' and r['physical_safety']['interval_mm'] is None and (r['physical_safety']['injury_probability'] is None)
        assert geo['nominal']['interval_mm'][0] == float(t['nominal_gap_mm'])
        assert geo['tool_envelope']['interval_mm'][0] == float(t['swept_gap_mm'])
        assert g['entry_mm'] == float(t['guide_entry_mm']) and g['apex_mm'] == float(t['guide_apex_mm'])
        assert g['angle_deg'] == float(t['guide_angle_deg']) and g['rotation_mm'] == float(t['guide_rotation_mm'])
        assert g['combined_displacement_mm'] == float(t['guide_combined_mm'])
        assert abs(g['combined_displacement_mm'] - (max(g['entry_mm'], g['apex_mm']) + g['rotation_mm'])) <= 1e-10
        assert c['interval_mm'] == [float(t['conditional_gap_lower_mm']), float(t['conditional_gap_upper_mm'])]
        assert c['class_vs_2mm'] == t['conditional_class']
        for (key, L) in [('nominal', q['length_mm']), ('revision_union', q['length_mm']), ('tool_envelope', q['length_mm'] + r['terms']['drill_extra_depth']['maximum_or_example_mm'])]:
            b = geo[key]['interval_mm']
            assert 0 <= b[1] - b[0] <= 1e-06
            assert geo[key]['rigorous_IEEE_enclosure'] == 'MISSING'
            by[site, key, L] = b
    ct = load(output / 'CONTROLS.json')
    assert len(ct['geometry']) == 88 and len(ct['guide']) == 16
    for v in ct['geometry']:
        b = by[v['site_id'], v['quantity'], v['length_mm']]
        assert b == v['candidate_interval_mm'], 'control candidate changed'
        check_distance(b, v['control_mm'])
        c = v['independent_control_enclosure']
        (lo, hi) = c['interval_mm']
        assert -1e-10 <= hi - lo <= 1e-06 and c['primal_feasibility_error'] <= 1e-10
        assert hi == v['control_mm']
        assert lo <= b[1] + 1e-06 and hi >= b[0] - 1e-06
        assert v['injected_plus1mm_rejected'] is True
    fresh_guide = recompute_controls(module, rows)['guide']
    assert len(ct['guide']) == len(fresh_guide)
    for (v, expected) in zip(ct['guide'], fresh_guide):
        assert v['guide'] == expected['guide']
        assert v.get('extra_mm') == expected.get('extra_mm')
        assert math.isfinite(v['error_mm']) and 0 <= v['error_mm'] <= 1e-10
        assert abs(v['error_mm'] - expected['error_mm']) <= 1e-10 and v['injected_plus1mm_rejected'] is True
        if 'independent_tail_bound' in expected:
            tail = v['independent_tail_bound']
            assert math.isfinite(tail) and 0 <= tail <= 0.05 + 1e-10
            assert abs(tail - expected['independent_tail_bound']) <= 1e-10
    assert ct['control_gate'] == 'PASS'
    assert len(load(output / 'CONTRACT_TESTS.json')) == 19
    assert all((r['result'] == 'PASS' for r in load(output / 'CONTRACT_TESTS.json')))
    suf = load(output / 'SUFFICIENCY.json')
    assert suf['identity_error_mm'] == 0 and suf['downstream_difference_mm'] == 1.0
    assert len(load(output / 'SOURCE_REPLAY.json')) == 44
    assert max((r['replay_max_error_mm'] for r in load(output / 'SOURCE_REPLAY.json'))) <= 1e-06
    return dict(status='PASS', raw_queries=176, raw_geometry_controls=88, source_replay_rows=44, physical_certificates=0)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output-dir', type=Path, default=ROOT / 'raw')
    a = p.parse_args()
    print(json.dumps(verify(a.output_dir)))
