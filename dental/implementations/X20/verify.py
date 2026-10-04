"""Actual data controls plus corruptions; no clinical acceptance from a green check."""
import argparse
import copy
import hashlib
import json
import numpy as np
from measure import HERE, ROOT, DATA, sha, dump, check_frozen, scalar_decide, trace
from capacity_map import load_local, candidate

def hash_accepts(data, expected):
    return hashlib.sha256(data).hexdigest() == expected

def distribution_gate(value, n, reference, tolerance=0.2):
    return 'UNKNOWN_MINIMUM_N' if n < 20 else 'PASS' if abs(value - reference) <= tolerance else 'FAIL'

def run(artifacts=False):
    for round in (1, 2, 3, 4):
        check_frozen(f'PREREG_R{round}')
    for name in ('FROZEN_PREDICTIONS', 'FROZEN_PREDICTIONS_R2', 'FROZEN_PREDICTIONS_R3', 'FROZEN_PREDICTIONS_R4'):
        check_frozen(name)
    for f in ('DATA_MANIFEST.json', 'DATA_MANIFEST_R2.json', 'DATA_MANIFEST_R4.json'):
        m = json.loads((HERE / f).read_text())
        arr = m['local_crops'] if isinstance(m, dict) else m
        for r in arr:
            if 'path' in r:
                assert sha(r['path']) == r['sha256'], r['path']
    rows = [json.loads(l) for l in (HERE / 'COMBINED_SITES.jsonl').read_text().splitlines()]
    controls = []
    r = next((r for r in rows if r['valid'] and r.get('cross_read_pass')))
    reference = r['independent_voxel_read']['bone_height_mm']
    measured = r['bone_height_mm']
    accepts = lambda v: abs(v - reference) <= 0.15
    controls.append({'control': 'R1_voxel_height_reread', 'source': f"{r['case']}:{r['fdi']}", 'original_accepts': accepts(measured), 'injected_height_mm': measured + 2, 'injected_rejected': not accepts(measured + 2)})
    v = json.loads((HERE / 'MIRROR_HELDOUT_PREDICTIONS.json').read_text())[0]
    axis = np.asarray(v['true_axis'])
    truth = np.asarray(v['true_crest'])
    normal = np.cross(axis, [1, 0, 0])
    if np.linalg.norm(normal) < 0.01:
        normal = np.cross(axis, [0, 1, 0])
    normal /= np.linalg.norm(normal)
    error = lambda p: float(np.linalg.norm(p - truth - np.dot(p - truth, axis) * axis))
    controls.append({'control': 'R2_external_label_pose', 'source': f"{v['case']}:{v['fdi']}", 'original_accepts': error(truth) <= 2, 'injected_transverse_error_mm': error(truth + 5 * normal), 'injected_rejected': error(truth + 5 * normal) > 2, 'actual_mirror_gate': json.loads((HERE / 'SUMMARY_R2.json').read_text())['location_gate']})
    for r in [r for r in rows if r['valid']]:
        (local, anchor, axis, normal, tangent) = load_local(r)
        c = candidate(local, np.asarray(r['crest_mm']), axis, normal, tangent, 4.0, 6.0)
        if c.get('control_pass'):
            value = c['margins_mm']['bone']
            ref = c['trilinear_control_margins_mm']['bone']
            controls.append({'control': 'R3_independent_field_lookup', 'source': f"{r['case']}:{r['fdi']}", 'original_accepts': abs(value - ref) <= 0.3696153, 'injected_margin_mm': value + 2, 'injected_rejected': abs(value + 2 - ref) > 0.3696153})
            local.occ['sinus'] = np.zeros_like(local.occ['sinus'])
            empty = trace(local, anchor, axis, normal)
            controls.append({'control': 'R3_missing_sinus_is_unknown', 'original_accepts': True, 'injected_rejected': not empty['valid'], 'actual_reason': empty['reason']})
            break
    edges = [json.loads(l) for l in (HERE / 'RAW_R4_LOCAL_EDGES.jsonl').read_text().splitlines()]
    e = next((e for e in edges if e['valid'] and e['proxy_box_pass']))
    offset = e['image_minus_label_offset_mm']
    controls.append({'control': 'R4_image_edge_box', 'source': f"{e['case']}:{e['fdi']}", 'original_accepts': abs(offset) <= 0.3, 'injected_offset_mm': offset + 2, 'injected_rejected': abs(offset + 2) > 0.3})
    controls.append({'control': 'R4_gradient_vs_halflevel', 'source': f"{e['case']}:{e['fdi']}", 'original_accepts': abs(offset - e['gradient_max_offset_mm']) <= 0.3, 'injected_rejected': abs(offset + 2 - e['gradient_max_offset_mm']) > 0.3})
    ref = json.loads((HERE / 'EXTERNAL_REFERENTS.json').read_text())['anatomy']['premolar']['height_lt8']
    controls.append({'control': 'independent_distribution_gate', 'original_accepts': distribution_gate(ref, 20, ref) == 'PASS', 'injected_rejected': distribution_gate(ref + 0.3, 20, ref) == 'FAIL', 'actual_small_sample_is_unknown': distribution_gate(ref, 3, ref) == 'UNKNOWN_MINIMUM_N'})
    data = (HERE / 'PREREG_R1.json').read_bytes()
    h = (HERE / 'PREREG_R1.sha256').read_text().strip()
    controls.append({'control': 'freeze_hash', 'original_accepts': hash_accepts(data, h), 'injected_rejected': not hash_accepts(data + b' ', h)})
    catalog = json.loads((ROOT / 'results/DESIGN_implant/catalog.json').read_text())
    pairs = {(D, L) for (D, L, n) in catalog['catalog_pairs']}
    assert (4.0, 6.0) in pairs
    controls.append({'control': 'actual_catalog_pair', 'original_accepts': (4.0, 6.0) in pairs, 'injected_rejected': (3.0, 6.0) not in pairs})
    assert len(controls) >= 9 and all((c['original_accepts'] and c['injected_rejected'] for c in controls)), controls
    assert len({(r['case'], r['fdi']) for r in rows}) == len(rows), 'overlap counted twice'
    assert all((r['clinical_decision'] == 'UNKNOWN' for r in rows))
    assert sum((p.stat().st_size for p in DATA.glob('*') if p.is_file())) < 3000000000
    if artifacts:
        for r in json.loads((HERE / 'ARTIFACT_MANIFEST.json').read_text())['files']:
            assert sha(r['path']) == r['sha256'], r['path']
        result = json.loads((HERE / 'results.json').read_text())
        assert result['external_referent']['refutes_us'] is True
        assert (HERE / 'ANATOMY_DECISIONS.png').stat().st_size > 10000
        print('Artifact hash and figure checks PASS')
    else:
        dump(HERE / 'VERIFICATION.json', {'corruption_tests': controls, 'all_corruptions_rejected': True, 'frozen_r1_r2_r3_r4_verified': True, 'dataset_crops_verified': True, 'clinical_unknown_n': len(rows), 'scientific_or_clinical_acceptance': False})
        print(json.dumps(controls, indent=2))
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--artifacts', action='store_true')
    args = p.parse_args()
    run(args.artifacts)
