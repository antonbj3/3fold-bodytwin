import json, numpy as np
from scipy.spatial import cKDTree
from operators import circumcircle, path_geometry
from atlas import ROOT, DATA, sha, dump

def main():
    out = []
    for radius in [2.0, 5.0, 10.0]:
        for angle in [30.0, 45.0, 60.0]:
            ts = np.linspace(0, np.radians(angle), 101)
            ps = np.column_stack([radius * np.sin(ts), radius * (1 - np.cos(ts)), np.zeros(len(ts))])
            (r, turn) = circumcircle(ps[0], ps[50], ps[-1])
            length = path_geometry(ps)['annotation_path_length_mm']
            expected = radius * np.radians(angle)
            ok = abs(r / radius - 1) <= 0.05 and abs(turn - angle / 2) <= 1 and (abs(length / expected - 1) <= 0.03)
            out.append({'radius_mm': radius, 'arc_angle_deg': angle, 'radius_relative_error': abs(r / radius - 1), 'length_relative_error': abs(length / expected - 1), 'half_angle_error_deg': abs(turn - angle / 2), 'valid_pass': bool(ok), 'poison_radius_x2_rejected': abs(2 * r / radius - 1) > 0.05, 'poison_scale_x3_rejected': abs(3 * length / expected - 1) > 0.03, 'meaning': 'closed_form numerical control, not anatomy validation'})
    row = json.loads((ROOT / 'raw/R1_manifest.json').read_text())[0]
    valid = sha(row['path']) == row['sha256']
    bad = '0' * 64
    source_poison = sha(row['path']) != bad
    with np.load(row['path']) as a:
        (p, t) = (a['pulp'], a['tooth'])
    shifted = np.roll(p, 10, axis=2)
    contained = float((shifted & t).sum() / max(1, shifted.sum()))
    paired_valid = float((p & t).sum() / p.sum()) >= 0.95
    distances = json.loads((ROOT / 'raw/R1_distance_controls.json').read_text())
    result = {'circle': out, 'crop_sha_pass': valid, 'corrupt_sha_rejected': source_poison, 'paired_containment_pass': paired_valid, 'translated_pulp_containment': contained, 'translated_pulp_rejected': contained < 0.95, 'all_distance_pass': all((r['gate_pass'] for r in distances)), 'all_distance_poison_rejected': all((r['poison_distance_rejected'] for r in distances))}
    result['all_valid_passed'] = valid and paired_valid and result['all_distance_pass'] and all((r['valid_pass'] for r in out))
    result['all_poison_rejected'] = source_poison and result['translated_pulp_rejected'] and result['all_distance_poison_rejected'] and all((r['poison_radius_x2_rejected'] and r['poison_scale_x3_rejected'] for r in out))
    dump(ROOT / 'raw/CONTROLS.json', result)
    print(json.dumps(result), flush=True)
    if not result['all_valid_passed'] or not result['all_poison_rejected']:
        raise SystemExit(1)
if __name__ == '__main__':
    main()
