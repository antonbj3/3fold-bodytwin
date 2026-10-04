"""Reopen actual labels and plans; verify all retained geometric certificates."""
import hashlib, json, resource, sys, time, zipfile
from pathlib import Path
import numpy as np
from full_geometry import project_cylinder
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(ROOT.parent.parent))) / 'cells/geometry'))
import tf2_io

def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()

def sha(p):
    return sha_bytes(Path(p).read_bytes())
start = time.perf_counter()
r5 = json.loads((ROOT / 'RAW_R5.json').read_text())
raw = Path(r5['raw_certificates']['path'])
assert sha(raw) == r5['raw_certificates']['sha256']
p2 = json.loads((ROOT / 'PREREG_R2.json').read_text())
assert sha(p2['input_path']) == p2['input_sha256']
assert sha(ROOT / 'GUIDE_PROFILES.json') == p2['guide_profiles_sha256']
groups = {}
for line in raw.read_text().splitlines():
    row = json.loads(line)
    groups.setdefault(row['case'], []).append(row)
checked_boxes = 0
checked_sites = 0
faults = 0
with zipfile.ZipFile(tf2_io.ZIP) as z:
    for (case, rr) in sorted(groups.items()):
        member = tf2_io.ROOT + '/labelsTr/' + case + '.mha'
        b = z.read(member)
        h = sha_bytes(b)
        (lab, sp, hdr) = tf2_io.read_mha_bytes(b)
        sp = np.array(sp)
        half = sp / 2
        direction = np.fromstring(hdr.get('TransformMatrix', '1 0 0 0 1 0 0 0 1'), sep=' ').reshape(3, 3)
        assert np.max(np.abs(direction.T @ direction - np.eye(3))) < 1e-08
        centers = np.argwhere(np.isin(lab, [3, 4])) * sp
        plan = Path(__import__('os').environ['DENTAL_PROJECT_ROOT']) / 'results/NV1_canals/per_case' / f'{case}.json'
        plan_hash = sha(plan)
        teeth = json.loads(plan.read_text())['teeth']
        for row in rr:
            assert row['source_member'] == member and row['source_member_sha256'] == h
            assert row['plan_file_sha256'] == plan_hash
            E = np.array(row['entry_zyx_mm'])
            a = np.array(row['axis_zyx'])
            L = row['length_mm']
            radius = row['radius_mm']
            tooth = teeth[str(row['fdi'])]
            expected_axis = -np.array(tooth['axis'])
            expected_axis /= np.linalg.norm(expected_axis)
            assert np.linalg.norm(E - np.array(tooth['crest_entry'])) < 1e-12
            assert np.linalg.norm(a - expected_axis) < 1e-12
            assert abs(np.linalg.norm(a) - 1) < 1e-12 and radius == 2.0
            dc = np.linalg.norm(centers - project_cylinder(centers, E, a, L, radius), axis=1)
            initial = float(dc.min())
            halfdiag = float(np.linalg.norm(half))
            assert abs(initial - row['initial_center_upper_mm']) < 1e-09
            assert abs(halfdiag - row['halfdiag_mm']) < 1e-12
            active = np.flatnonzero(dc - halfdiag <= initial + 1e-08)
            stored = np.array([box['voxel_index_in_center_array'] for box in row['records']])
            assert np.array_equal(active, stored)
            assert len(centers) == row['all_voxel_boxes']
            lo = []
            up = []
            for box in row['records']:
                c = centers[box['voxel_index_in_center_array']]
                assert np.linalg.norm(c - np.array(box['center_zyx_mm'])) < 1e-12
                p = np.array(box['witness_zyx_mm'])
                assert np.all(p >= c - half - 1e-10) and np.all(p <= c + half + 1e-10)
                v = p - project_cylinder(p, E, a, L, radius)
                d = float(np.linalg.norm(v))
                normal = v / d if d else np.zeros(3)
                lower = max(0.0, d + float(normal @ (np.where(normal >= 0, c - half, c + half) - p)))
                assert abs(lower - box['lower_mm']) < 1e-08
                assert abs(d - box['upper_mm']) < 1e-08
                lo.append(lower)
                up.append(d)
                checked_boxes += 1
            assert abs(min(lo) - row['lower_mm']) < 1e-08 and abs(min(up) - row['upper_mm']) < 1e-08
            box = row['records'][0]
            c = np.array(box['center_zyx_mm'])
            bad = c + 2 * sp
            assert not (np.all(bad >= c - half) and np.all(bad <= c + half))
            faults += 1
            checked_sites += 1
        del b, lab, centers
result = {'status': 'PASS', 'source_cases_reopened': len(groups), 'source_sites_checked': checked_sites, 'box_membership_and_supporting_planes_checked': checked_boxes, 'pruned_box_coverage_recomputed': True, 'label_and_plan_hashes_verified': True, 'injected_outside_box_witnesses_rejected': faults, 'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'scope': 'Actual source/artifact numerical audit; no independent clinical or anatomical ground truth'}
(ROOT / 'SOURCE_VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
