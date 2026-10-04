from common import *
import zipfile, sys, re, collections
import numpy as np
from report_parser import parse
sys.path.insert(0, str(X76 / 'code'))
import importlib.util
spec = importlib.util.spec_from_file_location('x76_parser', X76 / 'code/report_parser.py')
x76 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(x76)

def reports(ids, modality='intraoral-photo'):
    out = {c: [] for c in ids}
    with zipfile.ZipFile(ZIP) as z:
        for n in sorted(z.namelist()):
            c = n.split('/')[0]
            if c not in out or '/reports_' + modality + '_en/' not in n or (not n.endswith('.txt')):
                continue
            b = z.read(n)
            text = b.decode('utf8')
            p = parse(text)
            it = n.replace('_en/', '_it/')
            ib = z.read(it) if it in z.namelist() else None
            out[c].append(dict(member=n, sha256=digest(b), original_it_member=it if ib else None, original_it_sha256=digest(ib) if ib else None, parsed=p, source_text=text, missing_teeth=x76.parse(text)['observations']['missing_teeth'], modality=modality, resolution='PER_TOOTH', timescale='SIMULTANEOUS', temporal_alignment='Same dataset case, acquisition/report time not supplied'))
    return out

def treatment_rows(rr, geometry):
    rows = []
    ex = []
    for (c, reps) in sorted(rr.items()):
        if not reps:
            ex.append(dict(case_id=c, reason='NO_PHOTO_REPORT'))
            continue
        r = reps[0]
        assertions = r['parsed']['assertions']
        pos = {f: a for a in assertions if a['polarity'] for f in a['teeth']}
        neg = {f: a for a in assertions if not a['polarity'] for f in a['teeth']}
        glob = next((a for a in assertions if a['scope'] == 'ALL_REPORTED_TEETH'), None)
        if glob:
            for (f, g) in geometry.get(c, {}).items():
                if g['label_count'] >= 100:
                    neg[f] = glob
        for f in sorted(set(pos) | set(neg)):
            if c not in geometry or f not in geometry[c]:
                ex.append(dict(case_id=c, fdi=f, reason='NO_GEOMETRY'))
                continue
            a = pos.get(f) or neg[f]
            rows.append(dict(case_id=c, fdi=f, y=int(f in pos), kind=a['kind'], source_member=r['member'], source_sha256=r['sha256'], original_it_member=r['original_it_member'], original_it_sha256=r['original_it_sha256'], clause=a['clause'], scope=a['scope'], resolution='PER_TOOTH', timescale='SIMULTANEOUS'))
        for q in r['parsed']['excluded']:
            ex.append(dict(case_id=c, member=r['member'], **q))
        if not assertions and (not r['parsed']['excluded']):
            ex.append(dict(case_id=c, reason='SILENT_UNKNOWN'))
    return (rows, ex)

def geometry():
    result = {}
    manifest = []
    for p in sorted(X11.glob('F*/labels+landmarks.json')):
        j = read(p)
        c = p.parent.name
        if 'arches' not in j:
            continue
        contactp = X21 / 'raw/cases' / (c + '.json')
        contact = read(contactp) if contactp.exists() else {}
        pairs = contact.get('pairs', [])
        manifest.append(dict(case_id=c, landmarks_path=str(p), landmarks_sha256=sha(p), source_STL_hashes={a: v.get('input_sha256') for (a, v) in j['arches'].items()}, contact_path=str(contactp) if pairs else None, contact_sha256=sha(contactp) if pairs else None))
        teeth = {}
        for (jaw, arch) in j['arches'].items():
            if not arch.get('frame') or not arch.get('landmarks'):
                continue
            frame = arch['frame']
            R = np.array([frame[k] for k in ['right_unit', 'anterior_unit', 'superior_unit']]).T
            center = np.array(frame['center_mm'])
            lm = arch['landmarks']['teeth']
            counts = arch['label_counts']
            N = arch['vertex_count']
            for quadrant in [1, 2] if jaw == 'upper' else [3, 4]:
                for pos in range(1, 9):
                    f = 10 * quadrant + pos
                    t = lm.get(str(f))
                    n = counts.get(str(f), 0)
                    values = dict(label_count=float(n), count_fraction=n / max(N, 1), sample_count=float(t['predicted_sample_points']) if t else 0.0, position=float(pos), upper=float(jaw == 'upper'), right=float(quadrant in [1, 4]))
                    for v in range(1, 9):
                        values['tooth_position_' + str(v)] = float(pos == v)
                    if t:
                        p0 = np.array(t['centroid_xyz_mm'])
                        co = (p0 - center) @ R
                        for (k, v) in zip('ras', co):
                            values['centroid_' + k + '_mm'] = float(v)
                        a = np.array(t['mesial_surface_candidate_xyz_mm'])
                        b = np.array(t['distal_surface_candidate_xyz_mm'])
                        tip = np.array(t['occlusal_tip_candidate_xyz_mm'])
                        values['candidate_span_mm'] = float(np.linalg.norm(a - b))
                        values['tip_center_mm'] = float(np.linalg.norm(tip - p0))
                        values['tip_proximal_mid_mm'] = float(np.linalg.norm(tip - (a + b) / 2))
                        for direction in [-1, 1]:
                            np0 = pos + direction
                            neighbor = lm.get(str(10 * quadrant + np0)) if np0 in range(1, 9) else lm.get(str(21 if quadrant == 1 else 11 if quadrant == 2 else 41 if quadrant == 3 else 31)) if np0 == 0 else None
                            pre = 'mesial' if direction == -1 else 'distal'
                            if neighbor:
                                nc = np.array(neighbor['centroid_xyz_mm'])
                                values[pre + '_center_distance_mm'] = float(np.linalg.norm(nc - p0))
                                op = np.array(neighbor['distal_surface_candidate_xyz_mm' if direction == -1 else 'mesial_surface_candidate_xyz_mm'])
                                my = a if direction == -1 else b
                                values[pre + '_proximal_landmark_distance_mm'] = float(np.linalg.norm(my - op))
                    ps = [v for v in pairs if v[jaw + '_fdi'] == f]
                    values['opposition_degree'] = float(len(ps))
                    if ps:
                        values['opposition_gap_min_mm'] = min((v['minimum_projected_gap_mm'] for v in ps))
                        values['opposition_gap_max_mm'] = max((v['minimum_projected_gap_mm'] for v in ps))
                        values['projected_area_mm2'] = sum((v['projected_proximity_area_mm2'] for v in ps))
                    teeth[f] = values
        if teeth:
            result[c] = teeth
    write('raw/GEOMETRY_INPUT_MANIFEST.json', manifest)
    write('raw/TOOTH_FEATURES.json', {c: {str(f): v for (f, v) in row.items()} for (c, row) in result.items()})
    return result

def load_geometry():
    return {c: {int(f): v for (f, v) in row.items()} for (c, row) in read('raw/TOOTH_FEATURES.json').items()}

def prepare():
    st = time.perf_counter()
    cpu = time.process_time()
    g = geometry()
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    for split in ['train', 'calibration']:
        rr = reports(m[split])
        (rows, ex) = treatment_rows(rr, g)
        write('raw/R1_' + split + '_OBSERVATIONS.json', rr)
        write('raw/R1_' + split + '_ROWS.json', rows)
        write('raw/R1_' + split + '_EXCLUSIONS.json', ex)
    write('raw/PREPARATION_COST.json', cost(st, cpu))
    print('Geometry and train/cal reports prepared', flush=True)
if __name__ == '__main__':
    prepare()
