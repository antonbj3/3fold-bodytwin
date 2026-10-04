"""Actual isolated missing sites; no virtual healing. Coordinates: array (z,y,x), mm.

All clinical interpretations remain UNKNOWN. SDF supports the finite-body successor.
"""
from dental_release.paths import expand as _release_expand
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import zipfile
import numpy as np
from scipy import ndimage as ndi
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'cells/geometry'))
import tf2_io
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X20_short_implant_sinus'))
POSTERIOR = (14, 15, 16, 17, 24, 25, 26, 27)
STEP = 0.075

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def unit(a):
    a = np.asarray(a, float)
    return a / np.linalg.norm(a)

def check_frozen(name):
    assert sha(HERE / f'{name}.json') == (HERE / f'{name}.sha256').read_text().strip(), name

def select():
    cases = [json.loads(s) for s in (ROOT / 'results/BONE_HU/tf2_census.jsonl').read_text().splitlines()]
    sites = []
    for rec in cases:
        q = rec['counts']
        for t in POSTERIOR:
            if str(t) not in q and min(int(q.get(str(t - 1), 0)), int(q.get(str(t + 1), 0))) >= 500 and (sum((int(q.get(str(i), 0)) for i in (5, 6))) >= 500):
                sites.append({'case': rec['case'], 'fdi': t, 'stratum': 'PM' if t % 10 < 6 else 'M'})
    return (cases, sites)

def frame(case, fdi):
    p = ROOT / f'results/NV1_canals/per_case/{case}.json'
    j = json.loads(p.read_text())
    (l, r) = (j['teeth'][str(fdi - 1)], j['teeth'][str(fdi + 1)])
    anchor = (np.asarray(l['crest_entry']) + r['crest_entry']) / 2
    axis = -unit(unit(l['axis']) + unit(r['axis']))
    tangent = np.asarray(r['centroid']) - l['centroid']
    tangent = unit(tangent - np.dot(tangent, axis) * axis)
    normal = unit(np.cross(axis, tangent))
    return (anchor, axis, normal, tangent, sha(p))

class Local:

    def __init__(self, labels, spacing, anchor, axis):
        self.sp = np.asarray(spacing)
        ends = np.array([anchor - 14 * axis, anchor + 32 * axis])
        self.lo = np.maximum(np.floor((ends.min(0) - 15) / self.sp).astype(int), 0)
        self.hi = np.minimum(np.ceil((ends.max(0) + 15) / self.sp).astype(int) + 1, labels.shape)
        sl = tuple((slice(a, b) for (a, b) in zip(self.lo, self.hi)))
        self.lab = np.array(labels[sl], copy=True)
        self.occ = {k: v.astype(np.float32) for (k, v) in {'bone': self.lab == 2, 'sinus': np.isin(self.lab, (5, 6)), 'restoration': np.isin(self.lab, (8, 9, 10)), 'neighbor': (self.lab >= 11) & (self.lab <= 28)}.items()}
        self.sdf = {}

    def sample(self, key, points, order=1, fill=0):
        q = (np.asarray(points) / self.sp - self.lo).reshape(-1, 3).T
        arr = self.occ.get(key, self.sdf.get(key))
        return ndi.map_coordinates(arr, q, order=order, mode='constant', cval=fill).reshape(np.asarray(points).shape[:-1])

    def signed(self, key):
        m = self.occ[key].astype(bool)
        if not m.any():
            return None
        half = self.sp.min() / 2
        self.sdf[key + '_sdf'] = np.where(m, ndi.distance_transform_edt(m, sampling=self.sp) - half, -(ndi.distance_transform_edt(~m, sampling=self.sp) - half)).astype(np.float32)
        return self.sdf[key + '_sdf']

    def distance(self, key):
        m = self.occ[key].astype(bool)
        if not m.any():
            return None
        self.sdf[key + '_dist'] = (ndi.distance_transform_edt(~m, sampling=self.sp) - self.sp.min() / 2).astype(np.float32)
        return self.sdf[key + '_dist']

    def save(self, path, anchor, axis, normal, tangent):
        np.savez_compressed(path, lab=self.lab, lo=self.lo, spacing=self.sp, anchor=anchor, axis=axis, normal=normal, tangent=tangent)

def segments(v):
    q = np.diff(np.r_[False, v, False].astype(int))
    return list(zip(np.where(q == 1)[0], np.where(q == -1)[0]))

def boundary(values, coords, i, entry=True):
    if entry:
        (a, b) = (i - 1, i)
    else:
        (a, b) = (i - 1, i)
    if a < 0 or b >= len(values):
        return None
    den = values[b] - values[a]
    if abs(den) < 1e-12:
        return float((coords[a] + coords[b]) / 2)
    return float(coords[a] + (0.5 - values[a]) * (coords[b] - coords[a]) / den)

def trace(local, anchor, axis, normal, order=1):
    s = np.arange(-12, 30 + STEP / 2, STEP)
    pts = anchor + s[:, None] * axis
    b = local.sample('bone', pts, order=order)
    sin = local.sample('sinus', pts, order=order)
    runs = [(a, e) for (a, e) in segments(b >= 0.5) if (e - a) * STEP >= 1.0]
    if not runs:
        return {'valid': False, 'reason': 'no_contiguous_bone_crest'}
    (a, e) = runs[0]
    crest = boundary(b, s, a)
    bone_end = boundary(b, s, e, False)
    if crest is None or bone_end is None:
        return {'valid': False, 'reason': 'truncated_bone_trace'}
    ix = np.where((s > crest) & (sin >= 0.5))[0]
    if not len(ix):
        return {'valid': False, 'reason': 'no_sinus_within_30mm', 'crest_s_mm': crest, 'bone_height_mm': bone_end - crest}
    sinus_start = boundary(sin, s, ix[0])
    if sinus_start is None:
        return {'valid': False, 'reason': 'truncated_sinus_trace'}
    widths = {}
    r = np.arange(-12, 12 + STEP / 2, STEP)
    for depth in (1, 3, 7):
        p = anchor + (crest + depth) * axis + r[:, None] * normal
        w = local.sample('bone', p, order=order)
        selected = [(lo, hi) for (lo, hi) in segments(w >= 0.5) if r[lo] <= 0 <= r[hi - 1]]
        if not selected:
            widths[str(depth)] = None
            continue
        (lo, hi) = selected[0]
        (wlo, whi) = (boundary(w, r, lo), boundary(w, r, hi, False))
        widths[str(depth)] = None if wlo is None or whi is None else whi - wlo
    gap = sinus_start - bone_end
    out = {'valid': gap <= 0.9 and gap >= -0.3 and (widths['1'] is not None), 'reason': 'ok' if gap <= 0.9 and gap >= -0.3 and (widths['1'] is not None) else 'bone_to_sinus_gap_or_no_width', 'crest_s_mm': crest, 'crest_mm': (anchor + crest * axis).tolist(), 'bone_height_mm': bone_end - crest, 'crest_to_sinus_label_mm': sinus_start - crest, 'bone_to_sinus_gap_mm': gap, 'widths_mm': widths}
    return out

def restoration_near(local, anchor, axis):
    q = np.argwhere(local.occ['restoration'] > 0)
    if not len(q):
        return (False, 0)
    delta = (q + local.lo) * local.sp - anchor
    axial = delta @ axis
    radial = np.linalg.norm(delta - axial[:, None] * axis, axis=1)
    hit = (radial <= 3) & (axial >= -6) & (axial <= 8)
    return (bool(hit.any()), int(hit.sum()))

def scalar_decide(rec, half=0.6):
    if not rec['valid']:
        return 'uncertain'
    (h, w) = (rec['bone_height_mm'], rec['widths_mm']['1'])
    if h - half >= 7 and w - half >= 6:
        return 'short_geometry'
    if h + half < 7 and w - half >= 6:
        return 'height_augmentation_geometry'
    return 'uncertain'

def write_csv(rows, path):
    cols = ['case', 'fdi', 'stratum', 'valid', 'reason', 'bone_height_mm', 'crest_to_sinus_label_mm', 'bone_to_sinus_gap_mm', 'width1_mm', 'width3_mm', 'width7_mm', 'rule5_flag', 'scalar_nominal', 'scalar_interval', 'shortest_scalar_catalog_length_mm', 'clinical_decision']
    with open(path, 'w') as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        for row in rows:
            r = {k: row.get(k) for k in cols}
            for d in (1, 3, 7):
                r[f'width{d}_mm'] = row.get('widths_mm', {}).get(str(d))
            w.writerow(r)

def run_r1():
    start = time.monotonic()
    check_frozen('PREREG_R1')
    check_frozen('FROZEN_PREDICTIONS')
    DATA.mkdir(parents=True, exist_ok=True)
    assert sum((p.stat().st_size for p in DATA.glob('*') if p.is_file())) < 3000000000
    (cases, sites) = select()
    dump(HERE / 'SELECTED_SITES.json', {'census_n': len(cases), 'sinus_labeled_n': sum((sum((int(c['counts'].get(str(i), 0)) for i in (5, 6))) >= 500 for c in cases)), 'selected': sites})
    manifest = {'census': {'path': str(ROOT / 'results/BONE_HU/tf2_census.jsonl'), 'sha256': sha(ROOT / 'results/BONE_HU/tf2_census.jsonl')}, 'archive': {'path': tf2_io.ZIP, 'bytes': Path(tf2_io.ZIP).stat().st_size}, 'members': [], 'local_crops': []}
    rows = []
    with zipfile.ZipFile(tf2_io.ZIP) as z:
        for case in sorted(set((s['case'] for s in sites))):
            member = f'{tf2_io.ROOT}/labelsTr/{case}.mha'
            raw = z.read(member)
            (lab, spacing, header) = tf2_io.read_mha_bytes(raw)
            manifest['members'].append({'name': member, 'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw), 'zip_crc32': z.getinfo(member).CRC, 'spacing': spacing, 'transform': header.get('TransformMatrix'), 'offset': header.get('Offset')})
            for s in [q for q in sites if q['case'] == case]:
                rec = dict(s, clinical_decision='UNKNOWN')
                try:
                    (anchor, axis, normal, tangent, frame_hash) = frame(case, s['fdi'])
                    f = Local(lab, spacing, anchor, axis)
                    rec.update(trace(f, anchor, axis, normal))
                    rec['frame_sha256'] = frame_hash
                    rec['axis'] = axis.tolist()
                    rec['anchor'] = anchor.tolist()
                    (near, n) = restoration_near(f, anchor, axis)
                    rec['restoration_voxels_near'] = n
                    if near:
                        (rec['valid'], rec['reason']) = (False, 'restoration_at_proposed_gap')
                    h = rec.get('bone_height_mm')
                    rec['rule5_flag'] = 'lift_trigger' if h is not None and h < 5 else 'no_lift_trigger' if h is not None else 'unknown'
                    rec['scalar_nominal'] = scalar_decide(rec, half=0)
                    rec['scalar_interval'] = scalar_decide(rec)
                    rec['height_interval_mm'] = [max(0, h - 0.6), h + 0.6] if h is not None else None
                    w = rec.get('widths_mm', {}).get('1')
                    rec['width1_interval_mm'] = [max(0, w - 0.6), w + 0.6] if w is not None else None
                    rec['shortest_scalar_catalog_length_mm'] = 6 if rec['scalar_interval'] == 'short_geometry' else None
                    check = trace(f, anchor, axis, normal, order=0)
                    errors = {}
                    for k in ('bone_height_mm', 'crest_to_sinus_label_mm'):
                        if k in rec and k in check:
                            errors[k] = abs(rec[k] - check[k])
                    for d in ('1', '3', '7'):
                        if rec.get('widths_mm', {}).get(d) is not None and check.get('widths_mm', {}).get(d) is not None:
                            errors['width' + d] = abs(rec['widths_mm'][d] - check['widths_mm'][d])
                    rec['independent_voxel_read'] = check
                    rec['cross_read_errors_mm'] = errors
                    rec['cross_read_pass'] = all((v <= 0.15 + 1e-08 for v in errors.values())) and bool(errors)
                    path = DATA / f"{case}_{s['fdi']}.npz"
                    f.save(path, anchor, axis, normal, tangent)
                    manifest['local_crops'].append({'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size})
                except (KeyError, ValueError) as ex:
                    rec.update(valid=False, reason=f'frame_failure:{ex}', scalar_interval='uncertain')
                rows.append(rec)
                print(case, s['fdi'], rec['reason'], rec.get('bone_height_mm'), rec.get('widths_mm', {}).get('1'), rec['scalar_interval'], flush=True)
            del lab, raw
    dump(HERE / 'DATA_MANIFEST.json', manifest)
    with open(HERE / 'RAW_R1.jsonl', 'w') as f:
        for r in rows:
            f.write(json.dumps(r, allow_nan=False) + '\n')
    write_csv(rows, HERE / 'MEASUREMENTS_R1.csv')
    valid = [r for r in rows if r['valid']]
    counts = {k: sum((r.get('scalar_interval') == k for r in rows)) for k in ('short_geometry', 'height_augmentation_geometry', 'uncertain')}
    disagreements = [r for r in valid if r['scalar_interval'] == 'short_geometry' and r['rule5_flag'] == 'lift_trigger' or (r['scalar_interval'] == 'height_augmentation_geometry' and r['rule5_flag'] == 'no_lift_trigger')]
    summary = {'n_selected': len(rows), 'n_valid': len(valid), 'n_cases': len(set((r['case'] for r in rows))), 'counts': counts, 'coverage_gate': 'PASS' if len(valid) >= 8 else 'FAIL', 'definite_rule5_disagreements': len(disagreements), 'disagreement_denominator': len(valid), 'disagreement_ids': [f"{r['case']}:{r['fdi']}" for r in disagreements], 'measurement_abstentions': sum((r.get('scalar_interval') == 'uncertain' for r in valid)), 'semantic_or_selection_exclusions': len(rows) - len(valid), 'external_distribution_gate': 'UNKNOWN_MINIMUM_N', 'cross_read_pass_n': sum((r.get('cross_read_pass', False) for r in rows)), 'elapsed_s': time.monotonic() - start, 'data_bytes': sum((p.stat().st_size for p in DATA.glob('*.npz'))), 'external_referent': json.loads((HERE / 'EXTERNAL_REFERENTS.json').read_text())['anatomy']}
    dump(HERE / 'SUMMARY_R1.json', summary)
    dump(HERE / 'CURRENT_WORK_STATE.json', {'lane': 'X20-short-implant-sinus', 'milestone': 'R1_completed', 'latest_gate': summary['coverage_gate'], 'ongoing': 'Finite-body successor preregistration', 'next_operation': 'Replace scalar height/width by diameter-dependent all-surface field margins', 'summary': summary})
    print(json.dumps(summary, indent=2))
if __name__ == '__main__':
    run_r1()
