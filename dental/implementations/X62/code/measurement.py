"""Frozen bench and full-pose readers; no fitting to the submitted measurement."""
import argparse, csv, json, math
from fractions import Fraction
import numpy as np
import mpmath as mp
from common import ROOT, sha, save
from run_round1 import rodrigues

def finite(x):
    x = float(x)
    if not math.isfinite(x):
        raise ValueError('nonfinite measurement')
    return x

def ivfloat(x):
    (n, d) = float(x).as_integer_ratio()
    return mp.iv.mpf(n) / mp.iv.mpf(d)

def pose_enclosure(vertices, center, translation, angles):
    """Interval upper bound for every source vertex under specified finite rigid pose.

    Uses a bounding-box radius, exact binary input floats and interval pi/sine.
    The geodesic triangle inequality bounds the Euler product's rotation angle.
    """
    mp.iv.dps = 40
    bounds = []
    for j in range(3):
        a = ivfloat(float(vertices[:, j].min())) - ivfloat(center[j])
        b = ivfloat(float(vertices[:, j].max())) - ivfloat(center[j])
        hi = max(abs(float(a.a)), abs(float(a.b)), abs(float(b.a)), abs(float(b.b)))
        bounds.append(ivfloat(float(np.nextafter(hi, np.inf))))
    r = mp.iv.sqrt(sum((b * b for b in bounds)))
    tn = mp.iv.sqrt(sum((ivfloat(v) ** 2 for v in translation)))
    deg = sum((Fraction(abs(float(a))) for a in angles))
    if deg >= 180:
        rot = 2 * r
    else:
        angle_iv = mp.iv.mpf(deg.numerator) / mp.iv.mpf(deg.denominator) * mp.iv.pi / 360
        rot = 2 * r * mp.iv.sin(angle_iv)
    bound = tn + rot
    return float(np.nextafter(float(bound.b), np.inf))

def assess_pose(items, rows, meshes):
    if len(items) != 6 or sorted((int(x['FDI']) for x in items)) != [31, 32, 33, 41, 42, 43]:
        raise ValueError('exactly six unique FDI entries required')
    d = {int(x['FDI']): x for x in items}
    out = []
    for (row, m) in zip(rows, meshes):
        x = d[row['FDI']]
        tt = np.array([finite(x[k]) for k in ['mesiodistal_mm', 'occlusogingival_mm', 'buccolingual_mm']])
        aa = np.array([finite(x[k]) for k in ['torque_deg', 'rotation_deg', 'angulation_deg']])
        axes = [np.array(row[k]) for k in ['mesial_axis', 'long_axis', 'labial_axis']]
        t = sum((v * a for (v, a) in zip(tt, axes)))
        R = np.eye(3)
        for (a, angle) in zip(axes, aa):
            R = rodrigues(a, np.deg2rad(angle)) @ R
        c = np.array(row.get('pose_pivot_mm', row['crown_center_mm']))
        changed = (m.vertices - c) @ R.T + c + t
        observed = float(np.max(np.linalg.norm(changed - m.vertices, axis=1)))
        enclosure = pose_enclosure(m.vertices, c, t, aa)
        out.append(dict(FDI=row['FDI'], translation_components_mm=tt.tolist(), angles_deg=aa.tolist(), limits_pass=bool(np.all(np.abs(tt) <= 0.5) and np.all(np.abs(aa) <= 2)), whole_tooth_point_displacement_mm=observed, rigorous_rigid_enclosure_mm=enclosure, enclosure_pass=observed <= enclosure + 1e-12, resolution='PER_POINT', pose_pivot_mm=c.tolist(), pivot_origin=row.get('pivot_origin', 'LEGACY_VIRTUAL_CROWN_CENTER'), interpretation='rigid geometric implication; no achieved tooth movement predicted'))
    return dict(all_six_pass=all((x['limits_pass'] for x in out)), per_tooth=out, schema='full-pose-metrology-v1', measurement_kind='submitted relative errors; provenance must identify measured versus virtual', timescale='SIMULTANEOUS', registration='relative to nominal target in tooth frame; common global frame registration must be provided')

def target_rows(rows, targets):
    if len(targets) != 6 or sorted((int(x['FDI']) for x in targets)) != [31, 32, 33, 41, 42, 43]:
        raise ValueError('six unique frozen target pivots required')
    td = {int(x['FDI']): x for x in targets}
    out = []
    for row in rows:
        item = td[row['FDI']]
        if item.get('unit') != 'mm':
            raise ValueError('target pivot unit must be mm')
        xyz = item['pivot_mm']
        if len(xyz) != 3:
            raise ValueError('target pivot needs three coordinates')
        out.append(row | dict(pose_pivot_mm=[finite(x) for x in xyz], pivot_origin=item['origin']))
    return out

def read_targets(path):
    receipt = json.loads((ROOT / 'raw/FROZEN_TARGET_RECEIPT_R4.json').read_text())
    if sha(path) != receipt['sha256']:
        raise ValueError('frozen target hash mismatch; new targets need a new prediction round')
    return json.loads(open(path).read())['targets']

def assess_force(items, predicted_hash):
    path = ROOT / 'FROZEN_PREDICTIONS.json'
    if sha(path) != predicted_hash:
        raise ValueError('frozen prediction hash mismatch')
    expected = json.loads(path.read_text())['force_mode_predictions']
    indexed = {(int(x['mode']), int(x['FDI'])): x for x in items}
    if len(items) != 72 or len(indexed) != 72:
        raise ValueError('12 modes x 6 teeth required, without duplicates')
    out = []
    for mode in expected:
        for (fdi, f) in zip([33, 32, 31, 41, 42, 43], mode['force_N']):
            key = (mode['mode'], fdi)
            if key not in indexed:
                raise ValueError('mode/FDI missing')
            row = indexed[key]
            measured = [finite(row[k]) for k in ['Fx_N', 'Fy_N', 'Fz_N']]
            errors = np.abs(np.array(measured) - np.array(f))
            tol = np.maximum(0.02, 0.15 * np.maximum(np.abs(f), 0.1))
            out.append(dict(mode=key[0], FDI=fdi, error_N=errors.tolist(), tolerance_N=tol.tolist(), pass_gate=bool(np.all(errors <= tol))))
    return dict(pass_gate=all((x['pass_gate'] for x in out)), cells=out, prediction_sha256=predicted_hash, refit=False, physical_model_scope='conditional small-displacement rod; biological response excluded')

def main():
    from check_inputs import check
    check()
    a = argparse.ArgumentParser()
    a.add_argument('--force-csv')
    a.add_argument('--pose-csv')
    a.add_argument('--predictions-sha256')
    a.add_argument('--targets')
    a.add_argument('--out', required=True)
    args = a.parse_args()
    if bool(args.force_csv) == bool(args.pose_csv):
        raise ValueError('choose one measured file')
    source = args.force_csv or args.pose_csv
    with open(source) as f:
        items = list(csv.DictReader(f))
    if args.force_csv:
        r = assess_force(items, args.predictions_sha256)
    else:
        if not args.targets:
            raise ValueError('physical pose reader requires explicit frozen target pivots; no crown-center substitute')
        from patient import load
        (_, _, rows, meshes) = load()
        rows = target_rows(rows, read_targets(args.targets))
        r = assess_pose(items, rows, meshes)
        r['target_file'] = args.targets
        r['target_sha256'] = sha(args.targets)
    r.update(input_file=source, input_sha256=sha(source), measurement_origin='USER_SUPPLIED; independent physical provenance not automatically verified')
    save(args.out, r)
    print(json.dumps({'pass_gate': r.get('pass_gate', r.get('all_six_pass')), 'out': args.out}))
if __name__ == '__main__':
    main()
