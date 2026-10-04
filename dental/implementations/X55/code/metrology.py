"""Regional rigid scan/reference metrology. mm throughout; no scale fitting.

Normal projection is signed relative to the reference outward face normal.
It is neither an inside/outside test nor a seated crown/die gap.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from scipy.stats import t as student_t
REGIONS = ('marginal', 'intaglio', 'occlusal', 'axial')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def load_mesh(path, units):
    if units != 'mm':
        raise ValueError('Explicit millimetres required; convert and freeze inputs externally')
    m = trimesh.load_mesh(path, process=False)
    if not isinstance(m, trimesh.Trimesh) or not len(m.faces):
        raise ValueError('A nonempty triangular surface is required')
    if not np.isfinite(m.vertices).all() or np.any(m.area_faces <= 1e-14):
        raise ValueError('Nonfinite vertices or degenerate faces')
    return m

def load_regions(path, design_path, mesh):
    r = json.loads(Path(path).read_text())
    if r.get('design_sha256') != sha(design_path):
        raise ValueError('Region map does not match exact design hash')
    labels = np.asarray(r['face_regions'])
    if labels.shape != (len(mesh.faces),) or set(labels) != set(REGIONS):
        raise ValueError('Exactly four nonempty regions, one label per design face, required')
    return labels

def sample_surface(mesh, count, seed, labels=None):
    """Uniform area samples within region. Area weight preserves surface measure."""
    rng = np.random.default_rng(seed)
    if labels is None:
        labels = np.full(len(mesh.faces), 'all')
    groups = sorted(set(labels))
    (points, normals, ids, regions, weights) = ([], [], [], [], [])
    for (k, group) in enumerate(groups):
        faces = np.flatnonzero(labels == group)
        n = count // len(groups) + (k < count % len(groups))
        if n < 6:
            raise ValueError('At least six samples per region required')
        areas = mesh.area_faces[faces]
        f = rng.choice(faces, n, p=areas / areas.sum())
        uv = rng.random((n, 2))
        root = np.sqrt(uv[:, 0])
        bary = np.c_[1 - root, root * (1 - uv[:, 1]), root * uv[:, 1]]
        points.append(np.einsum('ij,ijk->ik', bary, mesh.triangles[f]))
        normals.append(mesh.face_normals[f])
        ids.append(f)
        regions.extend([group] * n)
        weights.extend([areas.sum() / n] * n)
    return (np.vstack(points), np.vstack(normals), np.concatenate(ids), np.asarray(regions), np.asarray(weights))

def closest(mesh, points):
    answers = [trimesh.proximity.closest_point(mesh, points[i:i + 512]) for i in range(0, len(points), 512)]
    return tuple((np.concatenate([a[k] for a in answers]) for k in range(3)))

def transform(points, T):
    return points @ T[:3, :3].T + T[:3, 3]

def kabsch(source, target):
    (a, b) = (np.asarray(source, float), np.asarray(target, float))
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3 or (len(a) < 3):
        raise ValueError('At least three matched 3D datum points required')
    (ca, cb) = (a.mean(0), b.mean(0))
    if np.linalg.matrix_rank(a - ca, tol=1e-10) < 2:
        raise ValueError('Collinear datum cannot fix rotation')
    (U, _, Vt) = np.linalg.svd((a - ca).T @ (b - cb))
    R = Vt.T @ np.diag([1, 1, np.linalg.det(Vt.T @ U.T)]) @ U.T
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = cb - R @ ca
    return T

def validate_transform(T):
    T = np.asarray(T, float)
    if T.shape != (4, 4) or not np.isfinite(T).all() or (not np.allclose(T[3], [0, 0, 0, 1], atol=1e-12)) or (not np.allclose(T[:3, :3].T @ T[:3, :3], np.eye(3), atol=1e-08)) or (abs(np.linalg.det(T[:3, :3]) - 1) > 1e-08):
        raise ValueError('A proper rigid homogeneous transform is required')
    return T

def register(design, scansites, labels, mode, initial=None, iterations=60, huber_mm=0.02, max_distance_mm=0.8):
    """Local robust point-to-plane ICP; region-balanced best fit or hard lock.

    Rank is an identifiability diagnostic, not a covariance/error guarantee.
    The caller must supply rough alignment; there is no global basin certificate.
    """
    T = np.eye(4) if initial is None else validate_transform(initial).copy()
    if initial is None:
        T[:3, 3] = design.centroid - scansites.mean(0)
    history = []
    converged = False
    for it in range(iterations):
        q = transform(scansites, T)
        (p, dist, f) = closest(design, q)
        n = design.face_normals[f]
        region = labels[f]
        active = dist < max_distance_mm
        if mode in ('marginal', 'intaglio'):
            active &= region == mode
        if active.sum() < 12:
            raise ValueError('Too little overlap for selected registration region')
        residual = np.einsum('ij,ij->i', q - p, n)
        weight = np.zeros(len(q))
        for group in np.unique(region[active]):
            mask = active & (region == group)
            weight[mask] = 1 / mask.sum()
        weight *= np.minimum(1, huber_mm / np.maximum(abs(residual), 1e-15))
        J = np.c_[np.cross(q, n), n]
        A = J[active] * np.sqrt(weight[active, None])
        b = -residual[active] * np.sqrt(weight[active])
        (delta, _, rank, singular) = np.linalg.lstsq(A, b, rcond=1e-10)
        factor = min(1.0, 0.15 / max(np.linalg.norm(delta[:3]), 1e-15), 0.5 / max(np.linalg.norm(delta[3:]), 1e-15))
        delta *= factor
        inc = np.eye(4)
        inc[:3, :3] = Rotation.from_rotvec(delta[:3]).as_matrix()
        inc[:3, 3] = delta[3:]
        T = inc @ T
        history.append(dict(iteration=it, active=int(active.sum()), rank=int(rank), weighted_rms_um=float(np.sqrt(np.average(residual[active] ** 2, weights=weight[active])) * 1000), step=float(np.linalg.norm(delta))))
        if np.linalg.norm(delta) < 1e-07:
            converged = True
            break
    return (T, dict(converged=converged, iterations=len(history), history=history, observability_rank=history[-1]['rank'], rigorous_pose_uncertainty_enclosure='MISSING', global_registration_certificate='MISSING'))

def summarize(fields_mm, distances_mm, regions, weights, coverage_mm):
    """Finite reference-site descriptors; CI counts repeats, never spatial points."""
    reports = {}
    for region in ('all', *REGIONS):
        chosen = np.ones(len(regions), bool) if region == 'all' else regions == region
        covered = np.all(distances_mm <= coverage_mm, axis=0) & chosen
        n = int(covered.sum())
        area_fraction = float(weights[covered].sum() / weights[chosen].sum())
        out = dict(resolution='PER_TOOTH' if region == 'all' else 'PER_SURFACE_REGION', sites=int(chosen.sum()), retained_sites=n, retained_area_fraction=area_fraction, rejected_fraction=1 - area_fraction, rejected_reason='reference site farther than coverage gate in at least one scan', status='AVAILABLE' if n >= 6 and area_fraction >= 0.95 else 'INSUFFICIENT_COVERAGE')
        if n < 6:
            reports[region] = out
            continue
        d = fields_mm[:, covered] * 1000
        w = weights[covered]
        w = w / w.sum()
        mean = d.mean(0)
        per_scan_bias = d @ w
        out.update(signed_bias_um=float(mean @ w), mean_field_rms_um=float(np.sqrt(mean ** 2 @ w)), per_scan_normal_rms_um=np.sqrt(d ** 2 @ w).tolist(), pooled_normal_rms_um=float(np.sqrt(np.mean(d ** 2 @ w))), pooled_euclidean_rms_um=float(np.sqrt(np.mean((distances_mm[:, covered] * 1000) ** 2 @ w))), abs_mean_field_p95_um=float(np.percentile(abs(mean), 95)), observed_abs_max_um=float(abs(d).max()))
        if len(d) >= 2:
            s2 = d.var(0, ddof=1)
            sr = float(np.sqrt(s2 @ w))
            pair2 = [(d[i] - d[j]) ** 2 @ w for i in range(len(d)) for j in range(i)]
            half = float(student_t.ppf(0.975, len(d) - 1) * per_scan_bias.std(ddof=1) / np.sqrt(len(d)))
            out.update(repeatability_normal_sd_um=sr, mean_pairwise_rms_um=float(np.mean(np.sqrt(pair2))), rms_pairwise_difference_um=float(np.sqrt(np.mean(pair2))), pairwise_identity_error_um=abs(float(np.sqrt(np.mean(pair2) / 2)) - sr), signed_bias_repeat_only_ci95_um=[float(mean @ w - half), float(mean @ w + half)], ci_assumption='Independent Gaussian repeated scan regional biases; systematic reference error excluded')
        else:
            out['precision_status'] = 'UNKNOWN_SINGLE_SCAN'
        reports[region] = out
    return reports

def analyze(design_path, scan_paths, regions_path, out, mode='balanced', units='mm', samples=2400, seed=5501, metadata=None, coverage_mm=0.3):
    start = time.perf_counter()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    design = load_mesh(design_path, units)
    labels = load_regions(regions_path, design_path, design)
    meta = {} if metadata is None else json.loads(Path(metadata).read_text())
    hashes = [sha(p) for p in scan_paths]
    if len(set(hashes)) != len(hashes):
        raise ValueError('Duplicate scan file is not an independent repeat')
    if len(scan_paths) > 1:
        rows = meta.get('scans', {})
        ids = [rows.get(Path(p).name, {}).get('specimen_id') for p in scan_paths]
        if any((v is None for v in ids)) or len(set(ids)) != 1:
            raise ValueError('Repeatability requires metadata identifying the same physical specimen')
        if meta.get('repeat_conditions') != 'same_scanner_operator_protocol':
            raise ValueError('Declare repeated measurement conditions explicitly')
    reference_kind = meta.get('reference_kind', 'design')
    specimen_ids = [meta.get('scans', {}).get(Path(p).name, {}).get('specimen_id', 'UNKNOWN') for p in scan_paths]
    if reference_kind not in ('design', 'independent_measurement'):
        raise ValueError('Unknown reference kind')
    if reference_kind == 'independent_measurement' and (not meta.get('reference_locator')):
        raise ValueError('Independent reference requires locator')
    (points, normals, faces, regions, area) = sample_surface(design, samples, seed, labels)
    (fields, distances, registration) = ([], [], [])
    for path in scan_paths:
        tick = time.perf_counter()
        scan = load_mesh(path, units)
        ratio = np.linalg.norm(scan.extents) / np.linalg.norm(design.extents)
        if not 0.5 < ratio < 2:
            raise ValueError('Gross size mismatch: units, wrong part or incomplete scan')
        row = meta.get('scans', {}).get(Path(path).name, {})
        if mode == 'datum':
            if 'scan_datum_points_mm' in row and 'reference_datum_points_mm' in row:
                T = kabsch(row['scan_datum_points_mm'], row['reference_datum_points_mm'])
                datum_residual = transform(np.asarray(row['scan_datum_points_mm']), T) - np.asarray(row['reference_datum_points_mm'])
            elif 'scan_to_reference' in row:
                T = validate_transform(row['scan_to_reference'])
                datum_residual = None
            else:
                raise ValueError('Datum mode requires frozen scan_to_reference transform')
            diagnostic = dict(converged=True, rigorous_pose_uncertainty_enclosure='MISSING_UNLESS_EXTERNAL_DATUM_BOUNDS_SUPPLIED', datum_rms_um=float(np.sqrt(np.mean(datum_residual ** 2)) * 1000) if datum_residual is not None else None, source='external frozen transform; no surface fit')
        else:
            sites = sample_surface(scan, samples, seed + 1)[0]
            (T, diagnostic) = register(design, sites, labels, mode, row.get('initial_scan_to_reference'))
        scan.apply_transform(T)
        (q, dist, _) = closest(scan, points)
        d = np.einsum('ij,ij->i', q - points, normals)
        fields.append(d)
        distances.append(dist)
        registration.append(dict(file=str(path), sha256=sha(path), scan_to_reference=T.tolist(), fit_and_measure_seconds=time.perf_counter() - tick, **diagnostic))
    (fields, distances) = (np.asarray(fields), np.asarray(distances))
    tick = time.perf_counter()
    reports = summarize(fields, distances, regions, area, coverage_mm)
    np.savez_compressed(out / 'point_fields.npz', xyz_mm=points, normals=normals, design_face=faces, regions=regions, area_weights_mm2=area, specimen_id=np.asarray(specimen_ids[0]), normal_deviation_mm=fields, euclidean_distance_mm=distances)
    result = dict(schema='K48-regional-metrology-v1', claim_type='capability', units='um except geometry and transforms in mm', registration_mode=mode, design_sha256=sha(design_path), regions_sha256=sha(regions_path), reference_kind=reference_kind, reference_locator=meta.get('reference_locator'), specimen_id=specimen_ids[0], estimand='manufactured_design_discrepancy' if reference_kind == 'design' else 'measurement_reference_bias', scanner_trueness_status='UNKNOWN_CAD_IS_NOT_ACTUAL_PART_REFERENCE' if reference_kind == 'design' else 'CONDITIONAL_ON_REFERENCE_AND_DATUM_UNCERTAINTY', scanner_systematic_uncertainty_um=meta.get('reference_systematic_uncertainty_um', 'UNKNOWN'), standard_scope='ISO 5725 terminology; ISO 12836 mounted-digitizer scope; no conformity certification', missing_uncertainty=['reference bias unless independently supplied', 'pose enclosure', 'tessellation', 'finite surface sampling maxima'], registration=registration, regions=reports, field_sha256=sha(out / 'point_fields.npz'), measurement_status='COMPUTED_CONDITIONAL' if all((r['converged'] for r in registration)) else 'UNKNOWN_REGISTRATION_NOT_CONVERGED', query_and_export_seconds=time.perf_counter() - tick, total_seconds=time.perf_counter() - start)
    write_json(out / 'report.json', result)
    return result

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design', required=True)
    p.add_argument('--scans', nargs='+', required=True)
    p.add_argument('--regions', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--units', required=True, choices=['mm'])
    p.add_argument('--metadata')
    p.add_argument('--registration', choices=['balanced', 'marginal', 'intaglio', 'datum'], default='balanced')
    p.add_argument('--samples', type=int, default=2400)
    p.add_argument('--coverage-mm', type=float, default=0.3)
    a = p.parse_args()
    if a.coverage_mm <= 0:
        p.error('positive coverage distance required')
    try:
        r = analyze(a.design, a.scans, a.regions, a.output, a.registration, a.units, a.samples, metadata=a.metadata, coverage_mm=a.coverage_mm)
    except (ValueError, KeyError) as e:
        p.exit(2, 'METROLOGY_REJECT: ' + str(e) + '\n')
    print(json.dumps(dict(output=a.output, estimand=r['estimand'], seconds=r['total_seconds'])))
if __name__ == '__main__':
    main()
