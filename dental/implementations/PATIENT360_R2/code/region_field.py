"""Native source-addressed region fields, not a fabricated common SDF.

IOS rays have signed vertical gaps, no volumetric inside predicate. CBCT has
piecewise-constant voxel ownership and a separately named grid-centre EDT.
Missing modalities, out-of-domain queries and ambiguous owners stay unknown.
"""
from dataclasses import dataclass
import copy
import numpy as np
from scipy.spatial import cKDTree
from common import *

def transform(points, affine):
    points = np.asarray(points, float)
    A = np.asarray(affine, float)
    return points @ A[:3, :3].T + A[:3, 3]

def affine_valid(A, rigid=False):
    A = np.asarray(A, float)
    if A.shape != (4, 4) or not np.isfinite(A).all() or (not np.array_equal(A[3], [0, 0, 0, 1])):
        raise ValueError('AFFINE_INVALID')
    if abs(np.linalg.det(A[:3, :3])) < 1e-12:
        raise ValueError('SINGULAR_FRAME')
    if rigid and (not np.allclose(A[:3, :3].T @ A[:3, :3], np.eye(3), rtol=0, atol=1e-10) or np.linalg.det(A[:3, :3]) < 0):
        raise ValueError('NONRIGID_REGISTRATION')
    return A

def read_stl(blob):
    n = int.from_bytes(blob[80:84], 'little')
    if len(blob) != 84 + 50 * n:
        raise ValueError('BINARY_STL_LENGTH')
    dt = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attribute', '<u2')])
    return np.frombuffer(blob, dt, offset=84, count=n)['vertices'].astype(float)

def ray_height(tri, xy, upper):
    """Plane solve + AABB screen; independent of X18's radius/barycentric code."""
    xy = np.asarray(xy, float)
    lo = tri[:, :, :2].min(1)
    hi = tri[:, :, :2].max(1)
    e1 = tri[:, 1] - tri[:, 0]
    e2 = tri[:, 2] - tri[:, 0]
    norm = np.cross(e1, e2)
    valid = np.abs(norm[:, 2]) > 1e-12
    height = np.full(len(xy), np.nan)
    faces = np.full(len(xy), -1, int)
    bary = np.full((len(xy), 3), np.nan)
    for (k, p) in enumerate(xy):
        ids = np.flatnonzero(valid & np.all(lo <= p + 1e-12, axis=1) & np.all(hi >= p - 1e-12, axis=1))
        if not len(ids):
            continue
        M = np.stack([e1[ids, :2], e2[ids, :2]], axis=2)
        uv = np.linalg.solve(M, (p - tri[ids, 0, :2])[..., None])[..., 0]
        ok = (uv[:, 0] >= -1e-10) & (uv[:, 1] >= -1e-10) & (uv.sum(1) <= 1 + 1e-10)
        ids = ids[ok]
        uv = uv[ok]
        if not len(ids):
            continue
        h = tri[ids, 0, 2] - (norm[ids, 0] * (p[0] - tri[ids, 0, 0]) + norm[ids, 1] * (p[1] - tri[ids, 0, 1])) / norm[ids, 2]
        j = np.argmin(h) if upper else np.argmax(h)
        height[k] = h[j]
        faces[k] = ids[j]
        bary[k] = [1 - uv[j].sum(), *uv[j]]
    return (height, faces, bary)

class RegionField:

    def __init__(self, descriptor, arrays):
        self.meta = descriptor
        self.arrays = arrays
        if descriptor.get('unit') != 'mm':
            raise ValueError('FIELD_UNIT_NOT_MM')
        self.A = affine_valid(descriptor['native_to_patient'])
        if not descriptor.get('patient_id') or not descriptor.get('frame_id') or (not descriptor.get('source_sha256')):
            raise ValueError('FIELD_PROVENANCE_MISSING')
        self.inv = np.linalg.inv(self.A)
        self.tri = transform(arrays['triangles'].reshape(-1, 3), self.A).reshape(-1, 3, 3) if descriptor['backend'] == 'IOS_SURFACE' else None

    @classmethod
    def open(cls, path):
        meta = load(path)
        file = Path(meta['arrays']['path'])
        if sha(file) != meta['arrays']['sha256']:
            raise ValueError('FIELD_ARRAY_HASH')
        return cls(meta, np.load(file))

    def query(self, points, *, patient_id, frame_id, unit='mm', theta=None):
        m = self.meta
        for (key, value) in [('patient_id', patient_id), ('frame_id', frame_id), ('unit', unit)]:
            if m[key] != value:
                raise ValueError('QUERY_MISMATCH:' + key)
        points = np.asarray(points, float)
        if points.ndim != 2 or points.shape[1] != 3 or (not np.isfinite(points).all()):
            raise ValueError('FINITE_XYZ_REQUIRED')
        n = len(points)
        theta = theta or {'pose_delta_mm': 0.0}
        delta = float(theta.get('pose_delta_mm', 0))
        if not np.isfinite(delta):
            raise ValueError('NONFINITE_THETA')
        if m['backend'] == 'IOS_SURFACE':
            (h, face, bary) = ray_height(self.tri, points[:, :2], m['jaw'] == 'upper')
            if m['jaw'] == 'upper':
                h = h + delta
            known = np.isfinite(h)
            labs = np.full((n, 3), -1, int)
            labs[known] = self.arrays['facet_vertex_labels'][face[known]]
            owner = np.where(np.all(labs == labs[:, 0, None], axis=1), labs[:, 0], -1)
            owner_known = known & (owner > 0)
            gap = h - points[:, 2] if m['jaw'] == 'upper' else points[:, 2] - h
            result = dict(region_id=owner, known_geometry=known, known_region=owner_known, source_address=face, barycentric=bary, vertical_gap_mm=gap, height_mm=h, quantity_kind='SIGNED_VERTICAL_GAP_NOT_SDF', volumetric_inside='UNKNOWN_OPEN_IOS_SURFACE')
        elif m['backend'] == 'CBCT_VOXELS':
            if delta != 0:
                raise ValueError('IOS_POSE_CANNOT_MOVE_CBCT')
            index = transform(points, self.inv)
            nearest = np.floor(index + 0.5).astype(np.int64)
            dims = np.array(self.arrays['owner'].shape)
            known = np.all((nearest >= 0) & (nearest < dims), axis=1)
            on_face = np.any(np.abs(index + 0.5 - np.round(index + 0.5)) < 1e-10, axis=1)
            known &= ~on_face
            owner = np.full(n, -1, int)
            owner[known] = self.arrays['owner'][tuple(nearest[known].T)]
            sdf = np.full(n, np.nan)
            if 'tooth_center_edt_mm' in self.arrays:
                sdf[known] = self.arrays['tooth_center_edt_mm'][tuple(nearest[known].T)]
            result = dict(region_id=owner, known_geometry=known, known_region=known, source_address=nearest, index_coordinate_zyx=index, grid_center_edt_mm=sdf, quantity_kind='VOXEL_OWNER_AND_NEAREST_GRID_CENTER_EDT_NOT_CONTINUOUS_SDF')
        else:
            raise ValueError('UNKNOWN_BACKEND')
        return dict(patient_id=patient_id, frame_id=frame_id, unit='mm', resolution='PER_POINT', time_scale='SIMULTANEOUS', theta_id=digest(theta), field_id=m['field_id'], source_sha256=m['source_sha256'], **result)

def query_modality(fields, modality, points, patient_id, frame_id):
    candidates = [f for f in fields if f.meta['modality'] == modality and f.meta['patient_id'] == patient_id]
    if not candidates:
        return dict(status='ABSTAIN_MISSING_MODALITY', patient_id=patient_id, modality=modality, known_region=[False] * len(points), region_id=[None] * len(points), binding_quantity='Same-subject acquired ' + modality + ' with explicit registration')
    return candidates[0].query(points, patient_id=patient_id, frame_id=frame_id)

def save_field(out, data_dir, meta, arrays):
    data_dir.mkdir(parents=True, exist_ok=True)
    budget()
    arraypath = data_dir / (meta['field_id'] + '.npz')
    np.savez_compressed(arraypath, **arrays)
    descriptor = {**meta, 'schema': 'Patient360-RegionField-v2', 'unit': 'mm', 'resolution': 'PER_POINT', 'time_scale': 'SIMULTANEOUS', 'arrays': source(arraypath), 'physical_accuracy': 'UNKNOWN'}
    path = out / (meta['field_id'] + '.json')
    dump(path, descriptor)
    return (RegionField(descriptor, arrays), path)

def validate_bridge(bridge, source_field, target_field):
    """Fit is supplied, held-out residuals checked here; neither is a uniform error proof."""
    (a, b) = (source_field.meta, target_field.meta)
    if a['patient_id'] != b['patient_id'] or bridge.get('patient_id') != a['patient_id']:
        raise ValueError('CROSS_PATIENT_JOIN')
    if bridge.get('unit') != 'mm':
        raise ValueError('BRIDGE_UNITS')
    for (name, field) in [('source', a), ('target', b)]:
        if bridge.get(name + '_frame_id') != field['frame_id']:
            raise ValueError('BRIDGE_FRAME')
        if bridge.get(name + '_source_sha256') != field['source_sha256']:
            raise ValueError('STALE_REGISTRATION')
    evidence = bridge.get('identity_evidence', {})
    if not evidence.get('locator') or not evidence.get('sha256'):
        raise ValueError('IDENTITY_EVIDENCE_REQUIRED')
    if bridge.get('evidence_kind') != 'our_own_fixture':
        path = evidence.get('path')
        if not path or not Path(path).is_file() or sha(path) != evidence['sha256']:
            raise ValueError('IDENTITY_ARTIFACT_REQUIRED')
        identity = load(path)
        if identity.get('patient_id') != a['patient_id'] or identity.get('source_sha256') != a['source_sha256'] or identity.get('target_sha256') != b['source_sha256']:
            raise ValueError('IDENTITY_CROSSWALK_MISMATCH')
        if identity.get('pairing_status') != 'SOURCE_DECLARED_SAME_SUBJECT':
            raise ValueError('UNVERIFIED_SUBJECT_PAIR')
    if bridge.get('time_scale') not in ['SIMULTANEOUS', 'HANDOVER'] or not bridge.get('acquisition_state'):
        raise ValueError('ACQUISITION_STATE_REQUIRED')
    A = affine_valid(bridge['source_to_target'], rigid=True)
    hold = bridge.get('heldout_landmarks', {})
    fit = set(bridge.get('fit_landmark_ids', []))
    ids = hold.get('ids', [])
    if len(ids) < 3 or fit.intersection(ids) or len(set(ids)) != len(ids):
        raise ValueError('INDEPENDENT_LANDMARKS_REQUIRED')
    p = np.asarray(hold['source_xyz_mm'], float)
    q = np.asarray(hold['target_xyz_mm'], float)
    if p.shape != q.shape or p.shape != (len(ids), 3) or (not np.isfinite(p).all()) or (not np.isfinite(q).all()):
        raise ValueError('LANDMARK_FORMAT')
    if np.linalg.matrix_rank(p - p.mean(0), tol=1e-08) < 2:
        raise ValueError('COLLINEAR_LANDMARKS')
    err = np.linalg.norm(transform(p, A) - q, axis=1)
    tol = bridge.get('heldout_acceptance_mm')
    if tol is None or not np.isfinite(tol) or tol < 0 or (np.max(err) > tol):
        raise ValueError('HELDOUT_REGISTRATION_RESIDUAL')
    bound = bridge.get('uniform_registration_error_bound_mm')
    if bound is not None and (not np.isfinite(bound) or bound < np.max(err) or (not bridge.get('bound_evidence_locator'))):
        raise ValueError('UNSUPPORTED_REGISTRATION_BOUND')
    return dict(status='VALID_COORDINATE_CONNECTOR', maximum_heldout_residual_mm=float(np.max(err)), uniform_registration_error_bound_mm=bound, physical_registration='UNKNOWN' if bound is None else 'CONDITIONAL_ON_BOUND_EVIDENCE', no_physical_pair_claim=bridge.get('evidence_kind') == 'our_own_fixture')

def registered_query(bridge, source_field, target_field, points):
    verdict = validate_bridge(bridge, source_field, target_field)
    result = target_field.query(transform(points, bridge['source_to_target']), patient_id=target_field.meta['patient_id'], frame_id=target_field.meta['frame_id'])
    result['registration'] = verdict
    return result
