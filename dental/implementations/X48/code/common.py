"""X48: row point arrays, column homogeneous poses, millimetres and degrees."""
from dental_release.paths import expand as _release_expand
import json, hashlib, datetime, os
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
ROOT = Path(__file__).resolve().parents[1]
DENTAL = ROOT.parents[1]
X37 = DENTAL / 'results/LANE_X37_ALIGNER_MOVEMENT'
X19 = DENTAL / 'results/LANE_X19_ALIGNER_FORCE'
DATA = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/AlignerMovement_Zenodo11280343'))
CACHE = Path(_release_expand('@DENTAL_WORK_ROOT@/X37_ALIGNER_MOVEMENT'))
OWN_DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X48_ALIGNER_MEASURED'))
REFS = ['all_crowns', 'bilateral_molars', 'left_molars', 'right_molars']

def read(p):
    return json.loads(Path(p).read_text())

def write(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(2 ** 20), b''):
            h.update(b)
    return h.hexdigest()

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def state(phase, gate, next_operation):
    write(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X48-aligner-measured', claim_type='information_link', updated_utc=now(), phase=phase, latest_gate=gate, next_operation=next_operation, review_state='PENDING_INDEPENDENT_REVIEW', resource_limits=dict(threads=4, max_intermediate_bytes=3000000000, gpu=False)))

def apply(H, p):
    return np.asarray(p) @ H[:3, :3].T + H[:3, 3]

def angle(H):
    return float(np.degrees(Rotation.from_matrix(H[:3, :3]).magnitude()))

def pose_index():
    rows = read(X37 / 'raw/R1_POSES.json')
    return (rows, {(r['patient'], r['jaw'], r['source'], r['week'], r['fdi']): r for r in rows})

def axes(ix, p, j, f):
    ids = sorted({k[4] for k in ix if k[:4] == (p, j, 'Sirona', 0)})
    c = np.array(ix[p, j, 'Sirona', 0, f]['centroid_base_mm'])
    (q, n) = divmod(f, 10)
    toward = [k for k in ids if k // 10 == q and k % 10 < n]
    other = (q + 1 if q % 2 else q - 1) * 10 + 1
    g = max(toward) if toward else other
    m = np.array(ix[p, j, 'Sirona', 0, g]['centroid_base_mm']) - c
    m[1] = 0
    m /= np.linalg.norm(m)
    b = np.array([-m[2], 0, m[0]])
    center = np.mean([ix[p, j, 'Sirona', 0, k]['centroid_base_mm'] for k in ids], axis=0)
    if np.dot(b, c - center) < 0:
        b = -b
    return {'mesial': m, 'buccal': b, 'intrusion': np.array([0.0, 1.0, 0.0])}

def delta(H1, H0, c, ax):
    t = apply(H1, c) - apply(H0, c)
    w = np.degrees(Rotation.from_matrix(H1[:3, :3] @ H0[:3, :3].T).as_rotvec())
    return ({'mesiodistal_translation': float(t @ ax['mesial']), 'buccolingual_translation': float(t @ ax['buccal']), 'intrusion_extrusion': float(t @ ax['intrusion']), 'buccolingual_tip': float(w @ ax['mesial']), 'mesiodistal_tip': float(w @ ax['buccal']), 'axial_rotation': float(w @ ax['intrusion'])}, t, w)

def ratio_box(a, p, ua, up):
    """Exact corner extrema for a/p on a rectangular input set avoiding p=0.
 Bounds are conditional on supplied deterministic input radii, not confidence.
 """
    (lo, hi) = (p - up, p + up)
    if lo <= 0 <= hi:
        return None
    vals = [x / y for x in [a - ua, a + ua] for y in [lo, hi]]
    return [float(min(vals)), float(max(vals))]
