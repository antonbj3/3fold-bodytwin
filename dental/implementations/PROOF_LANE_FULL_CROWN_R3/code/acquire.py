from common_r3 import *
import time
sys.path.insert(0, str(OLD / 'code'))
from score import VF, grid
st = time.perf_counter()
out = DATA / 'inputs'
out.mkdir(exist_ok=True)
rows = []
for rec in read(ROOT / 'FROZEN_COHORT.json')['rows']:
    key = rec['key']
    p = npz(PUBLIC / (key + '.npz'))
    a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
    tri = a['source_triangles']
    tri = tri[tri.mean(1)[:, 2] >= p['margin_z']]
    (xy, faces) = grid(tri)
    ceil = VF.height(a['antagonist_triangles'], xy, True)
    rg = ceil - VF.height(tri, xy)
    m = npz(OLD_DATA / 'R2/outer_first_joint' / key / 'mesh.npz')
    ext = m['vertices'][m['faces'][m['face_roles'] == 0]]
    pg = ceil - VF.height(ext, xy)
    np.savez_compressed(out / (key + '.npz'), **p, prior_vertices=m['vertices'], prior_faces=m['faces'], prior_roles=m['face_roles'], contact_xy=xy, contact_faces=faces, contact_gap_mm=np.clip(rg, -0.25, 0.35), contact_reference_gap_unclipped=rg, contact_ceiling_mm=ceil)
    rows.append(dict(rec, measurement_finite=int(np.isfinite(rg).sum()), measurement_scope='Digital extraction from same native scan; not a new lab measurement. Full axial source withheld. Unclipped gap retained only for validation; generator receives clipped field.'))
    print(key, len(xy), int(np.isfinite(rg).sum()), flush=True)
dump(out / 'RECORDS.json', rows)
freeze(ROOT / 'FROZEN_ACQUISITION_A.json', dict(files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()}, rows=rows, seconds=time.perf_counter() - st, role='Declared added input, not untouched validation data'))
