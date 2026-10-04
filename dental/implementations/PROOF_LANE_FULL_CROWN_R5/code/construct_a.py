from dental_release.paths import expand as _release_expand
import sys, os, time, json, hashlib, datetime, resource
from pathlib import Path
R = Path(__file__).resolve().parents[1]
BASE = R.parent
sys.path.insert(0, str(BASE / 'PROOF_LANE_FULL_CROWN_R4/code'))
from geometry import *
import score as scorer
D = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN_R5'))
D.mkdir(exist_ok=True)

def save(path, x):
    dump(path, x)

def checkpoint(phase, gate, nxt):
    save(R / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-full-crown-r5', claim_type='capability', phase=phase, latest_gate=gate, next_operation=nxt, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))

def inputs():
    lock = {r['key']: r for r in read(BASE / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']}
    return [(rec, lock[rec['key']]) for rec in read(BASE / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_PREDICTIONS_C.json')['rows'] if rec['method'] == 'exact_margin']

def normal_offset(v, f, roles):
    ids = np.unique(f[roles == 1])
    apex = ids[-1]
    out = v.copy()
    t = v[f[roles == 1]]
    nn = -np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    nn /= np.maximum(np.linalg.norm(nn, axis=1)[:, None], 1e-30)
    innerfaces = f[roles == 1]
    res = []
    for idx in ids:
        normals = nn[np.any(innerfaces == idx, axis=1)]
        g = 0.05 if idx == apex else 0.025
        delta = np.linalg.lstsq(normals, np.full(len(normals), g), rcond=None)[0]
        out[idx] += delta
        res.extend(abs(normals @ delta - g))
    return (out, dict(max_plane_offset_residual_mm=float(max(res)), maximum_vertex_offset_mm=float(np.linalg.norm(out[ids] - v[ids], axis=1).max()), internal_setting_mm=0.05, marginal_setting_mm=0.025))

def project(v, f, roles):
    out = v.copy()
    ext = trimesh.Trimesh(v, f[roles == 0], process=False)
    bound = loops(ext)[0]
    eids = np.unique(f[roles == 0])
    free = np.setdiff1d(eids, bound)
    inner = out[f[roles == 1]]
    nn = -np.cross(inner[:, 1] - inner[:, 0], inner[:, 2] - inner[:, 0])
    nn /= np.maximum(np.linalg.norm(nn, axis=1)[:, None], 1e-30)
    hist = []
    for it in range(12):
        (q, d, j) = closest(inner, out[free])
        dr = out[free] - q
        sg = np.einsum('ij,ij->i', dr, nn[j])
        direction = dr / np.maximum(d[:, None], 1e-30)
        inside = sg < 0
        direction[inside] = nn[j[inside]]
        active = (d < 0.6) | inside
        desired = q + 0.6 * direction
        out[free[active]] = desired[active]
        hist.append(dict(iteration=it, active_vertices=int(active.sum()), sampled_free_min_mm=float(d.min())))
        if not active.any():
            break
    assert np.array_equal(out[bound], v[bound])
    dis = np.linalg.norm(out[eids] - v[eids], axis=1)
    return (out, dict(history=hist, boundary_identity_error_mm=0.0, max_displacement_mm=float(dis.max()), p95_displacement_mm=float(np.quantile(dis, 0.95)), changed_vertex_fraction=float(np.mean(dis > 0))))

def run():
    st = time.perf_counter()
    rows = []
    for (rec, rr) in inputs():
        assert sha(rec['mesh_path']) == rec['mesh_sha256']
        m = npz(rec['mesh_path'])
        (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
        (ov, offset) = normal_offset(v, f, roles)
        for method in ['offset_only', 'locked_local_thickening']:
            t0 = time.perf_counter()
            vv = ov.copy()
            info = {}
            if method == 'locked_local_thickening':
                (vv, info) = project(vv, f, roles)
            dest = D / 'A' / method / (rec['key'] + '.npz')
            dest.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(dest, vertices=vv, faces=f, roles=roles, prep_triangles=v[f[roles == 1]])
            ext = trimesh.Trimesh(vv, f[roles == 0], process=False)
            curve = ext.vertices[loops(ext)[0]]
            (q, d, j) = closest(vv[f[roles == 1]], curve)
            i = int(d.argmin())
            row = dict(key=rec['key'], family=rec['family'], method=method, mesh_path=dest, mesh_sha256=sha(dest), offset=offset, repair=info, construction_seconds=time.perf_counter() - t0, boundary_witness=dict(outer_point=curve[i], inner_point=q[i], inner_triangle=vv[f[roles == 1]][j[i]], distance_mm=float(d[i]), inner_triangle_index=int(j[i])), resolution='PER_POINT')
            rows.append(row)
            save(R / 'raw/A_GENERATION.json', rows)
            print(rec['key'], method, 'rim distance', d[i], 'offset residual', offset['max_plane_offset_residual_mm'], flush=True)
        checkpoint('A_GENERATION', len(rows), 'Freeze predictions then score all36 requests without changing denominator')
    freeze(R / 'FROZEN_PREDICTIONS_A.json', dict(claim_type='capability', prereg_sha256=sha(R / 'PREREG_A.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
