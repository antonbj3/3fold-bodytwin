import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '4'
import json, hashlib, pathlib, time, datetime, resource, numpy as np, trimesh, argparse
from sweep import *
R = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def dump(p, j):
    p.write_text(json.dumps(j, indent=2, allow_nan=False) + '\n')

def check_hashes(co):
    for r in co['rows']:
        assert sha(r['mesh_path']) == r['mesh_sha256']
        assert sha(r['public_path']) == r['public_sha256']
    for r in co['exports']:
        assert sha(r['mesh_path']) == r['mesh_sha256']
        assert sha(r['stl_path']) == r['stl_sha256']
        assert sha(r['three_mf_path']) == r['three_mf_sha256']

def health(t):
    if not len(t):
        return dict(triangles=0, watertight=False)
    m = trimesh.Trimesh(t.reshape(-1, 3), np.arange(len(t) * 3).reshape(-1, 3), process=True)
    return dict(triangles=len(t), watertight=bool(m.is_watertight), winding_consistent=bool(m.is_winding_consistent))

def run(outdir, D=None):
    resource.setrlimit(resource.RLIMIT_AS, (3500 * 1024 ** 2, 3500 * 1024 ** 2))
    outdir.mkdir(parents=True, exist_ok=True)
    co = json.loads((R / 'FROZEN_COHORT.json').read_text())
    check_hashes(co)
    D = np.array(co['path']['translation_mm'] if D is None else D, float)
    tick = time.perf_counter()
    rows = []
    for rec in co['rows']:
        with np.load(rec['mesh_path']) as m, np.load(rec['public_path']) as p:
            v = m['vertices']
            f = m['faces']
            roles = m['roles']
            tri = v[f]
            results = {}
            statuses = []
            for side in ['mesial', 'distal']:
                obs = Obstacle(p[side])
                hp = health(p[side])
                order = np.r_[np.flatnonzero(roles == 0), np.flatnonzero(roles == 2), np.flatnonzero(roles == 1)]
                path = check_surface(tri[order], obs, D, face_ids=order, roles=roles[order], budget_s=180)
                endpoint = check_surface(tri[order], obs, D * 0, face_ids=order, roles=roles[order], budget_s=180)
                if path['status'] == 'COLLISION':
                    w = path['witness']
                    w['obstacle_region'] = side
                    path['inside_only_repair'] = 'IMPOSSIBLE_FOR_UNCHANGED_EXTERIOR' if w['source_role'] in [0, 2] else 'UNKNOWN until constructive repair'
                    (outdir / 'witnesses').mkdir(exist_ok=True)
                    np.savez_compressed(outdir / 'witnesses' / f"{rec['key']}_{side}.npz", source_triangle=np.array(w['source_triangle_mm']), obstacle_triangle=np.array(w['obstacle_triangle_mm']), point_mm=np.array(w['point_mm']), translation_mm=D, pose_s=np.array(w['s_float']))
                statuses.append(path['status'])
                results[side] = dict(path=path, endpoint=endpoint, obstacle_health=hp)
            status = 'COLLISION' if 'COLLISION' in statuses else 'UNKNOWN'
            row = dict(key=rec['key'], family=rec['family'], claim_type='capability', resolution='PER_POINT', timescale='SIMULTANEOUS', translation_mm=D.tolist(), crown_faces=len(f), mesh_sha256=rec['mesh_sha256'], full_path_status=status, neighbours=results, preparation_status='UNKNOWN separate prep solid/spacer absent; cone-boundary seated contact not called clearance', physical_insertion='UNKNOWN', inside_only_repair_status='IMPOSSIBLE_FOR_UNCHANGED_EXTERIOR' if any((q['path'].get('inside_only_repair') == 'IMPOSSIBLE_FOR_UNCHANGED_EXTERIOR' for q in results.values())) else 'UNKNOWN')
            rows.append(row)
            dump(outdir / 'COHORT.json', dict(rows=rows, elapsed_s=time.perf_counter() - tick, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, complete=len(rows) == 18, code_sha256={p.name: sha(p) for p in (R / 'code').glob('*.py')}))
            dump(R / 'CURRENT_WORK_STATE.json', dict(lane='X95-whole-insertion-path', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='COHORT_RUNNING', completed=len(rows), latest_gate={rec['key']: status}, next_operation='Finish18cases, preserve R1, change repair/path representation'))
            print(len(rows), rec['key'], status, [q['path']['status'] for q in results.values()], round(time.perf_counter() - tick, 2), flush=True)
    return rows
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(R / 'rounds/R1'))
    a = ap.parse_args()
    run(pathlib.Path(a.out))
