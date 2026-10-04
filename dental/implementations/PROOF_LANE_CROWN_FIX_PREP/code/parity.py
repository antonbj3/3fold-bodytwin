from core import *
rows = []
for tag in ['R2', 'R3']:
    for r in read(R / f'FROZEN_PREDICTIONS_{tag}.json')['rows']:
        if r['status'] != 'GENERATED':
            continue
        m = load(r['mesh_path'])
        p = sample(m['vertices'][m['faces'][m['roles'] == 1]], 64)
        q = trimesh.Trimesh(m['prep_vertices'], m['prep_faces'], process=False)
        a = Distance(q.vertices, q.faces).query(p)[0]
        (_, b, _) = trimesh.proximity.closest_point(q, p)
        err = float(abs(a - b).max())
        rows.append(dict(key=r['key'], material=r['material'], round=tag, points=64, max_backend_error_mm=err, pass_gate=err <= 1e-08))
        print(tag, r['key'], r['material'], err, flush=True)
dump(R / 'raw/PARITY.json', dict(rows=rows, all_pass=all((r['pass_gate'] for r in rows)), max_error_mm=max((r['max_backend_error_mm'] for r in rows)), count=len(rows), resolution='PER_POINT'))
