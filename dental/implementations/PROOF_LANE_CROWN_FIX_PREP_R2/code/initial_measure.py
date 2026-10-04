from common import *
rows = []
for rec in inputs():
    (v, f, info) = native(rec)
    old = read(PREV / 'raw/R3_GENERATION.json')
    rr = next((x for x in old if x['key'] == rec['key'] and x['status'] == 'GENERATED'), None)
    M = np.array(rr['construction']['M']) if rr else np.eye(3)
    S = load(rec['private_path'])['target']
    m = trimesh.Trimesh(v, f, process=False)
    bd = __import__('evaluate').boundary_edges(triangles_mesh(S)[1])
    (sv, sf) = triangles_mesh(S)
    z = (sv[np.unique(bd)] @ M.T)[:, 2]
    a = dict(key=rec['key'], family=rec['family'], native_volume_mm3=float(m.volume), native_euler=int(m.euler_number), native_components=len(m.split(only_watertight=False)), native_z_bounds_mm=(v @ M.T)[:, 2].min().item(), boundary_z_quantiles_mm=np.quantile(z, [0, 0.25, 0.5, 0.75, 1]), old_core_base_mm=rr['construction']['b'] if rr else None, old_top_mm=rr['construction']['H'] if rr else None, resolution='PER_TOOTH', source=info)
    rows.append(a)
    print(rec['key'], a['native_euler'], a['boundary_z_quantiles_mm'], flush=True)
dump(R / 'raw/INITIAL_LOCAL_MEASUREMENT.json', dict(timestamp=now(), rows=rows, scope='Observed cervical scan boundary; not measured CEJ or prepared finish line'))
