from construct_b import *
import geometry as geo
geo.closest = fast_nearest
scorer.closest = fast_nearest
rr = read(R / 'raw/B_GENERATION.json')
lock = {r['key']: r for (_, r) in inputs()}
rows = []
for r in rr:
    if 'local_thickening' not in r['method']:
        continue
    m = npz(r['mesh_path'])
    ext = m['vertices'][m['faces'][m['roles'] == 0]]
    target = npz(lock[r['key']]['private_path'])['target']
    met = scorer.metrics(ext, target)
    print(r['key'], met['p95_mm'], 'offseterr', r['offset']['max_plane_offset_residual_mm'], flush=True)
    rows.append(dict(key=r['key'], shape=met))
save(R / 'raw/B_EARLY_DIAGNOSTIC.json', rows)
