"""Whole-crown lift. Mounts preparation and frozen predictions, never original tooth."""
import sys, json, time
from pathlib import Path
import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay
from scipy import ndimage
from skimage.measure import marching_cubes
sys.path.insert(0, '/participants')
from legacy.geometry import export_stl

def run():
    spec = json.loads(Path('/runner/WHOLE_CONFIG.json').read_text())
    q = spec['parameters']
    records = json.loads(Path('/whole_inputs/RECORDS.json').read_text())
    names = json.loads(Path('/inputs/PARTICIPANTS.json').read_text())['participants']
    out = Path('/output')
    results = []
    for rec in records:
        key = rec['key']
        case = rec['case_key']
        fam = rec['family']
        ts = json.loads(Path('/inputs/tasks', case + '.json').read_text())
        i = next((i for (i, t) in enumerate(ts) if t['family'] == fam and t['level'] == spec['level']))
        t = ts[i]
        if rec['status'] != 'PREPARED' or t['status'] != 'READY':
            for name in names:
                results.append(dict(key=key, participant=name, status='UNKNOWN_PREPARATION', reason=rec.get('reason', t.get('site_error'))))
            continue
        with np.load(Path('/whole_inputs', key, 'preparation.npz')) as a:
            p = dict(a)
        with np.load(Path('/inputs', t['geometry_file'])) as a:
            s = dict(a)
        h = float(p['step'])
        origin = p['origin']
        shape = p['preparation'].shape
        axes = [origin[k] + h * np.arange(shape[k]) for k in range(3)]
        (xx, yy) = np.meshgrid(axes[0], axes[1], indexing='ij')
        xy = np.c_[xx.ravel(), yy.ravel()]
        zz = axes[2][None, None, :]
        mask = Delaunay(s['xy'] * q['side_scale']).find_simplex(xy) >= 0
        for name in names:
            start = time.perf_counter()
            rr = dict(key=key, case_key=case, family=fam, participant=name, status='FAILED')
            with np.load(Path('/roofs', name, case + '.npz')) as a:
                if a['status'][i] != 'DESIGN':
                    rr['status'] = str(a['status'][i])
                    results.append(rr)
                    continue
                z = a['outer_' + str(i)]
            try:
                roof = LinearNDInterpolator(s['xy'], z)(xy)
                bad = ~np.isfinite(roof)
                roof[bad] = NearestNDInterpolator(s['xy'], z)(xy[bad])
                roof = roof.reshape(xx.shape)
                occ = mask.reshape(xx.shape)[:, :, None] & (zz >= float(p['margin_z']) - 1.0) & (zz <= roof[:, :, None])
                outer = (ndimage.distance_transform_edt(~occ, sampling=h) - ndimage.distance_transform_edt(occ, sampling=h)).astype(np.float32)
                cut = np.broadcast_to(float(p['margin_z']) - zz, shape)
                field = np.maximum.reduce([outer, -p['cavity'], cut]).astype(np.float32)
                if field.min() >= q['iso_level_mm']:
                    raise ValueError('no crown volume above prescribed iso level')
                (v, f, _, _) = marching_cubes(field, q['iso_level_mm'], spacing=(h,) * 3, allow_degenerate=False)
                v += origin
                if np.einsum('ij,ij->i', v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]])).sum() < 0:
                    f = f[:, ::-1]
                pts = v[f].mean(1)
                ix = ((pts - origin) / h).T
                vals = np.stack([ndimage.map_coordinates(outer, ix, order=1), ndimage.map_coordinates(-p['cavity'], ix, order=1), pts[:, 2] * 0 + float(p['margin_z']) - pts[:, 2]])
                roles = np.argmax(vals, axis=0).astype(np.int8)
                dest = out / name / key
                dest.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(dest / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
                export_stl(dest / 'crown.stl', v @ p['source_R'].T + p['source_base'], f)
                rr.update(status='EXPORTED', vertices=len(v), faces=len(f), exterior_faces=int(np.sum(roles == 0)), intaglio_faces=int(np.sum(roles == 1)), margin_faces=int(np.sum(roles == 2)))
            except (ValueError, RuntimeError) as e:
                rr['reason'] = str(e)
            rr['seconds'] = time.perf_counter() - start
            results.append(rr)
        print('whole crowns', key, 'exports', sum((r['status'] == 'EXPORTED' for r in results)), 'requests', len(results), flush=True)
        (out / 'RECORDS.json').write_text(json.dumps(results, indent=2) + '\n')
if __name__ == '__main__':
    run()
