"""Prepared-surface information track. No original tooth/reference mount."""
import sys, json, time
from pathlib import Path
import numpy as np
from scipy import ndimage
from skimage.measure import marching_cubes
sys.path.insert(0, '/participants')
from legacy.geometry import export_stl

def run():
    cfg = json.loads(Path('/runner/PREP_CONFIG.json').read_text())
    q = cfg['parameters']
    records = json.loads(Path('/whole_inputs/RECORDS.json').read_text())
    out = Path('/output')
    rows = []
    for rec in records:
        if rec['status'] != 'PREPARED':
            for name in ['prep_uniform', 'prep_regional']:
                rows.append(dict(key=rec['key'], participant=name, status='UNKNOWN_PREPARATION'))
            continue
        with np.load(Path('/whole_inputs', rec['key'], 'preparation.npz')) as a:
            p = dict(a)
        prep = p['preparation']
        h = float(p['step'])
        origin = p['origin']
        z = origin[2] + h * np.arange(prep.shape[2])
        cut = np.broadcast_to(float(p['margin_z']) - z[None, None, :], prep.shape)
        grad = np.gradient(prep, h)
        nz = abs(grad[2]) / np.maximum(np.sqrt(sum((g * g for g in grad))), 1e-12)
        offsets = {'prep_uniform': q['uniform_dilation_mm'], 'prep_regional': q['axial_dilation_mm'] + (q['occlusal_dilation_mm'] - q['axial_dilation_mm']) * np.clip(nz, 0, 1) ** 2}
        for (name, delta) in offsets.items():
            st = time.perf_counter()
            rr = dict(key=rec['key'], case_key=rec['case_key'], family=rec['family'], participant=name, status='FAILED')
            outer = prep - delta
            field = np.maximum.reduce([outer, -p['cavity'], cut]).astype(np.float32)
            try:
                (v, f, _, _) = marching_cubes(field, q['iso_level_mm'], spacing=(h,) * 3, allow_degenerate=False)
                v += origin
                if np.einsum('ij,ij->i', v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]])).sum() < 0:
                    f = f[:, ::-1]
                pp = v[f].mean(1)
                ix = ((pp - origin) / h).T
                roles = np.argmax(np.stack([ndimage.map_coordinates(outer, ix, order=1), ndimage.map_coordinates(-p['cavity'], ix, order=1), float(p['margin_z']) - pp[:, 2]]), axis=0).astype(np.int8)
                dest = out / name / rec['key']
                dest.mkdir(exist_ok=True, parents=True)
                np.savez_compressed(dest / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
                export_stl(dest / 'crown.stl', v @ p['source_R'].T + p['source_base'], f)
                rr.update(status='EXPORTED', vertices=len(v), faces=len(f))
            except (ValueError, RuntimeError) as e:
                rr['reason'] = str(e)
            rr['seconds'] = time.perf_counter() - st
            rows.append(rr)
        print('prepared-field crowns', rec['key'], flush=True)
    (out / 'RECORDS.json').write_text(json.dumps(rows, indent=2) + '\n')
if __name__ == '__main__':
    run()
