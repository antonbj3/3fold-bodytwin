"""X18b source carrier + X1B volumetric virtual preparation, with explicit closures."""
import sys, time, resource, importlib.util
from pathlib import Path
import numpy as np
from scipy import ndimage
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'vendor'))
sys.path.insert(0, str(Path(__file__).parent / 'legacy/vendor'))
from util import *
from source_ops import pair
from legacy.geometry import height, export_stl, parse_stl, unit

def source_carrier(a, fdi):
    ids = np.flatnonzero(a['owner'] == fdi)
    if len(ids) < 100:
        raise ValueError('fewer than100 predicted source faces')
    ff = a['f'][ids]
    (vi, inv) = np.unique(ff, return_inverse=True)
    return (a['v'][vi].copy(), inv.reshape(-1, 3), ids, vi)

def build():
    import crown_fit_geometry as G
    spec = importlib.util.spec_from_file_location('gencad_x1b', ROOT / 'vendor/crown_design_geometry.py')
    CDG = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CDG)
    CDG.G = G
    prereg = read(ROOT / 'PREREG_R2.json')
    q = prereg['parameters']
    sel = read(PAYLOAD / 'private/PREP_SELECTION.json')['payload']
    cohort = read(PAYLOAD / 'private/COHORT.json')['payload']['cases']
    lookup = {r['case_key']: r for r in cohort}
    records = []
    start = time.perf_counter()
    for s in sel:
        key = s['case_key']
        dest = PAYLOAD / 'preparations' / key
        dest.mkdir(parents=True, exist_ok=True)
        recfile = dest / 'record.json'
        if recfile.exists():
            records.append(read(recfile))
            continue
        tic = time.perf_counter()
        record = dict(case_key=key, split=s['split'], fdi=36, resolution='PER_POINT', time_scale='HANDOVER')
        try:
            (p, sources) = pair(lookup[key])
            (v, f, face_ids, vertex_ids) = source_carrier(p['lower'], 36)
            frame = p['lower']['meta']['frame']
            z = unit(frame['occlusal_unit'])
            x = unit(frame['right_unit'])
            x = unit(x - z * (x @ z))
            y = unit(np.cross(z, x))
            R = np.stack([x, y, z], axis=1)
            center = v.mean(0)
            local = (v - center) @ R
            tri = local[f]
            h = q['grid_spacing_mm']
            lo = local.min(0) - h * 3
            hi = local.max(0) + h * 3
            grid = G.Grid(lo, hi, h)
            if np.prod(grid.shape) > q['max_grid_nodes']:
                raise ValueError('bounded grid node limit')
            ax = grid.axes()
            (xx, yy) = np.meshgrid(ax[0], ax[1], indexing='ij')
            xy = np.c_[xx.ravel(), yy.ravel()]
            roof = height(tri, xy).reshape(xx.shape)
            base = float(tri[:, :, 2].min())
            zz = ax[2][None, None, :]
            occupied = np.isfinite(roof)[:, :, None] & (zz >= base) & (zz <= roof[:, :, None])
            if not occupied.any():
                raise ValueError('empty virtual envelope')
            phi = (ndimage.distance_transform_edt(~occupied, sampling=h) - ndimage.distance_transform_edt(occupied, sampling=h)).astype(np.float32)
            margin = base + q['virtual_margin_above_base_mm']
            case = dict(grid=grid, phi_T=phi, z_m=margin)
            scratch = dest / 'construction'
            scratch.mkdir(exist_ok=True)
            CDG.SCRATCH = str(scratch)
            (prep, margin_pts) = CDG.preparation(case, q['occlusal_reduction_mm'], q['axial_reduction_mm'], q['spacer_mm'], cache=False)
            (pv, pf) = G.surface_mesh(prep, grid)
            volume = float(np.einsum('ij,ij->i', pv[pf[:, 0]], np.cross(pv[pf[:, 1]], pv[pf[:, 2]])).sum() / 6)
            if volume < 0:
                pf = pf[:, ::-1]
            world = pv @ R.T + center
            export_stl(dest / 'preparation.stl', world, pf)
            export_stl(dest / 'native_exterior.stl', v, f)
            np.savez_compressed(dest / 'fields.npz', source_vertices=v, source_faces=f, source_face_ids=face_ids, source_vertex_ids=vertex_ids, phi=phi, preparation=prep, origin=grid.origin, step=h, R=R, center=center, margin_pts=margin_pts, preparation_vertices=pv, preparation_faces=pf)
            for cache in scratch.glob('*.npz'):
                cache.unlink()
            scratch.rmdir()
            record.update(status='EXPORTED', source=sources, grid_shape=grid.shape, virtual_base_z_mm=base, virtual_margin_z_mm=margin, original_faces=len(f), preparation_faces=len(pf), closures=['X11 FDI unvalidated in target domain', 'Basal plane and insertion envelope', '3deg wall taper and design reductions', 'Grid distance reinitialization uses candidate-triangle approximation; continuous SDF error UNKNOWN'], physical_preparation_reference='UNKNOWN_NO_PREPARED_TOOTH_SCAN')
        except (ValueError, KeyError, RuntimeError) as e:
            record.update(status='FAILED_PREPARATION', reason=str(e))
        record['seconds'] = time.perf_counter() - tic
        dump(recfile, record)
        records.append(record)
        print('preparation', key, record['status'], round(record['seconds'], 2), flush=True)
        budget()
    files = {str(p.relative_to(PAYLOAD)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted((PAYLOAD / 'preparations').rglob('*')) if p.is_file()}
    freeze(ROOT / 'FROZEN_PREPARATIONS.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R2.json'), files=files, records=records, measurement_status='No physical measurement; independent numerical tests follow freeze'))
    dump(ROOT / 'raw/PREPARATION_BUILD_COST.json', dict(seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

def bad_edges(f):
    edges = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    return int(np.count_nonzero(np.unique(edges, axis=0, return_counts=True)[1] != 2))

def stl_tri(blob):
    n = int.from_bytes(blob[80:84], 'little')
    dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    if len(blob) != 84 + 50 * n:
        raise ValueError('invalid exported STL')
    return np.frombuffer(blob, dtype=dt, count=n, offset=84)['v'].astype(float)

def iso_error(world, a):
    local = (world - a['center']) @ a['R']
    ind = ((local - a['origin']) / float(a['step'])).T
    values = ndimage.map_coordinates(a['preparation'], ind, order=1, mode='constant', cval=100.0)
    return float(np.max(np.abs(values)))

def evaluate(bundle):
    frozen = bundle.json('FROZEN_PREPARATIONS.json')['payload']
    pr = bundle.json('PREREG_R2.json')
    rows = []
    for r in frozen['records']:
        row = {k: r[k] for k in ['case_key', 'split', 'status']}
        if r['status'] != 'EXPORTED':
            row['reason'] = r['reason']
            rows.append(row)
            continue
        prefix = 'payload/preparations/' + r['case_key'] + '/'
        a = bundle.npz(prefix + 'fields.npz')
        truth = a['source_vertices'][a['source_faces']]
        ext = stl_tri(bundle.bytes(prefix + 'native_exterior.stl'))
        tri = stl_tri(bundle.bytes(prefix + 'preparation.stl'))
        (verts, ix) = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
        f = ix.reshape(-1, 3)
        exterr = float(np.max(np.abs(truth - ext)))
        residual = iso_error(verts, a)
        nbad = bad_edges(f)
        volume = float(np.einsum('ij,ij->i', tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6)
        gates = dict(exterior=exterr <= pr['metrics']['exact_exterior_export_max_mm'], iso=residual <= pr['metrics']['preparation_stl_iso_residual_max_mm'], manifold=nbad == 0, positive_volume=volume > 0)
        fault = dict(shifted_exterior_rejected=float(np.max(np.abs(truth - (ext + [0, 0, 0.5])))) > pr['metrics']['exact_exterior_export_max_mm'], missing_facet_rejected=bad_edges(f[:-1]) > 0, shifted_preparation_rejected=iso_error(verts + [0, 0, 0.5], a) > pr['metrics']['preparation_stl_iso_residual_max_mm'])
        rows.append(dict(**row, gates=gates, fault_rejections=fault, exterior_export_max_mm=exterr, iso_residual_max_mm=residual, nonmanifold_edges=nbad, signed_volume_mm3=volume, physical_preparation='UNKNOWN', resolution='PER_POINT'))
    good = [r for r in rows if r['status'] == 'EXPORTED']
    passes = sum((all(r['gates'].values()) for r in good))
    out = dict(claim_type='capability', requested=len(rows), exported=len(good), geometric_pass=passes, failed=len(rows) - passes, exclusion_fraction=(len(rows) - len(good)) / len(rows), rows=rows, gate=passes >= pr['metrics']['minimum_successful_specimens'] and all((all(r['fault_rejections'].values()) for r in good)), external_referent=pr['external_referent'], scope='Virtual full-volume preparation plus exact native surface carrier; no measured preparation, physical fit or clinical restoration validation')
    dump(bundle.payload / 'evaluation/PREPARATION_SCORES.json', out)
    return out
if __name__ == '__main__':
    build()
