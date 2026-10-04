"""Full-volume crown track consuming the actual R3 preparation field."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'vendor'))
import numpy as np
from scipy import ndimage
from skimage.measure import marching_cubes
from util import *
from preparation import bad_edges, stl_tri
from legacy.geometry import export_stl

def build():
    import crown_fit_geometry as G
    pr = read(ROOT / 'PREREG_R4.json')
    par = pr['parameters']
    source = read(ROOT / 'FROZEN_PREPARATIONS_R3.json')['payload']
    prior = read(ROOT / 'rounds/R3.json')
    passed = {r['case_key'] for r in prior['rows'] if all(r.get('gates', {'missing': False}).values())}
    records = []
    for row in source['records']:
        key = row['case_key']
        record = dict(case_key=key, split=row['split'], status='UNAVAILABLE_PREPARATION')
        if key not in passed:
            records.append(record)
            continue
        tic = time.perf_counter()
        p = PAYLOAD / 'preparations_r3' / key / 'fields.npz'
        if sha(p) != source['files']['preparations_r3/' + key + '/fields.npz']['sha256']:
            raise ValueError('R3 input changed')
        a = dict(np.load(p, allow_pickle=False))
        h = float(a['step'])
        grid = G.Grid(a['origin'], a['origin'] + h * (np.array(a['phi'].shape) - 1), h)
        cavity = G.design_cavity(a['preparation'], grid, np.full(a['preparation'].shape, par['spacer_mm']), levels=(0.0,))
        z = grid.axes()[2][None, None, :]
        field = np.maximum.reduce([a['phi'], -cavity, np.broadcast_to(row['virtual_margin_z_mm'] - z, a['phi'].shape)]).astype(np.float32)
        (v, f, _, _) = marching_cubes(field, level=par['iso_level_mm'], spacing=(h,) * 3, allow_degenerate=False)
        v += grid.origin
        volume = np.einsum('ij,ij->i', v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]])).sum() / 6
        if volume < 0:
            f = f[:, ::-1]
        world = v @ a['R'].T + a['center']
        dest = PAYLOAD / 'crowns' / key
        dest.mkdir(parents=True, exist_ok=True)
        export_stl(dest / 'crown.stl', world, f)
        np.savez_compressed(dest / 'prediction.npz', field=field, cavity=cavity, vertices=v, faces=f, margin_z=row['virtual_margin_z_mm'])
        record.update(status='EXPORTED', seconds=time.perf_counter() - tic)
        records.append(record)
        print('crown', key, flush=True)
    files = {str(p.relative_to(PAYLOAD)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted((PAYLOAD / 'crowns').rglob('*')) if p.is_file()}
    freeze(ROOT / 'FROZEN_CROWNS.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R4.json'), source_preparations_sha256=sha(ROOT / 'FROZEN_PREPARATIONS_R3.json'), files=files, records=records, measurement_status='NO_PHYSICAL_MEASUREMENT; source-distance evaluation after freeze'))

def field_values(field, world, a, cval=100.0):
    local = (world - a['center']) @ a['R']
    idx = ((local - a['origin']) / float(a['step'])).T
    return ndimage.map_coordinates(field, idx, order=1, mode='constant', cval=cval)

def evaluate(bundle):
    import trimesh
    pr = bundle.json('PREREG_R4.json')
    par = pr['parameters']
    m = pr['metrics']
    fr = bundle.json('FROZEN_CROWNS.json')['payload']
    rows = []
    start = time.perf_counter()
    for record in fr['records']:
        key = record['case_key']
        row = dict(record)
        if record['status'] != 'EXPORTED':
            rows.append(row)
            continue
        a = bundle.npz('payload/preparations_r3/' + key + '/fields.npz')
        b = bundle.npz('payload/crowns/' + key + '/prediction.npz')
        tri = stl_tri(bundle.bytes('payload/crowns/' + key + '/crown.stl'))
        (world, ix) = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
        f = ix.reshape(-1, 3)
        local = (world - a['center']) @ a['R']
        residual = float(np.max(np.abs(field_values(b['field'], world, a) - par['iso_level_mm'])))
        nbad = bad_edges(f)
        volume = float(np.einsum('ij,ij->i', tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6)
        native_value = field_values(a['phi'], world, a)
        outer = (abs(native_value - par['iso_level_mm']) <= par['outer_role_field_tolerance_mm']) & (local[:, 2] > float(b['margin_z']) + par['source_comparison_above_margin_mm'])
        pts = world[outer]
        mesh = trimesh.Trimesh(a['source_vertices'], a['source_faces'], process=False)
        distance = []
        shifted = []
        for chunk in np.array_split(pts, max(1, int(np.ceil(len(pts) / 256)))):
            if len(chunk):
                distance.extend(trimesh.proximity.closest_point(mesh, chunk)[1])
                shifted.extend(trimesh.proximity.closest_point(mesh, chunk + [0, 0, 2.0])[1])
        if len(pts):
            j = np.linspace(0, len(pts) - 1, min(32, len(pts)), dtype=int)
            exact = trimesh.proximity.closest_point_naive(mesh, pts[j])[1]
            parity = float(np.max(abs(exact - np.asarray(distance)[j])))
            p95 = float(np.quantile(distance, 0.95))
            bad95 = float(np.quantile(shifted, 0.95))
        else:
            parity = p95 = bad95 = None
        gates = dict(manifold=nbad == 0, positive_volume=volume > 0, iso=residual <= m['field_iso_max_mm'], exterior=p95 is not None and p95 <= m['exterior_p95_distance_max_mm'], distance_control=parity is not None and parity <= 1e-08)
        faults = dict(missing_facet=bad_edges(f[:-1]) > 0, shifted_source=bad95 is not None and bad95 > m['exterior_p95_distance_max_mm'], shifted_field=float(np.max(abs(field_values(b['field'], world + [0, 0, 0.5], a) - par['iso_level_mm']))) > m['field_iso_max_mm'])
        row.update(gates=gates, fault_rejections=faults, outer_vertices_checked=len(pts), exterior_p95_mm=p95, independent_distance_max_error_mm=parity, field_iso_residual_max_mm=residual, nonmanifold_edges=nbad, signed_volume_mm3=volume, resolution='PER_POINT', physical_seated_gap_mm=None, physical_strength_N=None)
        rows.append(row)
    usable = [r for r in rows if r['status'] == 'EXPORTED']
    good = [r for r in usable if all(r['gates'].values())]
    out = dict(claim_type='capability', requested=len(rows), exported=len(usable), geometric_pass=len(good), failed_exported=len(usable) - len(good), unavailable=len(rows) - len(usable), exclusion_fraction=(len(rows) - len(usable)) / len(rows), rows=rows, gate=len(good) >= m['minimum_successful_pairs'] and all((all(r['fault_rejections'].values()) for r in usable)), external_referent=pr['external_referent'], scope=pr['limits'], seconds=time.perf_counter() - start)
    dump(bundle.payload / 'evaluation/CROWN_SCORES.json', out)
    return out
if __name__ == '__main__':
    build()
