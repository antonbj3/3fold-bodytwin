"""R6 preserves triangles while splitting disconnected vertex fans."""
from source_solid import *
from topology import split_fans
from source_solid import cap_patch as original_cap

def cap_patch(v, f, apex_z):
    (vv, ff, diagnostic) = split_fans(v, f)
    return original_cap(vv, ff, apex_z)

def build():
    import crown_fit_geometry as G
    spec = importlib.util.spec_from_file_location('gencad_x1b_r5', ROOT / 'vendor/crown_design_geometry.py')
    CDG = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CDG)
    CDG.G = G
    pr = read(ROOT / 'PREREG_R6.json')
    q = pr['parameters']
    prior = read(ROOT / 'FROZEN_PREPARATIONS_R3.json')['payload']
    meta = {r['case_key']: r for r in prior['records']}
    records = []
    for key in pr['cases']:
        start = time.perf_counter()
        record = dict(case_key=key, split=meta[key]['split'], status='FAILED_SOURCE_SOLID')
        try:
            source_path = PAYLOAD / 'preparations_r3' / key / 'fields.npz'
            if sha(source_path) != prior['files']['preparations_r3/' + key + '/fields.npz']['sha256']:
                raise ValueError('source drift')
            a = dict(np.load(source_path, allow_pickle=False))
            local = (a['source_vertices'] - a['center']) @ a['R']
            base = float(local[:, 2].min())
            (v, f, loops) = cap_patch(local, a['source_faces'], base - q['basal_apex_below_source_mm'])
            mesh = trimesh.Trimesh(v, f, process=False)
            h = q['grid_spacing_mm']
            grid = G.Grid(v.min(0) - q['padding_mm'], v.max(0) + q['padding_mm'], h)
            if np.prod(grid.shape) > q['max_grid_nodes']:
                raise ValueError('bounded grid node limit')
            occupied = np.zeros(grid.shape, bool)
            (xs, ys, zs) = grid.axes()
            (xx, yy) = np.meshgrid(xs, ys, indexing='ij')
            xy = np.c_[xx.ravel(), yy.ravel()]
            for (k, z) in enumerate(zs):
                pts = np.c_[xy, np.full(len(xy), z)]
                vals = []
                for block in np.array_split(pts, max(1, int(np.ceil(len(pts) / 512)))):
                    vals.extend(mesh.contains(block))
                occupied[:, :, k] = np.asarray(vals).reshape(xx.shape)
            if not occupied.any():
                raise ValueError('empty capped-mesh volume')
            phi = (ndimage.distance_transform_edt(~occupied, sampling=h) - ndimage.distance_transform_edt(occupied, sampling=h)).astype(np.float32)
            margin = base + q['virtual_margin_above_min_source_mm']
            case = dict(grid=grid, phi_T=phi, z_m=margin)
            dest = PAYLOAD / 'preparations_r6' / key
            dest.mkdir(parents=True, exist_ok=True)
            scratch = dest / 'construction'
            scratch.mkdir(exist_ok=True)
            CDG.SCRATCH = str(scratch)
            (prep, mp) = CDG.preparation(case, q['occlusal_reduction_mm'], q['axial_reduction_mm'], q['spacer_mm'], cache=False)
            (pv, pf) = G.surface_mesh(prep, grid)
            if np.einsum('ij,ij->i', pv[pf[:, 0]], np.cross(pv[pf[:, 1]], pv[pf[:, 2]])).sum() < 0:
                pf = pf[:, ::-1]
            export_stl(dest / 'preparation.stl', pv @ a['R'].T + a['center'], pf)
            export_stl(dest / 'native_exterior.stl', a['source_vertices'], a['source_faces'])
            np.savez_compressed(dest / 'fields.npz', source_vertices=a['source_vertices'], source_faces=a['source_faces'], phi=phi, preparation=prep, origin=grid.origin, step=h, R=a['R'], center=a['center'], preparation_vertices=pv, preparation_faces=pf, closed_source_vertices=v, closed_source_faces=f)
            for p in scratch.glob('*.npz'):
                p.unlink()
            scratch.rmdir()
            cavity = G.design_cavity(prep, grid, np.full(prep.shape, q['spacer_mm']), levels=(0.0,))
            z = grid.axes()[2][None, None, :]
            field = np.maximum.reduce([phi, -cavity, np.broadcast_to(margin - z, phi.shape)]).astype(np.float32)
            (cv, cf, _, _) = marching_cubes(field, q['iso_level_mm'], spacing=(h,) * 3, allow_degenerate=False)
            cv += grid.origin
            if np.einsum('ij,ij->i', cv[cf[:, 0]], np.cross(cv[cf[:, 1]], cv[cf[:, 2]])).sum() < 0:
                cf = cf[:, ::-1]
            cd = PAYLOAD / 'crowns_r6' / key
            cd.mkdir(parents=True, exist_ok=True)
            export_stl(cd / 'crown.stl', cv @ a['R'].T + a['center'], cf)
            np.savez_compressed(cd / 'prediction.npz', field=field, cavity=cavity, vertices=cv, faces=cf, margin_z=margin)
            record.update(status='EXPORTED', boundary_loops=len(loops), source_faces=len(a['source_faces']), closure_faces=len(f) - len(a['source_faces']), grid_shape=grid.shape, virtual_margin_z_mm=margin)
        except (ValueError, RuntimeError, IndexError) as e:
            record.update(reason=str(e))
        record['seconds'] = time.perf_counter() - start
        records.append(record)
        print('source-solid', key, record['status'], record.get('reason', ''), round(record['seconds'], 2), flush=True)
        budget()
    files = {str(p.relative_to(PAYLOAD)): {'sha256': sha(p), 'bytes': p.stat().st_size} for d in ['preparations_r6', 'crowns_r6'] for p in sorted((PAYLOAD / d).rglob('*')) if p.is_file()}
    freeze(ROOT / 'FROZEN_CROWNS_R6.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R6.json'), files=files, records=records, measurement_status='Source measurement held until frozen construction'))
if __name__ == '__main__':
    build()
