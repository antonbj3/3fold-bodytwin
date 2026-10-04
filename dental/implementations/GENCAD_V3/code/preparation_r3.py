"""R3: change the morphology domain boundary; retain R2 source and failures."""
from preparation import *

def build():
    import crown_fit_geometry as G
    spec = importlib.util.spec_from_file_location('gencad_x1b', ROOT / 'vendor/crown_design_geometry.py')
    CDG = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CDG)
    CDG.G = G
    prereg = read(ROOT / 'PREREG_R3.json')
    q = prereg['parameters']
    sel = read(PAYLOAD / 'private/PREP_SELECTION.json')['payload']
    cohort = read(PAYLOAD / 'private/COHORT.json')['payload']['cases']
    lookup = {r['case_key']: r for r in cohort}
    records = []
    start = time.perf_counter()
    for s in sel:
        key = s['case_key']
        dest = PAYLOAD / 'preparations_r3' / key
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
            lo = local.min(0) - q['padding_mm']
            hi = local.max(0) + q['padding_mm']
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
    files = {str(p.relative_to(PAYLOAD)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted((PAYLOAD / 'preparations_r3').rglob('*')) if p.is_file()}
    freeze(ROOT / 'FROZEN_PREPARATIONS_R3.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R3.json'), files=files, records=records, measurement_status='No physical measurement; independent numerical tests follow freeze'))
    dump(ROOT / 'raw/PREPARATION_BUILD_COST_R3.json', dict(seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    build()
