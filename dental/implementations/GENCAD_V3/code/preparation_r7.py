"""R7 separates virtual closure from the source-supported crown region."""
from source_solid import *

def support_metrics(a, guard, above):
    v = a['closed_source_vertices']
    f = a['closed_source_faces']
    n = len(a['source_faces'])
    caps = v[f[n:]]
    if not len(caps):
        raise ValueError('no explicit virtual cap partition')
    capmax = float(caps[:, :, 2].max())
    margin = capmax + guard
    tri = (a['source_vertices'][a['source_faces']] - a['center']) @ a['R']
    area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    retained = np.all(tri[:, :, 2] > margin + above, axis=1)
    return dict(virtual_margin_z_mm=margin, virtual_cap_max_z_mm=capmax, closure_separation_mm=margin - capmax, retained_source_area_fraction=float(area[retained].sum() / area.sum()), retained_source_triangles=int(retained.sum()), source_triangles=len(tri), source_area_mm2=float(area.sum()), retained_source_area_mm2=float(area[retained].sum()))

def build():
    import crown_fit_geometry as G
    spec = importlib.util.spec_from_file_location('gencad_x1b_r7', ROOT / 'vendor/crown_design_geometry.py')
    CDG = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CDG)
    CDG.G = G
    pr = read(ROOT / 'PREREG_R7.json')
    q = pr['parameters']
    prior = read(ROOT / 'FROZEN_CROWNS_R6.json')['payload']
    meta = {r['case_key']: r for r in prior['records']}
    records = []
    for key in pr['cases']:
        start = time.perf_counter()
        record = dict(case_key=key, split=meta[key]['split'], status='REFUSED_SOURCE_COVERAGE')
        path = PAYLOAD / 'preparations_r6' / key / 'fields.npz'
        if sha(path) != prior['files']['preparations_r6/' + key + '/fields.npz']['sha256']:
            raise ValueError('R6 input drift')
        a = dict(np.load(path, allow_pickle=False))
        support = support_metrics(a, q['closure_guard_mm'], q['source_comparison_above_margin_mm'])
        record.update(support)
        if support['retained_source_area_fraction'] < pr['metrics']['minimum_retained_source_area_fraction']:
            record['reason'] = 'insufficient measured exterior above every virtual closure'
            records.append(record)
            print(key, record['status'], support['retained_source_area_fraction'], flush=True)
            continue
        h = float(a['step'])
        grid = G.Grid(a['origin'], a['origin'] + h * (np.array(a['phi'].shape) - 1), h)
        margin = support['virtual_margin_z_mm']
        case = dict(grid=grid, phi_T=a['phi'], z_m=margin)
        dest = PAYLOAD / 'preparations_r7' / key
        dest.mkdir(parents=True, exist_ok=True)
        scratch = dest / 'construction'
        scratch.mkdir(exist_ok=True)
        CDG.SCRATCH = str(scratch)
        (prep, mp) = CDG.preparation(case, q['occlusal_reduction_mm'], q['axial_reduction_mm'], q['spacer_mm'], cache=False)
        (pv, pf) = G.surface_mesh(prep, grid)
        if np.einsum('ij,ij->i', pv[pf[:, 0]], np.cross(pv[pf[:, 1]], pv[pf[:, 2]])).sum() < 0:
            pf = pf[:, ::-1]
        export_stl(dest / 'preparation.stl', pv @ a['R'].T + a['center'], pf)
        a.update(preparation=prep, preparation_vertices=pv, preparation_faces=pf)
        np.savez_compressed(dest / 'fields.npz', **a)
        for p in scratch.glob('*.npz'):
            p.unlink()
        scratch.rmdir()
        cavity = G.design_cavity(prep, grid, np.full(prep.shape, q['spacer_mm']), levels=(0.0,))
        z = grid.axes()[2][None, None, :]
        field = np.maximum.reduce([a['phi'], -cavity, np.broadcast_to(margin - z, prep.shape)]).astype(np.float32)
        (cv, cf, _, _) = marching_cubes(field, q['iso_level_mm'], spacing=(h,) * 3, allow_degenerate=False)
        cv += grid.origin
        if np.einsum('ij,ij->i', cv[cf[:, 0]], np.cross(cv[cf[:, 1]], cv[cf[:, 2]])).sum() < 0:
            cf = cf[:, ::-1]
        cd = PAYLOAD / 'crowns_r7' / key
        cd.mkdir(parents=True, exist_ok=True)
        export_stl(cd / 'crown.stl', cv @ a['R'].T + a['center'], cf)
        np.savez_compressed(cd / 'prediction.npz', field=field, cavity=cavity, vertices=cv, faces=cf, margin_z=margin)
        record.update(status='EXPORTED', seconds=time.perf_counter() - start)
        records.append(record)
        budget()
    files = {str(p.relative_to(PAYLOAD)): {'sha256': sha(p), 'bytes': p.stat().st_size} for d in ['preparations_r7', 'crowns_r7'] for p in sorted((PAYLOAD / d).rglob('*')) if p.is_file()}
    freeze(ROOT / 'FROZEN_CROWNS_R7.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R7.json'), source_freeze_sha256=sha(ROOT / 'FROZEN_CROWNS_R6.json'), files=files, records=records, measurement_status='Eligibility evaluated under frozen rule; crown comparison after prediction freeze'))
if __name__ == '__main__':
    build()
