from construct_b import *
sys.path.insert(0, str(BASE / 'LANE_X49_DESIGN_GATE/code'))
import design_gate as x49
from export import three_mf, load3mf

def run(tag):
    st = time.perf_counter()
    fr = read(R / f'FROZEN_PREDICTIONS_{tag}.json')
    idx = {r['key']: (rec, r) for (rec, r) in inputs()}
    rows = []
    for rec in fr['rows']:
        if 'local_thickening' not in rec['method'] or 'mesh_path' not in rec:
            continue
        m = npz(rec['mesh_path'])
        p = npz(idx[rec['key']][1]['public_path'])
        folder = D / 'x49' / tag / rec['key']
        folder.mkdir(parents=True, exist_ok=True)
        meshes = {'crown': trimesh.Trimesh(m['vertices'], m['faces'], process=False), 'prep': mesh(m['prep_triangles']), 'antagonist': mesh(p['antagonist'])}
        paths = {}
        for (k, v) in meshes.items():
            path = folder / (k + '.stl')
            v.export(path)
            paths[k] = str(path)
        contract = dict(units='mm', common_frame_confirmed=False, physical_geometry_status='virtual preparation and raw unverified bite; digital nominal geometry only', ifu_profile='katana-ht', indication='anterior' if rec['family'] == 'anterior' else 'posterior', input_sha256={k: sha(v) for (k, v) in paths.items()}, surface_error_mm={k: 0.0 for k in paths}, error_source='zero at encoded nominal STL only; scanner/process uncertainty UNKNOWN', regions={'exterior': np.flatnonzero(m['roles'] == 0), 'intaglio': np.flatnonzero(m['roles'] == 1), 'preparation': np.arange(len(m['prep_triangles']))})
        save(folder / 'contract.json', contract)
        try:
            report = x49.check(paths, contract=str(folder / 'contract.json'), units='mm')
            error = None
        except Exception as exc:
            report = None
            error = repr(exc)
        save(folder / 'REPORT.json', dict(report=report, error=error))
        rr = dict(key=rec['key'], family=rec['family'], method=rec['method'], report_path=folder / 'REPORT.json', sha256=sha(folder / 'REPORT.json'), error=error, rules={k: v['status'] for (k, v) in report['rules'].items()} if report else {})
        rows.append(rr)
        save(R / f'raw/X49_{tag}.json', rows)
        print(tag, rec['key'], rr['rules'], error, flush=True)
    save(R / f'X49_{tag}.json', dict(rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, scope='Actual unchanged X49. Inherited0.5mm wall independently retained. Physical pose unknown; metadata does not fabricate registration.'))
if __name__ == '__main__':
    run(sys.argv[1])
