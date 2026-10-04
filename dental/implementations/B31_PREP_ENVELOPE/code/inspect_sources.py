from geometry import *

def main():
    rows = []
    tic = time.perf_counter()
    for r in read(ROOT / 'FROZEN_COHORT.json')['records']:
        assert sha(r['private_path']) == r['private_sha256']
        assert sha(r['public_path']) == r['public_sha256']
        t = load_np(r['private_path'])['target']
        p = load_np(r['public_path'])
        base = np.floor((p['margin_curve'][:, 2].min() - 0.2) / 0.2) * 0.2
        rr = {'key': r['key'], 'family': r['family'], 'base_mm': base, 'source_triangles': len(t), 'native_margin_z_span_mm': np.ptp(p['margin_curve'][:, 2]), 'native_top_z_mm': t[:, :, 2].max(), 'resolution': 'PER_TOOTH'}
        try:
            (m, rim, poly) = closed_source(t, base)
            rr.update(status='VIRTUAL_SOURCE_CLOSED', projected_margin_area_mm2=poly.area, virtual_volume_mm3=m.volume, closure_watertight=bool(m.is_watertight), closure_self_intersection='UNKNOWN_NOT_YET_CHECKED', rim_points=len(rim))
            path = DATA / 'source_closures' / (r['key'] + '.npz')
            path.parent.mkdir(exist_ok=True)
            np.savez_compressed(path, vertices=m.vertices, faces=m.faces, original=t, rim=rim, base=base)
            rr.update(path=path, sha256=sha(path))
        except Exception as e:
            rr.update(status='UNKNOWN_SOURCE_CLOSURE', reason=str(e))
        rows.append(rr)
        dump(ROOT / 'raw/SOURCE_INSPECTION.json', rows)
        print(rr['key'], rr['status'], rr.get('reason', ''), flush=True)
    dump(ROOT / 'raw/SOURCE_INSPECTION_COST.json', {'seconds': time.perf_counter() - tic, **usage()})
    state('SOURCE_INSPECTED', {'closed': sum((r['status'] == 'VIRTUAL_SOURCE_CLOSED' for r in rows)), 'requested': 18}, 'Build eroded full-box field and D carrier dictionary; unavailable closures remain UNKNOWN')
if __name__ == '__main__':
    main()
