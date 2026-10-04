from common import *
from contact import compare, gates, enclosure, spatial, Q
from reuse import VF, grid
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import collections, resource

def proximal(tri, ext):
    """Connected near patches on the measured neighbour mesh, PL vertex distance."""
    if not len(tri):
        return dict(status='UNKNOWN_NO_NEIGHBOUR')
    lo = ext.min((0, 1)) - 0.2
    hi = ext.max((0, 1)) + 0.2
    keep = np.all(tri.max(1) >= lo, axis=1) & np.all(tri.min(1) <= hi, axis=1)
    tt = tri[keep]
    if not len(tt):
        return dict(status='NO_NEAR_CANDIDATE', count=0, area_mm2=0.0, regions=[], cropped_neighbour_faces=0, total_neighbour_faces=len(tri))
    (v, inv) = np.unique(tt.reshape(-1, 3), axis=0, return_inverse=True)
    f = inv.reshape(-1, 3)
    d = VF.distance(VF.make_mesh(ext), v)
    areas = []
    mom = []
    idx = []
    links = []
    edges = {}
    for (fi, ids) in enumerate(f):
        if d[ids].min() > 0.1:
            continue
        poly = Q.clip(np.eye(3), 0.1 - d[ids])
        pts = poly @ v[ids]
        if len(pts) < 3:
            continue
        a = np.linalg.norm(np.cross(pts[1:-1] - pts[0], pts[2:] - pts[0]), axis=1) / 2
        ar = float(a.sum())
        if ar <= 1e-12:
            continue
        j = len(idx)
        idx.append(fi)
        areas.append(ar)
        mom.append(np.sum(a[:, None] * (pts[0] + pts[1:-1] + pts[2:]) / 3, axis=0))
        for (x, y) in zip(ids, np.roll(ids, -1)):
            if min(d[x], d[y]) < 0.1 or d[x] == d[y] == 0.1:
                k = tuple(sorted((int(x), int(y))))
                if k in edges:
                    links.append((j, edges[k]))
                else:
                    edges[k] = j
    regs = []
    if idx:
        ij = np.array(links, int).reshape(-1, 2)
        (n, lab) = connected_components(coo_matrix((np.ones(len(ij)), (ij[:, 0], ij[:, 1])), shape=(len(idx), len(idx))).tocsr(), directed=False)
        aa = np.bincount(lab, weights=areas)
        mm = np.array([np.bincount(lab, weights=np.array(mom)[:, k]) for k in range(3)]).T
        regs = [dict(area_mm2=float(ar), centroid_xyz_mm=mm[k] / ar) for (k, ar) in enumerate(aa)]
    return dict(status='UNSIGNED_PL_PROXIMITY', count=len(regs), area_mm2=float(sum(areas)), regions=regs, cropped_neighbour_faces=len(tt), total_neighbour_faces=len(tri), vertex_distance_min_mm=float(d.min()), physical_contact='UNKNOWN', distance_interpolation_enclosure='MISSING')

def run():
    st = time.perf_counter()
    pr = read(ROOT / 'PREREG_R2.json')
    metrics = pr['metrics']
    meta = {r['key']: r for r in read(V4 / 'payload/whole_inputs/RECORDS.json')}
    cache = {}
    out = []
    ctrl = []
    DATA.mkdir(parents=True, exist_ok=True)
    for (n, r0) in enumerate(pr['cohort']):
        r = dict(r0)
        key = r['key']
        md = meta[key]
        r.update(dataset=md['dataset'], source_fdi=md['source_fdi'], resolution='PER_TOOTH')
        if r['status'] != 'AVAILABLE':
            r['L2'] = 'FAIL'
            out.append(r)
            continue
        if key not in cache:
            a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
            p = npz(V4 / 'payload/whole_inputs' / key / 'preparation.npz')
            tri = a['source_triangles']
            tri = tri[tri.mean(1)[:, 2] >= p['margin_z']]
            (xy, ff) = grid(tri)
            ceiling = VF.height(a['antagonist_triangles'], xy, True)
            rg = ceiling - VF.height(tri, xy)
            cache[key] = (a, xy, ff, ceiling, rg, tri)
        (a, xy, ff, ceiling, rg, tri) = cache[key]
        m = npz(r['mesh_path'])
        ext = m['vertices'][m['faces'][m['face_roles'] == 0]]
        if not len(ext):
            r.update(status='NO_EXTERIOR', L2='FAIL')
            out.append(r)
            continue
        pg = ceiling - VF.height(ext, xy)
        c = compare(xy, ff, pg, rg)
        gg = gates(c, metrics)
        b = enclosure(xy, ff, pg, rg, metrics['gap_uncertainty_mm'])
        robust = 'PASS' if b.get('symdiff_upper_mm2', float('inf')) <= 1 + 1e-09 else 'FAIL' if b.get('symdiff_lower_mm2', 0) > 1 + 1e-09 else 'UNKNOWN'
        prox = {side: proximal(a[side + '_triangles'], ext) for side in ['mesial', 'distal']}
        r.update(status='SCORED', contact=c, gates=gg, pattern_enclosure=b, robust_pattern=robust, proximal=prox, L2='PASS' if all(gg.values()) and r['parent_geometry_pass'] else 'UNKNOWN' if c['status'] != 'SCORED' else 'FAIL')
        artifact = DATA / 'contact_fields' / (r['uid'] + '.npz')
        artifact.parent.mkdir(exist_ok=True, parents=True)
        np.savez_compressed(artifact, xy=xy, faces=ff, pred_gap=pg, reference_gap=rg, ceiling=ceiling)
        r.update(field_path=str(artifact), field_sha256=sha(artifact))
        out.append(r)
        if n < 3 or (key == pr['cohort'][0]['key'] and r['participant'] in FRONTIER and (r['track'] == 'V5B_R3')):
            q = Q.contact_map(xy, ff, pg, rg)
            cc = dict(uid=r['uid'], area_parity_mm2=abs(q['contact_symdiff_mm2'] - c['symdiff_mm2']) if q['contact_symdiff_mm2'] is not None else None, self_identity=compare(xy, ff, rg, rg).get('symdiff_mm2'))
            if r['participant'] in FRONTIER:
                (x2, f2) = grid(tri, 0.125)
                c2 = VF.height(a['antagonist_triangles'], x2, True)
                m2 = compare(x2, f2, c2 - VF.height(ext, x2), c2 - VF.height(tri, x2))
                cc.update(fine_grid_contact=m2, fine_L2=all(gates(m2, metrics).values()) and r['parent_geometry_pass'], coarse_L2=r['L2'] == 'PASS')
            ctrl.append(cc)
        if n % 12 == 0:
            dump(ROOT / 'raw/R2_CHECKPOINT.json', out)
            state('R2_RUNNING', str(n + 1) + ' complete crown queries', 'Remaining geometry, uncertainty and neighbour patches')
            print('R2', n + 1, len(pr['cohort']), round(time.perf_counter() - st, 1), flush=True)
    by = collections.defaultdict(list)
    for r in out:
        if r['track'] == 'V5B_R3' and r['participant'] in FRONTIER:
            by[r['key']].append(r)
    items = [dict(key=k, participants={r['participant']: r['L2'] for r in rr}, discriminates=len({r['L2'] == 'PASS' for r in rr}) > 1, all_scored=all((r['status'] == 'SCORED' for r in rr))) for (k, rr) in by.items()]
    result = dict(requested=len(out), scored=sum((r['status'] == 'SCORED' for r in out)), L2_counts=dict(collections.Counter((r['L2'] for r in out))), robust_pattern_counts=dict(collections.Counter((r.get('robust_pattern', 'NOT_SCORED') for r in out))), frontier_items=items, frontier_discriminating=sum((x['discriminates'] for x in items)), dropout=dict(reasons=dict(collections.Counter((r['status'] for r in out if r['status'] != 'SCORED'))), fraction=sum((r['status'] != 'SCORED' for r in out)) / len(out)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'raw/R2_ROWS.json', out)
    dump(ROOT / 'raw/R2_CONTROLS.json', ctrl)
    dump(ROOT / 'rounds/R2.json', result)
    state('R2_DECIDED', str(result['frontier_discriminating']) + '/' + str(len(items)) + ' complete-crown sites discriminate', 'Replace source tooth in full arch contact graph and evaluate load sets')
    (ROOT / 'history/HANDOFF_R2.md').write_text(json.dumps(clean(result), indent=2) + '\nNext: full arch load-set consumer. No physical pressure validation achieved.\n')
    return result
if __name__ == '__main__':
    print(clean(run()))
