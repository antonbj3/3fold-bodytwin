from construct_c import *

def continuous_cut(v, f, values):
    verts = v.tolist()
    out = []
    cache = {}

    def intersect(a, b):
        key = tuple(sorted((int(a), int(b))))
        if key not in cache:
            (aa, bb) = key
            t = values[aa] / (values[aa] - values[bb])
            cache[key] = len(verts)
            verts.append(((1 - t) * v[aa] + t * v[bb]).tolist())
        return cache[key]
    for face in f:
        poly = []
        for (a, b) in zip(face, np.roll(face, -1)):
            if values[a] >= 0:
                poly.append(int(a))
            if (values[a] >= 0) != (values[b] >= 0):
                poly.append(intersect(a, b))
        if len(poly) >= 3:
            out.extend([[poly[0], poly[k], poly[k + 1]] for k in range(1, len(poly) - 1)])
    return compact(np.array(verts), np.array(out, int))

def make_e(m, key):
    cache = D / 'C_IMPLICIT_CACHE' / (key + '.npz')
    if not cache.exists():
        try:
            implicit_inner(m, key + '_E')
        except ValueError:
            pass
        cache = D / 'C_IMPLICIT_CACHE' / (key + '_E.npz')
    data = npz(cache)
    (v, f, pv, pf) = (data['vertices'], data['faces'], data['prep_vertices'], data['prep_faces'])
    nt = int(data['original_cap_faces'])
    cap = pv[pf[:nt]]
    extension = pv[pf[nt:]]
    (_, dc, _) = fast_nearest(cap, v)
    (_, de, _) = fast_nearest(extension, v)
    (iv, inf) = continuous_cut(v, f, de - dc - 1e-06)
    obj = trimesh.Trimesh(iv, inf, process=False)
    ll = loops(obj)
    if len(ll) != 1:
        raise ValueError('continuous cap has ' + str(len(ll)) + ' boundary loops')
    if len(obj.split(only_watertight=False)) != 1:
        raise ValueError('continuous cap has multiple components')
    for _ in range(8):
        (sd, idx, q, n) = igl.signed_distance(iv, pv, pf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        iv += (boundary_gap(q, data['rim']) - sd)[:, None] * n
    tri = iv[inf]
    (sd, idx, q, n) = igl.signed_distance(tri.mean(1), pv, pf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    res = abs(sd - boundary_gap(q, data['rim']))
    normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    if np.einsum('ij,ij->i', normal, n).sum() > 0:
        inf = inf[:, ::-1]
    (_, cover, _) = fast_nearest(iv[inf], sample(cap, 8192))
    return (iv, inf, dict(max_plane_offset_residual_mm=float(res.max()), metric='centroid implicit field residual', preparation_coverage_sampled_max_mm=float(cover.max()), preparation_coverage_scope='8192 source cap probes; no continuous coverage proof', boundary_loop_count=1, cache_path=cache, cache_sha256=sha(cache)))

def run():
    st = time.perf_counter()
    rows = []
    for (rec, rr) in inputs():
        m = npz(rec['mesh_path'])
        t = time.perf_counter()
        try:
            (iv, inf, off) = make_e(m, rec['key'])
            error = None
        except Exception as exc:
            error = repr(exc)
        for method in ['continuous_offset_only', 'continuous_local_thickening']:
            row = dict(key=rec['key'], family=rec['family'], method=method, status='REJECTED', resolution='PER_POINT')
            if error:
                row['reason'] = error
            else:
                try:
                    (ov, of) = compact(m['vertices'], m['faces'][m['roles'] == 0])
                    info = {}
                    if method == 'continuous_local_thickening':
                        (ov, of, info) = outer_repair(m, iv[inf])
                    (v, f, roles) = join(ov, of, iv, inf)
                    dest = D / 'E' / method / (rec['key'] + '.npz')
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(dest, vertices=v, faces=f, roles=roles, prep_triangles=m['vertices'][m['faces'][m['roles'] == 1]])
                    row.update(status='GENERATED', mesh_path=dest, mesh_sha256=sha(dest), offset=off, repair=info)
                except Exception as exc:
                    row['reason'] = repr(exc)
            row['construction_seconds'] = time.perf_counter() - t
            rows.append(row)
            save(R / 'raw/E_GENERATION.json', rows)
            print(rec['key'], method, row['status'], row.get('reason', off if not error else ''), flush=True)
        checkpoint('E_GENERATION', len(rows), 'Continuous margin-side extraction on all same18 source caps')
    freeze(R / 'FROZEN_PREDICTIONS_E.json', dict(claim_type='capability', prereg_sha256=sha(R / 'PREREG_E.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
