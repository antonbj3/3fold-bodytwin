from core import *

def supports(rec):
    key = rec['key']
    tri = load(rec['private_path'])['target']
    (sv, sf) = triangles_mesh(tri)
    (tv, tf, tc) = close(sv, sf, key + '_native')
    x = load(rec['outer']['mesh_path'])
    (ev, ef) = (x['vertices'], x['faces'])
    if 'roles' in x:
        (ev, ef) = compact(ev, ef[x['roles'] == 0])
    (ov, of, oc) = close(ev, ef, key + '_outer')
    return ((tv, tf, tc), (ov, of, oc), tri)

def search(rec, mat, T, E, roundtag):
    pr = read(R / f'PREREG_{roundtag}.json')
    s = pr['search']
    rows = []
    Tdist = Distance(T[0], T[1])
    Edist = Distance(E[0], E[1])
    best = None
    for (mi, M) in enumerate(frame_candidates(rec['key'])):
        tt = T[0] @ M.T
        ee = E[0] @ M.T
        low = max(tt[:, 2].min(), ee[:, 2].min())
        high = min(tt[:, 2].max(), ee[:, 2].max())
        zspan = high - low
        center = (tt.min(0) + tt.max(0)) / 2
        for bf in s['base_fraction']:
            b = low + bf * zspan
            for ex in s['top_reduction_extra_mm']:
                H = high - mat['top_reduction_mm'] - ex
                if H - b < s['min_height_mm']:
                    continue
                zs = np.linspace(b, H, max(2, int(np.ceil((H - b) / s['axis_sampling_mm'])) + 1))
                mid = (zs[:-1] + zs[1:]) / 2
                cover = np.diff(zs) / 2
                off = offset(zs[:-1], H)
                for dx in s['center_offsets_mm']:
                    for dy in s['center_offsets_mm']:
                        c = center[:2] + [dx, dy]
                        pts = np.c_[np.tile(c, (len(mid), 1)), mid] @ M
                        ds = Tdist.query(pts)[0]
                        de = Edist.query(pts)[0]
                        ar = float(min(np.min(ds - cover - off - mat['source_reduction_mm'] - 0.025), np.min(de - cover - off - mat['wall_mm'] - G - 0.025)))
                        inside = False
                        if ar >= s['min_flat_top_radius_mm']:
                            inside = bool(Tdist.inside(pts).all() and Edist.inside(pts).all())
                        good = ar >= s['min_flat_top_radius_mm'] and inside
                        a = np.floor(ar * 1000000.0) / 1000000.0
                        volume = float(np.pi * np.trapezoid((a + offset(zs, H)) ** 2, zs)) if good else None
                        row = dict(axis_id=mi, M=M, c=c, b=b, H=H, flat_radius_mm=a, feasible=good, volume_proxy_mm3=volume, reject_reason=None if good else 'RADIUS_TOO_SMALL' if ar < s['min_flat_top_radius_mm'] else 'AXIS_OUTSIDE')
                        rows.append(row)
                        if good and (best is None or volume > best['volume_proxy_mm3']):
                            best = row
    dump(D / f"{rec['key']}_{mat['id']}_{roundtag}_SEARCH.json", rows)
    return (best, dict(requested=len(rows), rejected=sum((not x['feasible'] for x in rows)), reasons={r: sum((x['reject_reason'] == r for x in rows)) for r in ['RADIUS_TOO_SMALL', 'AXIS_OUTSIDE']}, path=str(D / f"{rec['key']}_{mat['id']}_{roundtag}_SEARCH.json")))

def construct(rec, mat, T, E, tri, ch, roundtag):
    key = rec['key'] + '_' + mat['id'] + '_' + roundtag
    M = np.array(ch['M'])
    c = np.array(ch['c'])
    a = ch['flat_radius_mm']
    b = ch['b']
    H = ch['H']
    (qv, qf, qzs, _) = body(c, a, b, H, M)
    (kv, kf, kzs, _) = body(c, a, b, H + G, M, rho=RHO + G)
    cq = certify(T[2]['path'], query_ball_rows(qv, qzs, w=mat['source_reduction_mm']), key + '_Q_native')
    cqe = certify(E[2]['path'], query_ball_rows(qv, qzs, w=mat['wall_mm'] + G), key + '_Q_outer')
    ck = certify(E[2]['path'], query_ball_rows(kv, kzs, w=mat['wall_mm']), key + '_K_outer')
    if not all((x['all_pass'] for x in [cq, cqe, ck])):
        raise ValueError('EXACT_COVER_FAILURE')
    box = trimesh.creation.box(extents=[100, 100, 100])
    box.vertices += np.array([c[0], c[1], b - 50])
    box.vertices = box.vertices @ M
    bp = D / (key + '_cutbox.mesh')
    meshwrite(bp, box.vertices, box.faces)
    trim = D / (key + '_trimmed.mesh')
    clipinfo = boolean(E[2]['path'], bp, trim)
    (ev, ef) = meshread(trim)
    (cv0, cf0, _, _) = body(c, a, b, H + G, M, rho=RHO + G, bottom_override=b - 10)
    cav = D / (key + '_cutter.mesh')
    meshwrite(cav, cv0, cf0)
    crown = D / (key + '_crown.mesh')
    bo = boolean(trim, cav, crown)
    (v, f) = meshread(crown)
    pts = v[f].mean(1)
    ld = pts @ M.T
    innerd = Distance(cv0, cf0).query(pts)[0]
    roles = np.where(abs(ld[:, 2] - b) < 1e-07, 2, np.where(innerd < 1e-07, 1, 0))
    path = D / (key + '.npz')
    np.savez_compressed(path, vertices=v, faces=f, roles=roles, prep_vertices=qv, prep_faces=qf, cavity_vertices=kv, cavity_faces=kf, target=tri, M=M, base=b, top=H, center=c, flat_radius=a)
    qpath = D / (key + '_prep.stl')
    trimesh.Trimesh(qv, qf, process=False).export(qpath)
    return dict(status='GENERATED', mesh_path=str(path), mesh_sha256=sha(path), prep_stl=str(qpath), prep_stl_sha256=sha(qpath), construction=ch, certificates=dict(Q_dilated_inside_native=cq, Q_plus_wall_gap_inside_outer=cqe, cavity_wall_inside_outer=ck), boolean=bo, clip=clipinfo, source_closure=T[2], outer_closure=E[2], prediction=dict(wall_min_mm=mat['wall_mm'], internal_gap_mm=G, gap_max_error_mm=0.01, isolated_insertion='PASS_BY_MONOTONE_CONSTRUCTION', shape='UNMEASURED', material_shoulder='UNMEASURED'), scope='Q is compact virtual core; full tooth shoulder and whole-arch suitability not inferred')

def run(tag='R1', limit=None):
    st = time.perf_counter()
    pr = read(R / f'PREREG_{tag}.json')
    records = read(R / 'INPUTS.json')['records']
    records = records if limit is None else records[:limit]
    rows = []
    for rec in records:
        pilot = R / f'FROZEN_PREDICTIONS_{tag}_PILOT.json'
        if limit is None and pilot.exists():
            reused = [x for x in read(pilot)['rows'] if x['key'] == rec['key']]
            if len(reused) == len(pr['materials']):
                rows.extend(reused)
                print('REUSE_FROZEN_PILOT', rec['key'], flush=True)
                continue
        try:
            (T, E, tri) = supports(rec)
            err = None
        except Exception as e:
            err = repr(e)
        for mat in pr['materials']:
            tick = time.perf_counter()
            r = dict(key=rec['key'], case_key=rec['case_key'], family=rec['family'], fdi=rec['fdi'], material=mat['id'], material_contract=mat, round=tag, outer_route=rec['outer'], resolution='PER_TOOTH')
            try:
                if err:
                    raise ValueError(err)
                (ch, sr) = search(rec, mat, T, E, tag)
                r['search'] = sr
                if ch is None:
                    raise ValueError('NO_FEASIBLE_COMPACT_TAPER_IN_FROZEN_FAMILY')
                r.update(construct(rec, mat, T, E, tri, ch, tag))
            except Exception as e:
                r.update(status='REJECTED', reason=repr(e))
            r['generation_seconds'] = time.perf_counter() - tick
            rows.append(r)
            dump(R / f'raw/{tag}_GENERATION.json', rows)
            state(tag + '_GENERATING', dict(key=rec['key'], material=mat['id'], status=r['status'], reason=r.get('reason')), 'Freeze outputs, independently evaluate all chain links')
            print(rec['key'], mat['id'], r['status'], r.get('reason', '')[:100], flush=True)
    name = f'FROZEN_PREDICTIONS_{tag}' + ('_PILOT' if limit else '') + '.json'
    lock(name, dict(prereg_sha256=sha(R / f'PREREG_{tag}.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, physical_predictions=None, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'R1', 3 if '--pilot' in sys.argv else None)
