from ray_crown import *

def select_oriented(v, f, native, key):
    pr = read(R / 'PREREG_D.json')
    S = v[f[native]].copy()
    rows = []
    for tilt in pr['search']['tilt_degrees']:
        for az in [0] if tilt == 0 else pr['search']['azimuth_degrees']:
            (t, a) = np.deg2rad([tilt, az])
            d = np.array([np.sin(t) * np.cos(a), np.sin(t) * np.sin(a), np.cos(t)])
            e = np.array([np.cos(t) * np.cos(a), np.cos(t) * np.sin(a), -np.sin(t)])
            b = np.cross(d, e)
            M = np.stack([e, b, d])
            vv = v @ M.T
            ss = S @ M.T
            zmin = float(vv[:, 2].min())
            zmax = float(vv[:, 2].max())
            for rel in pr['search']['top_z_relative_max_mm']:
                H = zmax + rel
                zs = np.linspace(zmin, H, max(2, int(np.ceil((H - zmin) / 0.05)) + 1))
                cover = float(np.diff(zs).max() / 2)
                for x in pr['search']['center_offsets_mm']:
                    for y in pr['search']['center_offsets_mm']:
                        pts = np.c_[np.full(len(zs), x), np.full(len(zs), y), zs]
                        dd = fast_nearest(ss, pts)[1]
                        lb = float(dd.min() - cover)
                        r = np.floor((lb - 0.65) * 1000000.0) / 1000000.0
                        if r >= 0.5:
                            wn = float(igl.winding_number(v, f, np.array([[x, y, H]]) @ M)[0])
                        else:
                            wn = None
                        good = r >= 0.5 and abs(wn) > 0.5
                        rows.append(dict(tilt_degrees=tilt, azimuth_degrees=az, world_to_local=M.tolist(), center_xy=[x, y], axis_top_mm=H, base_mm=zmin - 2, radius_mm=float(r), axis_clearance_lower_mm=lb, axis_sample_min_mm=float(dd.min()), axis_cover_mm=cover, inside_top_winding=wn, feasible=bool(good), volume_proxy_mm3=float(np.pi * r * r * (H - zmin) + 2 * np.pi * r ** 3 / 3) if good else None, reject_reason=None if good else 'RADIUS_LT0P5' if r < 0.5 else 'AXIS_TOP_OUTSIDE'))
    dump(R / 'raw' / ('D_SEARCH_' + key + '.json'), rows)
    good = [q for q in rows if q['feasible']]
    if not good:
        raise ValueError('NO_FEASIBLE_ORIENTED_RAY_IN3300_CANDIDATES')
    return (max(good, key=lambda q: q['volume_proxy_mm3']), rows)

def generate_d(rec):
    st = time.perf_counter()
    (v, f, native, ci) = close(rec)
    (ch, search) = select_oriented(v, f, native, rec['key'])
    c = np.array(ch['center_xy'])
    H = ch['axis_top_mm']
    r = ch['radius_mm']
    M = np.array(ch['world_to_local'])
    (qv, qf) = capsule(c, H, r, ch['base_mm'])
    (kv, kf) = capsule(c, H, r + 0.05, ch['base_mm'])
    qv = qv @ M
    kv = kv @ M
    support = D / (rec['key'] + '_D_support.mesh')
    cavity = D / (rec['key'] + '_D_cavity.mesh')
    out = D / (rec['key'] + '_D_crown.mesh')
    meshwrite(support, v, f)
    meshwrite(cavity, kv, kf)
    proc = subprocess.run([str(D / 'boolean'), str(support), str(cavity), str(out)], capture_output=True, text=True)
    if proc.returncode:
        raise ValueError('BOOLEAN_FAILURE ' + str(proc.returncode) + ' ' + proc.stdout[:200] + ' ' + proc.stderr[:200])
    bo = json.loads(proc.stdout)
    (cv, cf) = meshread(out)
    S = v[f[native]].copy()
    A = v[f[~native]].copy()
    src = {tuple(sorted((tuple(x) for x in t))) for t in S}
    kept = np.array([tuple(sorted((tuple(x) for x in t))) in src for t in cv[cf]])
    roles = np.zeros(len(cf), int)
    pts = cv[cf[~kept]].mean(1)
    roles[~kept] = np.where(fast_nearest(kv[kf], pts)[1] < fast_nearest(A, pts)[1], 1, 2)
    path = D / (rec['key'] + '_D.npz')
    np.savez_compressed(path, vertices=cv, faces=cf, roles=roles, prep_vertices=qv, prep_faces=qf, cavity_vertices=kv, cavity_faces=kf, native_triangles=S, closure_triangles=A, support_vertices=v, support_faces=f)
    result = dict(key=rec['key'], family=rec['family'], status='GENERATED', round='D', mesh_path=str(path), mesh_sha256=sha(path), chosen=ch, search_requested=len(search), search_rejected=sum((not q['feasible'] for q in search)), native_facets=len(S), retained_native_facets=int(kept.sum()), boolean=bo, closure=ci, wall_analytic_lower_mm=ch['axis_clearance_lower_mm'] - r - 0.05, seconds=time.perf_counter() - st)
    print('D', rec['key'], 'r', r, 'tilt', ch['tilt_degrees'], 'native', int(kept.sum()), len(S), 'seconds', result['seconds'], flush=True)
    return result

def run():
    records = read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']
    rows = []
    st = time.perf_counter()
    if '--first' in sys.argv:
        records = records[:3]
    for rec in records:
        try:
            out = generate_d(rec)
        except Exception as e:
            out = dict(key=rec['key'], family=rec['family'], status='REJECTED', reason=repr(e))
            print(out, flush=True)
        rows.append(out)
        dump(R / 'raw/D_GENERATION.json', rows)
        dump(R / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-crown-root-cause', phase='D_CONSTRUCTING', latest_gate=out, next_operation='Complete unchanged18-case denominator and independent scores; physical limits remain explicit', review_state='PENDING_INDEPENDENT_REVIEW'))
    p = R / ('FROZEN_PREDICTIONS_D_FIRST.json' if '--first' in sys.argv else 'FROZEN_PREDICTIONS_D.json')
    assert not p.exists()
    dump(p, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(R / 'PREREG_D.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, physical_predictions=None))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
if __name__ == '__main__':
    run()
