from geometry import *

def collar(m, width=0.8):
    ll = loops(m)
    if len(ll) != 1:
        raise ValueError('native mesh has%d boundary loops' % len(ll))
    cv = m.vertices[ll[0]]
    d = cKDTree(cv).query(m.vertices, workers=1)[0]
    sel = np.all(d[m.faces] <= width, axis=1)
    if sel.sum() < 30:
        raise ValueError('too few collar faces')
    return (cv, m.triangles[sel])

def run():
    st = time.perf_counter()
    lock = read(ROOT / 'FROZEN_TEST_SET.json')
    rows = []
    sources = lock['sources']
    for s in sources:
        assert sha(s['upper_path']) == s['upper_sha256']
        (v, f) = read_obj(s['upper_path'])
        v = v @ np.array(s['frame']).T
        lab = np.array(read(s['label_path'])['labels'])
        (lv, lf) = read_obj(s['lower_path'])
        lv = lv @ np.array(s['frame']).T
        for (fdi, family) in [(11, 'anterior'), (14, 'premolar'), (16, 'molar')]:
            rec = dict(key=s['case_key'] + '_' + family, case_key=s['case_key'], fdi=fdi, family=family, status='PREPARED')
            try:
                (t, dt, nt) = component(v[f[np.all(lab[f] == fdi, axis=1)]])
                (d, dd, nd) = component(v[f[np.all(lab[f] == fdi + 10, axis=1)]])
                if max(dt, dd) > 0.05:
                    raise ValueError('component loss>5%')
                (curve, tc) = collar(t)
                (dc, donorc) = collar(d)
                origin = curve.mean(0)
                center = curve[:, :2].mean(0)
                apex = np.r_[center, t.vertices[:, 2].max() - 1.5]
                sides = [v[f[np.all(lab[f] == i, axis=1)]] for i in ([21, 12] if fdi == 11 else [fdi - 1, fdi + 1])]
                lo = t.bounds[0] - 3
                hi = t.bounds[1] + 3
                lt = lv[lf]
                keep = np.all(lt[:, :, :2].max(1) >= lo[:2], axis=1) & np.all(lt[:, :, :2].min(1) <= hi[:2], axis=1)
                p = dict(donor=d.triangles - origin, donor_collar=donorc - origin, collar=tc - origin, margin_curve=curve - origin, cavity_center=center - origin[:2], cavity_apex=apex - origin, cavity_scale=np.array(0.65), mesial=sides[0] - origin, distal=sides[1] - origin, antagonist=lt[keep] - origin, origin_world=origin)
                for (folder, arr) in [('public_b', p), ('private_b', dict(target=t.triangles - origin))]:
                    (DATA / folder).mkdir(exist_ok=True)
                    np.savez_compressed(DATA / folder / (rec['key'] + '.npz'), **arr)
                rec.update(public_path=DATA / 'public_b' / (rec['key'] + '.npz'), private_path=DATA / 'private_b' / (rec['key'] + '.npz'), target_drop=dt, donor_drop=dd, target_faces=len(t.faces), public_collar_faces=len(tc), collar_area_fraction=float(mesh(tc).area / t.area), native_margin_z_span_mm=float(np.ptp(curve[:, 2])), antagonist_faces=int(keep.sum()))
                rec.update(public_sha256=sha(rec['public_path']), private_sha256=sha(rec['private_path']))
            except Exception as e:
                rec.update(status='PREPARATION_REJECTED', reason=repr(e))
            rows.append(rec)
            print(rec['key'], rec['status'], rec.get('reason', ''), flush=True)
    bank = []
    for r in read(BASE / 'PROOF_LANE_FULL_CROWN_DIAG/raw/ANNOTATED_PAIRS.json')['rows']:
        family = r['family'].replace('_crown', '')
        src = DATA.parent / 'PROOF_LANE_FULL_CROWN_DIAG' / (r['key'] + '_annotated.npz')
        tri = npz(src)['target']
        (m, drop, nc) = component(tri)
        try:
            (cv, tc) = collar(m)
        except Exception as e:
            bank.append(dict(key=r['key'], status='REJECTED', reason=repr(e)))
            continue
        out = DATA / 'templates' / (r['key'] + '.npz')
        out.parent.mkdir(exist_ok=True)
        np.savez_compressed(out, triangles=m.triangles, collar=tc)
        bank.append(dict(key=r['key'], family=family, status='AVAILABLE', path=out, sha256=sha(out), source_path=src, source_sha256=sha(src)))
    freeze(ROOT / 'FROZEN_INPUTS_B.json', dict(records=rows, templates=bank, sources=sources, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    state('B_INPUTS_READY', dict(prepared=sum((r['status'] == 'PREPARED' for r in rows)), requested=18), 'Generate B and rigid-control surfaces; freeze predictions')
if __name__ == '__main__':
    run()
