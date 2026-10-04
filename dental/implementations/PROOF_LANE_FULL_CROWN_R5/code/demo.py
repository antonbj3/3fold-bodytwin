"""Rebuild one frozen-case representative per family, seal, then score and verify."""
from construct_b import *
from score_a import exact_failure
import argparse

def run():
    ap = argparse.ArgumentParser()
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--method', choices=['A', 'B', 'C'], default='A')
    args = ap.parse_args()
    start = time.perf_counter()
    for row in read(R / 'DEPENDENCIES.json')['files']:
        assert sha(row['path']) == row['sha256'], row['path']
    dest = D / 'replays' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    dest.mkdir(parents=True)
    chosen = inputs() if args.full else inputs()[:3]
    pred = []
    for (rec, rr) in chosen:
        m = npz(rec['mesh_path'])
        t = time.perf_counter()
        r = dict(key=rec['key'], family=rec['family'], status='REJECTED')
        try:
            if args.method == 'A':
                (vv, off) = normal_offset(m['vertices'], m['faces'], m['roles'])
                (v, repair) = project(vv, m['faces'], m['roles'])
                (f, roles) = (m['faces'], m['roles'])
            elif args.method == 'B':
                (iv, inf, off, _, _) = make_inner(m)
                (ov, of, repair) = outer_repair(m, iv[inf])
                (v, f, roles) = join(ov, of, iv, inf)
            else:
                from construct_c import implicit_inner
                (iv, inf, off) = implicit_inner(m, rec['key'] + '_replay')
                (ov, of, repair) = outer_repair(m, iv[inf])
                (v, f, roles) = join(ov, of, iv, inf)
            path = dest / (rec['key'] + '.npz')
            np.savez_compressed(path, vertices=v, faces=f, roles=roles)
            r.update(status='GENERATED', path=path, sha256=sha(path), offset=off, repair=repair)
        except Exception as exc:
            r['reason'] = repr(exc)
        r['seconds'] = time.perf_counter() - t
        pred.append(r)
    freeze(dest / 'FROZEN_PREDICTIONS.json', dict(method=args.method, rows=pred, prereg_sha256=sha(R / f'PREREG_{args.method}.json'), code_sha256=sha(Path(__file__))))
    import geometry as geo
    old = geo.closest
    geo.closest = fast_nearest
    scorer.closest = fast_nearest
    rows = []
    thresholds = read(R / 'PREREG_A.json')['metrics']['shape_threshold_mm']
    for (rec, (source, rr)) in zip(pred, chosen):
        out = dict(rec)
        if rec['status'] == 'GENERATED':
            m = npz(rec['path'])
            ext = m['vertices'][m['faces'][m['roles'] == 0]]
            inner = m['vertices'][m['faces'][m['roles'] == 1]]
            target = npz(rr['private_path'])['target']
            p = npz(rr['public_path'])
            form = scorer.metrics(ext, target)
            wall = scorer.wall(ext, inner)
            (q, d, j) = fast_nearest(inner, ext.mean(1))
            i = int(d.argmin())
            wit = exact_failure(ext[i], inner[j[i]])
            om = trimesh.Trimesh(m['vertices'], m['faces'][m['roles'] == 0], process=False)
            margin = curve_error(om.vertices[loops(om)[0]], p['margin_curve'])
            mm = trimesh.Trimesh(m['vertices'], m['faces'], process=False)
            out.update(shape=form, wall=wall, exact_wall_witness=wit, margin_mm=margin, complete=bool(form['p95_mm'] <= thresholds[rec['family']] and wall['continuous_lower_mm'] >= 0.5 and (margin <= 0.025) and mm.is_watertight and (rec['offset']['max_plane_offset_residual_mm'] <= 0.01)))
        rows.append(out)
    geo.closest = old
    from verify_exclusion import run as verify
    verification = verify()
    from controls import run as checks
    control = checks()
    result = dict(rows=rows, complete_count=sum((r.get('complete', False) for r in rows)), exact_exclusion_verified=verification['verified_rational_witnesses'], controls_pass=control['all_pass'], seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, source='same frozen R4 six-jaw cohort, first jaw by frozen order for default3-family demo', limitations='virtual preparation, no clinical function or physical cement validation', prediction_manifest_sha256=sha(dest / 'FROZEN_PREDICTIONS.json'))
    save(dest / 'REPLAY_RESULTS.json', result)
    save(R / 'raw/LAST_DEMO.json', dict(path=dest / 'REPLAY_RESULTS.json', sha256=sha(dest / 'REPLAY_RESULTS.json'), seconds=result['seconds'], peak_rss_MiB=result['peak_rss_MiB'], complete_count=result['complete_count'], controls_pass=result['controls_pass']))
    print(json.dumps(clean(dict(replay=dest, complete=result['complete_count'], controls=result['controls_pass'], exact_witnesses=result['exact_exclusion_verified'], seconds=result['seconds'], peak_rss_MiB=result['peak_rss_MiB'])), indent=2))
if __name__ == '__main__':
    run()
