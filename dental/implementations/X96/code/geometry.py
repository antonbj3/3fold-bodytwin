from common import *
import zipfile, gc, traceback

def run(depths, round_name):
    check()
    diskcheck()
    start = time.perf_counter()
    dst = DATA / (round_name + '_GEOMETRY.jsonl')
    old = {}
    if dst.exists():
        cp = ROOT / 'raw' / f'{round_name}_COST.json'
        if cp.exists():
            assert sha(dst) == load(cp)['artifact']['sha256'], 'Completed geometry artifact changed'
        for l in dst.read_text().splitlines():
            r = json.loads(l)
            old[r['case'], r['fdi'], r['extra_mm']] = r
    ss = certs()
    pairs = load(X73 / 'raw/LINEAGE_SITES.json')
    pp = {(r['case'], r['fdi']): r for r in pairs}
    groups = {}
    for s in ss:
        groups.setdefault(s['case'], []).append(s)
    fails = []
    processed = 0
    with zipfile.ZipFile(tf2_io.ZIP) as z, dst.open('a') as out:
        for (case, sites) in sorted(groups.items()):
            todo = [s for s in sites if any(((case, s['fdi'], h) not in old for h in depths))]
            if not todo:
                continue
            try:
                ps = pp.get((case, todo[0]['fdi']))
                casepair = next((pp[case, s['fdi']] for s in sites if (case, s['fdi']) in pp), None)
                provenance = []
                if casepair:
                    arts = [casepair['point_artifact']] + [r['point_artifact'] for r in casepair['alias_source_bindings']]
                    pointlists = []
                    tf2 = None
                    for art in arts:
                        assert sha(art['path']) == art['sha256']
                        provenance.append(art)
                        with np.load(art['path']) as a:
                            for k in a.files:
                                if k.endswith('_voxels_zyx'):
                                    pt = a[k].astype(float) * 0.3
                                    pointlists.append(pt)
                                    if k == 'tf2_voxels_zyx' and tf2 is None:
                                        tf2 = pt
                    union = np.unique(np.vstack(pointlists), axis=0)
                    image_group = casepair['image_group']
                else:
                    b = z.read(todo[0]['source_member'])
                    assert hashlib.sha256(b).hexdigest() == todo[0]['source_member_sha256']
                    (lab, sp, hdr) = tf2_io.read_mha_bytes(b)
                    assert tuple(sp) == (0.3, 0.3, 0.3) and hdr.get('Offset') == '0 0 0' and (hdr.get('TransformMatrix') == '1 0 0 0 1 0 0 0 1')
                    tf2 = np.argwhere(np.isin(lab, [3, 4])).astype(float) * np.array(sp)
                    union = tf2
                    image_group = case
                    provenance = [dict(archive=tf2_io.ZIP, member=todo[0]['source_member'], sha256=todo[0]['source_member_sha256'])]
                    del b, lab
                for s in todo:
                    e = np.array(s['entry_zyx_mm'])
                    a = np.array(s['axis_zyx'])
                    L = s['length_mm']
                    r = s['radius_mm']
                    assert abs(np.linalg.norm(a) - 1) < 1e-12
                    baseline = pp.get((case, s['fdi']))
                    assert sha(DENT / 'results/NV1_canals/per_case' / f'{case}.json') == s['plan_file_sha256']
                    for h in depths:
                        key = (case, s['fdi'], h)
                        if key in old:
                            continue
                        br = voxel_cylinder_bracket(union, np.full(3, 0.3), e, a, L + h, r)
                        assert br['gap_mm'] <= 0.0001
                        ref = voxel_cylinder_bracket(tf2, np.full(3, 0.3), e, a, L + h, r) if casepair else br
                        if h == 0:
                            assert max((abs(ref[q] - s[q]) for q in ['lower_mm', 'upper_mm'])) <= 0.0001
                        baseunion = baseline['union']['lower_mm'] if baseline else s['lower_mm']
                        assert br['lower_mm'] >= max(0, baseunion - h) - 1e-07 and br['upper_mm'] <= baseunion + 0.0001
                        rr = dict(case=case, fdi=s['fdi'], image_group=image_group, paired_revisions=bool(casepair), extra_mm=h, entry_zyx_mm=e.tolist(), axis_zyx=a.tolist(), length_mm=L, radius_mm=r, nominal_tf2_lower_mm=s['lower_mm'], nominal_tf2_upper_mm=s['upper_mm'], nominal_union_lower_mm=baseunion, tf2_lower_mm=ref['lower_mm'], tf2_upper_mm=ref['upper_mm'], union_lower_mm=br['lower_mm'], union_upper_mm=br['upper_mm'], numeric_gap_mm=br['gap_mm'], witness=br['witness'], point_sources=provenance, tool='FULL_RADIUS_ENVELOPE', resolution='PER_TOOTH', time_scale='SIMULTANEOUS')
                        out.write(json.dumps(rr) + '\n')
                        out.flush()
                        old[key] = rr
                    processed += 1
                del tf2, union
                if casepair:
                    del pointlists
                gc.collect()
                if processed % 50 < len(todo):
                    print(round_name, case, 'new sites', processed, 'records', len(old), 'seconds', round(time.perf_counter() - start, 1), flush=True)
                    state(round_name + '_RUNNING', 'IN_PROGRESS', 'Complete source-bound geometry', records=len(old), new_sites=processed)
            except Exception as exc:
                fails.append(dict(case=case, error=repr(exc), traceback=traceback.format_exc()))
                dump(ROOT / 'raw' / f'{round_name}_FAILURES.json', fails)
                raise
    costfile = ROOT / 'raw' / f'{round_name}_COST.json'
    previous = load(costfile) if costfile.exists() else {}
    elapsed = time.perf_counter() - start
    history = ROOT / 'raw' / f'{round_name}_COST_HISTORY.json'
    calls = load(history) if history.exists() else [previous] if previous else []
    cost = dict(wall_seconds=elapsed, cold_geometry_wall_seconds=previous.get('cold_geometry_wall_seconds', previous.get('wall_seconds', elapsed)), cumulative_geometry_wall_seconds=previous.get('cumulative_geometry_wall_seconds', previous.get('wall_seconds', 0)) + elapsed, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, threads=1, gpu=False, records=len(old), new_sites=processed, scratch_bytes=diskcheck(), failures=fails, artifact=dict(path=str(dst), sha256=sha(dst), bytes=dst.stat().st_size))
    calls.append(cost)
    dump(history, calls)
    dump(costfile, cost)
    state(round_name + '_GEOMETRY_COMPLETE', 'GEOMETRY_COVERAGE_COMPLETE', 'Summarize and execute independent controls', records=len(old))
    return list(old.values())
if __name__ == '__main__':
    run(load(ROOT / 'PREREG_R1.json')['selection']['depths_mm'], 'R1')
