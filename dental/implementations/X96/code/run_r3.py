from common import *
from stats import collapse, proportion_ci
from full_geometry import project_cylinder
from scipy.optimize import minimize_scalar
import zipfile, gc

def segment_box(center, e, a, h, half=0.15):
    lo = center - half
    hi = center + half
    breaks = [0.0, h]
    for j in range(3):
        if a[j] != 0:
            breaks.extend((float(np.clip((v - e[j]) / a[j], 0, h)) for v in (lo[j], hi[j])))
    breaks = sorted(set(breaks))
    best = float('inf')
    bestt = None
    for (low, high) in zip(breaks[:-1], breaks[1:]):
        mid = (low + high) / 2
        p = e + mid * a
        active = (p < lo) | (p > hi)
        b = np.where(p < lo, lo, hi)
        den = float(a[active] @ a[active])
        t = float(np.clip(-a[active] @ (e - b)[active] / den, low, high)) if den else mid
        v = e + t * a
        d = float(np.linalg.norm(v - np.clip(v, lo, hi)))
        if d < best:
            best = d
            bestt = t
    if h == 0:
        v = e
        best = float(np.linalg.norm(v - np.clip(v, lo, hi)))
        bestt = 0.0
    return (best, bestt)

def axis_distance(points, e, a, h):
    dc = np.linalg.norm(points - project_cylinder(points, e, a, h, 0), axis=1)
    upper = float(dc.min())
    ix = np.where(dc - 0.15 * np.sqrt(3) <= upper + 1e-07)[0]
    vals = [segment_box(points[i], e, a, h) for i in ix]
    j = int(np.argmin([v[0] for v in vals]))
    return (vals[j][0], dict(center_zyx_mm=points[ix[j]].tolist(), t_mm=vals[j][1]), ix)

def main():
    check()
    diskcheck()
    start = time.perf_counter()
    r2 = [json.loads(l) for l in (DATA / 'R2_GEOMETRY.jsonl').read_text().splitlines()]
    assert len(r2) == 4 * len(certs())
    baseline = {(r['case'], r['fdi']): r for r in map(json.loads, (DATA / 'R1_GEOMETRY.jsonl').read_text().splitlines()) if r['extra_mm'] == 0}
    groups = {}
    for r in r2:
        groups.setdefault(r['case'], []).append(r)
    selected = sorted({r['image_group'] for r in r2}, key=lambda g: hashlib.sha256(('X96-control' + g).encode()).hexdigest())[:8]
    done_control = set()
    controls = []
    out = []
    with zipfile.ZipFile(tf2_io.ZIP) as z:
        for (i, (case, rr)) in enumerate(sorted(groups.items()), 1):
            arrays = []
            for src in rr[0]['point_sources']:
                if 'path' in src:
                    assert sha(src['path']) == src['sha256']
                    with np.load(src['path']) as a:
                        arrays.extend((a[k].astype(float) * 0.3 for k in a.files if k.endswith('_voxels_zyx')))
                else:
                    b = z.read(src['member'])
                    assert hashlib.sha256(b).hexdigest() == src['sha256']
                    (lab, sp, hdr) = tf2_io.read_mha_bytes(b)
                    arrays = [np.argwhere(np.isin(lab, [3, 4])) * 0.3]
                    del b, lab
            pt = np.unique(np.vstack(arrays), axis=0)
            for r in rr:
                e = np.array(r['entry_zyx_mm']) + r['length_mm'] * np.array(r['axis_zyx'])
                a = np.array(r['axis_zyx'])
                h = r['extra_mm']
                (d, w, active) = axis_distance(pt, e, a, h)
                b0 = baseline[r['case'], r['fdi']]
                small = min(b0['union_upper_mm'], d)
                assert r['union_lower_mm'] <= small + 1e-06, 'Tool containment failed'
                out.append(dict(case=case, fdi=r['fdi'], image_group=r['image_group'], paired=r['paired_revisions'], extra_mm=h, accepted=r['nominal_tf2_lower_mm'] >= 2, axis_extension_gap_mm=d, minimal_tool_gap_mm=small, full_envelope_lower_mm=r['union_lower_mm'], full_envelope_upper_mm=r['union_upper_mm'], all_shapes_breach=small < 2, any_shape_breach=r['union_upper_mm'] < 2, shape_ambiguous=small >= 2 and r['union_upper_mm'] < 2, axis_witness=w, resolution='PER_TOOTH', time_scale='SIMULTANEOUS'))
                cg = (r['image_group'], h)
                if r['image_group'] in selected and cg not in done_control:
                    vals = []
                    for j in active:
                        lo = pt[j] - 0.15
                        hi = pt[j] + 0.15
                        fn = lambda t: float(np.sum((e + t * a - np.clip(e + t * a, lo, hi)) ** 2))
                        opt = minimize_scalar(fn, bounds=(0, h), method='bounded', options={'xatol': 1e-12})
                        vals.append(np.sqrt(min(fn(0), fn(h), opt.fun)))
                    c = float(min(vals))
                    controls.append(dict(image_group=cg[0], extra_mm=h, error_mm=abs(c - d), injected_plus1mm_rejected=abs(c - (d + 1)) > 1e-06))
                    done_control.add(cg)
            if i % 40 == 0:
                dump(ROOT / 'raw/R3_PARTIAL.json', dict(cases=i, rows=len(out)))
                state('R3_RUNNING', 'IN_PROGRESS', 'Complete axis geometry', cases=i, rows=len(out))
                print('R3', i, len(out), round(time.perf_counter() - start, 1), flush=True)
            del pt, arrays
            gc.collect()
    dest = DATA / 'R3_TOOL_SHAPE.jsonl'
    dest.write_text(''.join((json.dumps(r) + '\n' for r in out)))
    dump(ROOT / 'raw/R3_CONTROLS.json', controls)
    tables = []
    for h in [0.5, 1.0, 1.2, 1.3]:
        rows = [r for r in out if r['extra_mm'] == h]
        unique = {}
        for r in rows:
            unique.setdefault((r['image_group'], r['fdi']), []).append(r)
        cc = []
        for ((g, f), rr) in unique.items():
            acc = [r for r in rr if r['accepted']]
            cc.append(dict(image_group=g, fdi=f, accepted=bool(acc), lower=any((r['all_shapes_breach'] for r in acc)), upper=any((r['any_shape_breach'] for r in acc)), ambiguous=any((r['any_shape_breach'] for r in acc)) and (not any((r['all_shapes_breach'] for r in acc))), paired=any((r['paired'] for r in rr))))
        tables.append(dict(extra_mm=h, all=dict(lower=proportion_ci(cc, 'lower'), upper=proportion_ci(cc, 'upper'), ambiguous=proportion_ci(cc, 'ambiguous')), paired=dict(lower=proportion_ci([r for r in cc if r['paired']], 'lower'), upper=proportion_ci([r for r in cc if r['paired']], 'upper'), ambiguous=proportion_ci([r for r in cc if r['paired']], 'ambiguous')), scope='identified shape interval for explicit family and maximum-depth scenario; not true physical rate'))
    p = np.array([[12.0, 3.5, 0]])
    e = np.zeros(3)
    a = np.array([1.0, 0, 0])
    nom = float(np.linalg.norm(p - project_cylinder(p, e, a, 10, 2)))
    full = float(np.linalg.norm(p - project_cylinder(p, e, a, 12, 2)))
    axis = float(np.linalg.norm(p - project_cylinder(p, np.array([10.0, 0, 0]), a, 2, 0)))
    suff = dict(summary={'nominal_gap_mm': [nom, nom], 'extra_depth_mm': [2.0, 2.0]}, identity_error_mm=0.0, downstream_gap_mm=[min(nom, axis), full], downstream_difference_mm=abs(min(nom, axis) - full), decision_changed=min(nom, axis) >= 2 and full < 2, smallest_extension='radial occupancy profile along extra segment or frozen-query signed tool clearance; axis/max-depth alone insufficient', external_referent={'kind': 'our_own_fixture', 'locator': 'code/run_r3.py sufficiency witness', 'compared_quantity': 'identical-summary tool-shape counterexample', 'refutes_us': True})
    dump(ROOT / 'raw/SUFFICIENCY_R3.json', suff)
    res = dict(round='R3', claim_type='information_link', outcome='EXTRA_DEPTH_ALONE_DOES_NOT_IDENTIFY_BREACH_FREQUENCY_TOOL_SHAPE_INTERVAL_EXECUTED', tables=tables, gates=dict(coverage=len(out) == 4 * len(certs()), controls=all((r['error_mm'] <= 1e-06 and r['injected_plus1mm_rejected'] for r in controls)), containment=True, decisive_shape_ambiguity=all((t['all']['ambiguous']['fraction'] >= 0.05 for t in tables)), summary_sufficiency=not suff['decision_changed']), external_referent=load(ROOT / 'PREREG_R3.json')['external_referent'], cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, threads=1, gpu=False, scratch_bytes=diskcheck()), artifact=dict(path=str(dest), sha256=sha(dest), bytes=dest.stat().st_size), formal_floating_point_enclosure='MISSING')
    dump(ROOT / 'rounds/R3.json', res)
    state('R3_COMPLETE', res['gates'], 'Prepare prospective tool-profile and signed-clearance measurement port')
    (ROOT / 'HANDOFF_R3.md').write_text('# R3 borrformens informationsintervall\n\n' + json.dumps(res, indent=2) + "\n\nAxialsvep and full radius housing limit an explicit drilling family. Real profile and registered drill pose must be measured; extreme shapes are logical control objects. Next design: independent tool silhouette in the same length date and signed channel distance.\n")
    print(json.dumps(dict(tables=[dict(h=t['extra_mm'], lower=t['all']['lower']['fraction'], upper=t['all']['upper']['fraction'], ambiguous=t['all']['ambiguous']['fraction']) for t in tables], gates=res['gates']), indent=2))
if __name__ == '__main__':
    main()
