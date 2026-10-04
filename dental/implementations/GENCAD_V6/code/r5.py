from common import *
from contact import Q, compare
import csv, gzip, collections, itertools

def linear_integral(poly, xy, values):
    if len(poly) < 3:
        return 0.0
    q = poly @ xy
    v = poly @ values
    areas = abs(np.cross(q[1:-1] - q[0], q[2:] - q[0])) / 2
    return float(np.sum(areas * (v[0] + v[1:-1] + v[2:]) / 3))

def trajectory(xy, faces, pg, rg, b=0.1):
    valid = np.isfinite(pg[faces]).all(1) & np.isfinite(rg[faces]).all(1)
    ff = faces[valid]
    value = 0.0
    support = 0.0
    for ids in ff:
        e = pg[ids] - rg[ids]
        q = xy[ids]
        a = Q.area(np.eye(3), q)
        support += a
        if np.all(e >= b) or np.all(e <= -b):
            value += 2 * b * a
            continue
        if np.all(e >= 0) and np.all(e <= b):
            value += 2 * a * e.mean()
            continue
        if np.all(e <= 0) and np.all(e >= -b):
            value -= 2 * a * e.mean()
            continue
        pos = Q.clip(Q.clip(np.eye(3), e), b - e)
        neg = Q.clip(Q.clip(np.eye(3), -e), b + e)
        z = linear_integral(pos, q, e) - linear_integral(neg, q, e)
        if e.max() > b:
            z += b * Q.area(Q.clip(np.eye(3), e - b), q)
        if e.min() < -b:
            z += b * Q.area(Q.clip(np.eye(3), -e - b), q)
        value += 2 * z
    return dict(J_mm3=value if support else None, support_mm2=support, normalized_disagreement=value / (2 * b * support) if support else None, interpretation='Analytic all-offset contact-filtration difference, not physical swept volume')

def controls():
    xy = np.array([[0.0, 0], [1, 0], [0, 1]])
    f = np.array([[0, 1, 2]])
    pg = np.array([-0.12, 0.07, 0.18])
    rg = np.array([-0.02, -0.04, 0.13])
    a = trajectory(xy, f, pg, rg)['J_mm3']
    ss = np.linspace(-0.4, 0.4, 4001)
    vv = [Q.contact_map(xy, f, pg + s, rg + s)['contact_symdiff_mm2'] for s in ss]
    numeric = float(np.trapz(vv, ss))
    shift = trajectory(xy, f, pg + 4, rg + 4)['J_mm3']
    selfv = trajectory(xy, f, rg, rg)['J_mm3']
    xx = np.concatenate([xy, xy + [2, 0]])
    ff = np.array([[0, 1, 2], [3, 4, 5]])
    ga = np.array([0.0, 0, 0, 0.2, 0.2, 0.2])
    gb = np.array([0.2, 0.2, 0.2, 0, 0, 0])
    rr = np.zeros(6)
    ja = trajectory(xx, ff, ga, rr)['J_mm3']
    jb = trajectory(xx, ff, gb, rr)['J_mm3']
    down = Q.contact_map(xx, ff, ga, gb)['contact_symdiff_mm2']
    suff = dict(summary='integrated all-offset contact disagreement J', A_J_mm3=ja, B_J_mm3=jb, identity_error_mm3=abs(ja - jb), downstream_current_mask_difference_mm2=down, minimum_extension='location-labelled activation field; scalar J insufficient for local load/repair query', resolution='PER_SURFACE_REGION', external_referent=dict(kind='our_own_fixture', locator='code/r5.py::controls', compared_quantity='J-summary sufficiency, not physical validation', refutes_us=True))
    dump(ROOT / 'raw/SUFFICIENCY_TRAJECTORY.json', suff)
    out = dict(analytic_J_mm3=a, numerical_J_mm3=numeric, quadrature_error_mm3=abs(a - numeric), shared_shift_error_mm3=abs(a - shift), self_identity_error_mm3=selfv, wrong_integral_rejected=abs(numeric - (a + 1)) > 1e-08, summary_identity_error_mm3=abs(ja - jb), downstream_mask_difference_mm2=down, all_pass=abs(a - numeric) < 1e-08 and abs(a - shift) < 1e-12 and (selfv == 0.0) and (ja == jb) and (down == 1.0))
    if not out['all_pass']:
        raise AssertionError(out)
    return out

def run():
    st = time.perf_counter()
    p = read(ROOT / 'PREREG_R5.json')['metric']
    ctrl = controls()
    roof = []
    raw = read(ROOT / 'raw/R1_ROWS.json')
    cases = sorted({r['case_key'] for r in raw})
    old = {(r['task_id'], r['participant']): r for r in raw}
    for key in cases:
        tasks = read(V4 / 'payload/public/tasks' / f'{key}.json')
        ref = npz(V4 / 'payload/private/references' / f'{key}.npz')
        pred = {m: npz(V4 / 'payload/predictions' / m / f'{key}.npz') for m in FRONTIER}
        cache = {}
        for (i, t) in enumerate(tasks):
            for name in FRONTIER:
                o = old[t['task_id'], name]
                row = dict(task_id=t['task_id'], case_key=key, participant=name, status='UNKNOWN', L1=o['L1'], native_band_nonempty=o.get('contact', {}).get('reference', {}).get('area_mm2', 0) > 1e-12)
                if f'outer_{i}' in pred[name]:
                    if t['family'] not in cache:
                        cache[t['family']] = npz(V4 / 'payload/public' / t['geometry_file'])
                    z = cache[t['family']]
                    a = trajectory(z['xy'], z['faces'], z['ceiling'] - pred[name][f'outer_{i}'], z['ceiling'] - ref[t['family']])
                    row.update(a)
                    if a['J_mm3'] is not None:
                        row['status'] = 'PASS' if a['J_mm3'] <= p['integrated_contact_difference_max_mm3'] + p['numeric_tolerance_mm3'] and o['contact']['coverage'] >= 0.95 else 'FAIL'
                roof.append(row)
    whole = []
    for o in read(ROOT / 'raw/R2_ROWS.json'):
        r = dict(uid=o['uid'], key=o['key'], track=o['track'], participant=o['participant'], status='UNKNOWN', native_band_nonempty=o.get('contact', {}).get('reference', {}).get('area_mm2', 0) > 1e-12)
        if o.get('field_path'):
            z = npz(o['field_path'])
            a = trajectory(z['xy'], z['faces'], z['pred_gap'], z['reference_gap'])
            r.update(a)
            if a['J_mm3'] is not None:
                r['status'] = 'PASS' if a['J_mm3'] <= p['integrated_contact_difference_max_mm3'] + p['numeric_tolerance_mm3'] and o['contact']['coverage'] >= 0.95 else 'FAIL'
        whole.append(r)

    def summarize(rows, key):
        by = collections.defaultdict(list)
        for r in rows:
            by[r[key]].append(r)
        items = [dict(key=k, discriminates=len({r['status'] == 'PASS' for r in rr}) > 1, native_band_nonempty=any((r['native_band_nonempty'] for r in rr)), statuses={r['participant']: r['status'] for r in rr}) for (k, rr) in by.items()]
        return dict(items=items, requested_items=len(items), discriminating=sum((r['discriminates'] for r in items)), discriminating_nonempty_reference=sum((r['discriminates'] and r['native_band_nonempty'] for r in items)), pass_count=sum((r['status'] == 'PASS' for r in rows)), dropout=sum((r['status'] == 'UNKNOWN' for r in rows)), scored_rows=len(rows) - sum((r['status'] == 'UNKNOWN' for r in rows)))
    out = dict(roof=summarize(roof, 'task_id'), whole=summarize([r for r in whole if r['participant'] in FRONTIER and r['track'] == 'V5B_R3'], 'key'), whole_all_methods=dict(collections.Counter((r['status'] for r in whole))), controls=ctrl, seconds=time.perf_counter() - st)
    dump(ROOT / 'raw/R5_ROOF_ROWS.json', roof)
    dump(ROOT / 'raw/R5_WHOLE_ROWS.json', whole)
    dump(ROOT / 'rounds/R5.json', out)
    state('R5_DECIDED', 'Contact filtration integrated analytically across all offsets', 'Report which questions remain unsupported by physical measurements')
    (ROOT / 'history/HANDOFF_R5.md').write_text('R5 changes the observation operation to the entire ideal axial contact filtration. The equivalence to clipped absolute gap error is explicitly derived; this is not a new optimization algorithm or physiological motion. Earlier snapshot failures remain unchanged.\n')
    return {k: {x: v for (x, v) in r.items() if x != 'items'} if isinstance(r, dict) else r for (k, r) in out.items()}
if __name__ == '__main__':
    print(clean(run()))
