from common import *
from loadsets import X, C, M, pair_intervals, classes, edges, F
import zipfile, collections, resource

def arch(sel, base, R):
    meta = read(sel['metadata'])
    prov = read(V4 / 'payload/private/provenance' / f"{sel['case_key']}.json")['sources']
    arches = []
    for (jaw, src) in zip(['upper', 'lower'], prov):
        with zipfile.ZipFile(src['zip_path']) as z:
            blob = z.read(src['member'])
        if hashlib.sha256(blob).hexdigest() != src['sha256']:
            raise ValueError('ARCH_HASH_MISMATCH')
        (T, lab, coeff, ids, stats) = X.prepare(blob, meta['arches'][jaw], jaw, R, base)
        arches.append((T, lab, stats))
    return arches

def run():
    st = time.perf_counter()
    pr = read(ROOT / 'PREREG_R3.json')
    rows = []
    native = []
    checks = []
    meta = {x['key']: x for x in read(V4 / 'payload/whole_inputs/RECORDS.json')}
    for sel in pr['selection']:
        ck = sel['case_key']
        prep = npz(V4 / 'payload/whole_inputs' / (ck + '_molar_crown') / 'preparation.npz')
        (base, R) = (prep['source_base'], prep['source_R'])
        ((U, ul, us), (L, ll, ls)) = arch(sel, base, R)
        (nb, nc) = pair_intervals(U, L, ul, ll, True)
        nclass = classes(nb)
        nrob = classes(nb, 0.05)
        (uc, _) = X.affine(U)
        (lc, _) = X.affine(L)
        raster = X.make_map(U, L, ul, ll, uc, lc, h=0.1)
        active = (raster['gap'] >= 0) & (raster['gap'] <= 0.1)
        areas = {str(t): float(np.sum(active & ((raster['upper_fdi'] == t) | (raster['lower_fdi'] == t))) * 0.01) for t in M.FDI}
        total = float(active.sum() * 0.01)
        shares = {t: 100 * a / total if total else None for (t, a) in areas.items()}
        record = dict(case_key=ck, dataset=sel['dataset'], source_case=sel['case'], edges=nb, classes=nclass, robust_classes=nrob, upper_stats=us, lower_stats=ls, kernel=nc, area_control=dict(projected_per_tooth_mm2=areas, per_arch_total_mm2=total, per_tooth_share_pp=shares, band_mm=[0, 0.1], grid_mm=0.1, physical_pressure_law='UNIFORM_PRESSURE_UNCALIBRATED'))
        native.append(record)
        checks.append(dict(case=ck, kind='NATIVE', **nc['control']))
        e = edges(nb)
        gap = [(x[2] + x[3]) / 2 for x in e]
        compliance = {t: F(1, 750) for t in set((v for edge in e for v in edge[:2]))}
        qp = M.exact_point(e, gap, F(205), compliance)
        rational = M.rational_json(qp)
        good = M.replay_point(e, rational)
        wrong = dict(rational)
        wrong['total_N'] = '206'
        rejected = not M.replay_point(e, wrong)
        checks.append(dict(case=ck, kind='KKT', nominal_pass=good, injected_rejected=rejected, point_scenario='exact gap-box midpoint, T=205N, k=750N/mm; not calibrated'))
        dump(ROOT / 'raw/qp' / (ck + '.json'), rational)
        for (n, r0) in enumerate((x for x in pr['cohort'] if x['key'].startswith(ck + '_'))):
            r = dict(r0)
            t = meta[r['key']]['source_fdi']
            r.update(case_key=ck, dataset=sel['dataset'], source_fdi=t, resolution='PER_TOOTH')
            if r['status'] != 'AVAILABLE':
                r['load_status'] = 'UNKNOWN_NO_MESH'
                rows.append(r)
                continue
            m = npz(r['mesh_path'])
            tri = m['vertices'][m['faces'][m['face_roles'] == 0]]
            p = npz(V4 / 'payload/whole_inputs' / r['key'] / 'preparation.npz')
            tri = np.ascontiguousarray((tri @ p['source_R'].T + p['source_base'] - base) @ R)
            (new, kern) = pair_intervals(U, tri, ul, np.full(len(tri), t, np.int32), n == 0)
            if kern['control']:
                checks.append(dict(case=ck, uid=r['uid'], kind='CROWN', **kern['control']))
            combined = [v for v in nb if v['lower'] != t] + new
            cls = classes(combined)
            rob = classes(combined, 0.05)
            flips = [f for f in cls if cls[f]['classification'] != nclass[f]['classification']]
            rflips = [f for f in rob if rob[f]['classification'] != nrob[f]['classification']]
            r.update(load_status='CONDITIONAL_AXIAL_SETS', classes=cls, robust_classes=rob, changed_teeth=flips, robust_changed_teeth=rflips, incident_edges=new, all_edge_count=len(combined), kernel=kern, physical_force='UNKNOWN', graph_scope='ALL_RETAINED_LABELLED_NONDEGENERATE_PROJECTED_TRIANGLES; unlabelled anatomy excluded')
            rows.append(r)
            if n % 12 == 0:
                print('R3', ck, n, round(time.perf_counter() - st, 1), flush=True)
                dump(ROOT / 'raw/R3_CHECKPOINT.json', rows)
        dump(ROOT / 'raw/R3_NATIVE.json', native)
        state('R3_RUNNING', ck + ' full arch consumed', 'Continue remaining crown replacements')
    by = collections.defaultdict(list)
    for r in rows:
        if r['track'] == 'V5B_R3' and r['participant'] in FRONTIER:
            by[r['key']].append(r)
    items = []
    for (k, rr) in by.items():
        good = [r for r in rr if 'classes' in r]

        def signature(r, field):
            return tuple((r[field][str(t)]['classification'] for t in M.FDI))
        items.append(dict(key=k, requested=len(rr), scored=len(good), discriminates=len({signature(r, 'classes') for r in good}) > 1, robust_discriminates=len({signature(r, 'robust_classes') for r in good}) > 1, signatures={r['participant']: signature(r, 'classes') for r in good}))
    summary = dict(requested=len(rows), scored=sum(('classes' in r for r in rows)), crowns_changing_native_classes=sum((bool(r.get('changed_teeth')) for r in rows)), crowns_changing_robust_native_classes=sum((bool(r.get('robust_changed_teeth')) for r in rows)), class_counts=dict(collections.Counter((a['classification'] for r in rows for a in r.get('classes', {}).values()))), frontier_items=items, frontier_discriminating=sum((r['discriminates'] for r in items)), frontier_robust_discriminating=sum((r['robust_discriminates'] for r in items)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, physical_validation='UNKNOWN; no matched loaded acquisition', dropout=dict(requested=len(rows), rejected=sum(('classes' not in r for r in rows)), reasons=dict(collections.Counter((r['load_status'] for r in rows if 'classes' not in r)))))
    dump(ROOT / 'raw/R3_ROWS.json', rows)
    dump(ROOT / 'raw/R3_CONTROLS.json', checks)
    dump(ROOT / 'rounds/R3.json', summary)
    state('R3_DECIDED', str(summary['frontier_discriminating']) + ' load-class sites discriminate', 'External force comparisons, joint sufficiency probes and demo figure')
    return summary
if __name__ == '__main__':
    print(clean(run()))
