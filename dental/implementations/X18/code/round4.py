"""Measured same-patient morphology and measured antagonist: factorial input test."""
import sys, time, datetime, resource
from geometry import *
from experiment import contact, score, local_tri
from obstacle import constraints, lp_control

def reflected_prior(arch, fdi, s, n):
    opposite = fdi + 10 if fdi // 10 == 3 else fdi - 10
    donor = site(arch, opposite, n)
    xy = donor['center'] + s['uv'] * np.array([-1, 1]) * donor['half']
    tri = arch['tri'][arch['owner'] == opposite]
    (height, faces) = query_height(tri, xy, False)
    ok = np.isfinite(height)
    if ok.mean() < 0.85:
        raise ValueError('Contralateral measured roof coverage<85%')
    if not ok.all():
        (distance, ix) = cKDTree(xy[ok]).query(xy[~ok])
        if distance.max() > 0.8:
            raise ValueError('Contralateral prior hole requires>0.8mm interpolation')
        height[~ok] = height[ok][ix]
    apex = np.quantile(arch['v'][arch['labels'] == opposite, 2], 0.95)
    return (s['z_top'] + height - apex, dict(fdi=opposite, measured_vertex_fraction=float(ok.mean()), interpolation_status='Nearest observed donor roof for small mask holes; constitutive closure', donor_site_center_mm=donor['center'], donor_source_faces=faces))

def informed(pr, U, s, A, b, prior, ceiling):
    first = None
    nonempty = None
    attempts = []
    for lift in pr['lift_candidates_mm']:
        t = time.perf_counter()
        (z, info) = lp_control(A, b, prior + lift)
        if z is None:
            attempts.append(dict(lift_mm=lift, LP=info, seconds=time.perf_counter() - t))
            continue
        (z, rep) = repair(U, s, z, 0.025)
        co = contact(ceiling - z, z, s['xy'], s['index'], s['area_weight'], 0.1, 0.1)
        sel = pr['selection']
        eligible = sel['patch_min'] <= co['patch_count'] <= sel['patch_max'] and sel['area_min_mm2'] <= co['area_mm2'] <= sel['area_max_mm2']
        attempts.append(dict(lift_mm=lift, LP=info, area_mm2=co['area_mm2'], patch_count=co['patch_count'], selection_gate=eligible, repair=rep, seconds=time.perf_counter() - t))
        if first is None:
            first = (z, lift)
        if nonempty is None and co['patch_count']:
            nonempty = (z, lift)
        if eligible:
            return (z, lift, attempts)
    ans = nonempty or first
    return (ans[0], ans[1], attempts) if ans else (None, None, attempts)

def predict():
    pr = json.loads((H / 'PREREG_R4.json').read_text())
    start = time.perf_counter()
    generic = np.load(D / 'donor_template.npz')
    rows = []
    for case in pr['test_cases']:
        state('R4_DESIGNING', 'R2 contact recovery FAIL; R3 force uncertainty quantified', 'Measured homologous morphology ' + str(case), completed=len(rows))
        (data, man) = pair(case)
        arch = data['lower']
        for fdi in pr['tooth_fdi']:
            t = time.perf_counter()
            try:
                s = site(arch, fdi, 25)
                (personal, donor) = reflected_prior(arch, fdi, s, 25)
                prior = s['z_top'] + generic[str(fdi % 10)]
                U = local_tri(data['upper']['tri'], s['xy'])
                (ceiling, uf) = query_height(U, s['xy'], True)
                (A, b) = constraints(U, s['xy'], s['faces'], 0.025)
                (g, gl, ga) = informed(pr, U, s, A, b, prior, ceiling)
                (p, pl, pa) = informed(pr, U, s, A, b, personal, ceiling)
                if g is None or p is None:
                    raise ValueError('All LP lift candidates failed for one prior')
                path = D / (str(case) + '_' + str(fdi) + '_R4_pred.npz')
                np.savez_compressed(path, xy=s['xy'], uv=s['uv'], faces=s['faces'], border=s['border'], index=s['index'], weights=s['area_weight'], ceiling=ceiling, z_generic_without=prior, z_generic_with=g, z_personal_without=personal, z_personal_with=p)
                for (arm, z) in [('generic_without', prior), ('generic_with', g), ('personal_without', personal), ('personal_with', p)]:
                    (v, f) = shell(s['xy'], z, s['faces'], 1.2)
                    write_stl(D / (str(case) + '_' + str(fdi) + '_R4_' + arm + '.stl'), v, f)
                rows.append(dict(case=case, fdi=fdi, status='DESIGNED', file=str(path), sha256=sha(path), source=man, donor=donor, lift_generic_mm=gl, lift_personal_mm=pl, generic_attempts=ga, personal_attempts=pa, constraints=A.shape[0], design_seconds=time.perf_counter() - t))
            except (ValueError, IndexError) as e:
                rows.append(dict(case=case, fdi=fdi, status='UNKNOWN_SITE_OR_PRIOR', reason=str(e)))
            dump(H / 'raw/PREDICTIONS_R4_checkpoint.json', rows)
        print('R4 designed', case, len(rows), flush=True)
    dump(H / 'raw/PREDICTIONS_R4.json', rows)
    dump(H / 'FROZEN_PREDICTIONS_R4.json', dict(round='R4', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R4.json'), predictions_sha256=sha(H / 'raw/PREDICTIONS_R4.json'), files=[dict(path=r['file'], sha256=r['sha256']) for r in rows if r['status'] == 'DESIGNED'], code={p.name: sha(p) for p in (H / 'code').glob('*.py')}, target_contacts_read=False))
    dump(H / 'raw/PREDICT_R4_COST.json', dict(wall_seconds=time.perf_counter() - start, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    state('R4_FROZEN', 'New original contacts still unread', 'Evaluate factorial contact errors')

def evaluate():
    pr = json.loads((H / 'PREREG_R4.json').read_text())
    fr = json.loads((H / 'FROZEN_PREDICTIONS_R4.json').read_text())
    assert sha(H / 'PREREG_R4.json') == fr['prereg_sha256']
    assert sha(H / 'raw/PREDICTIONS_R4.json') == fr['predictions_sha256']
    preds = json.loads((H / 'raw/PREDICTIONS_R4.json').read_text())
    rows = []
    current = None
    start = time.perf_counter()
    for r in preds:
        if r['status'] != 'DESIGNED':
            rows.append(r)
            continue
        if current != r['case']:
            (data, _) = pair(r['case'])
            current = r['case']
            state('R4_EVALUATING', 'Four-arm predictions frozen', 'Original contacts ' + str(current))
        assert sha(r['file']) == r['sha256']
        z = np.load(r['file'])
        xy = z['xy']
        C = data['lower']['tri'][data['lower']['owner'] == r['fdi']]
        (original, faces) = query_height(C, xy, False)
        U = local_tri(data['upper']['tri'], xy)
        raw = dict(xy=xy, original_z=original, original_faces=faces, original_gap=z['ceiling'] - original)
        row = dict(case=r['case'], fdi=r['fdi'], status='SCORED', type={4: 'P1', 5: 'P2', 6: 'M1', 7: 'M2'}[r['fdi'] % 10], arms={}, continuous={}, design_seconds=r['design_seconds'])
        for band in [0.1, 0.05, 0.2]:
            key = str(band)
            truth = contact(z['ceiling'] - original, original, xy, z['index'], z['weights'], band, 0.1)
            row['arms'][key] = {'original': truth}
            for arm in pr['arms']:
                zz = z['z_' + arm]
                co = contact(z['ceiling'] - zz, zz, xy, z['index'], z['weights'], band, 0.1)
                co['errors'] = score(co, truth, xy)
                row['arms'][key][arm] = co
            for (arm, co) in row['arms'][key].items():
                raw[key + '_' + arm + '_contact'] = co.pop('mask')
        for arm in pr['arms']:
            row['continuous'][arm] = signed_gap(U, xy, z['z_' + arm], z['faces'])
        np.savez_compressed(D / (str(r['case']) + '_' + str(r['fdi']) + '_R4_reference.npz'), **raw)
        rows.append(row)
    valid = [r for r in rows if r['status'] == 'SCORED']
    summary = {}
    for kind in ['all', 'P1', 'P2', 'M1', 'M2']:
        rr = [r for r in valid if kind == 'all' or r['type'] == kind]
        s = {'n': len(rr), 'arms': {}}
        for arm in pr['arms']:
            s['arms'][arm] = {q: float(np.median([r['arms']['0.1'][arm]['errors'][q] for r in rr])) if rr else None for q in ['IoU', 'absolute_area_error_mm2', 'absolute_patch_count_error', 'centroid_error_mm']}
            s['arms'][arm]['max_penetration_mm'] = max([r['continuous'][arm]['maximum_penetration_mm'] for r in rr], default=0)
        for base in ['generic_without', 'personal_without']:
            s['personal_with_vs_' + base] = dict(median_IoU_gain=float(np.median([r['arms']['0.1']['personal_with']['errors']['IoU'] - r['arms']['0.1'][base]['errors']['IoU'] for r in rr])) if rr else None, area_error_ratio=s['arms']['personal_with']['absolute_area_error_mm2'] / max(s['arms'][base]['absolute_area_error_mm2'], 0.1) if rr else None)
        summary[kind] = s
    a = summary['all']
    contrast = a['personal_with_vs_generic_without']
    gates = dict(nonpenetration=a['arms']['personal_with']['max_penetration_mm'] <= 1e-07, IoU=contrast['median_IoU_gain'] is not None and contrast['median_IoU_gain'] >= 0.15, area=contrast['area_error_ratio'] is not None and contrast['area_error_ratio'] <= 0.8, coverage=all((summary[t]['n'] >= 6 for t in ['P1', 'P2', 'M1', 'M2'])))
    decision = 'UNKNOWN_COVERAGE' if not gates['coverage'] else 'PASS' if all(gates.values()) else 'FAIL'
    dump(H / 'raw/RESULTS_R4_ROWS.json', rows)
    dump(H / 'rounds/R4.json', dict(round='R4', claim_type='information_link', external_referent=pr['external_referent'], summary=dict(by_type=summary, gates=gates, decision=decision), cost=dict(wall_seconds=time.perf_counter() - start, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
    state('R4_DECIDED', decision, 'Whole-preparation source and independent pressure/landmark measurement')
    print(json.dumps(dict(decision=decision, gates=gates, all=summary['all']), indent=2))
if __name__ == '__main__':
    {'predict': predict, 'evaluate': evaluate}[sys.argv[1]]()
