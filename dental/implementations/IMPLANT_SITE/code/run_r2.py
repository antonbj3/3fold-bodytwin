from common_local import *
import time, resource, itertools
from scipy.optimize import minimize
start = time.perf_counter()
check_frozen()
pr = load(P / 'PREREG_R2.json')
DATA.mkdir(parents=True, exist_ok=True)
geom = module('x8_geometry', X8 / 'full_geometry.py')
guide = module('x87_budget', X87 / 'science.py')
register_inputs([X8 / 'full_geometry.py'])
keys = {(r['case'], r['fdi']) for r in load(P / 'raw/DENSE_BONE_INTERSECTION.json')}
sites = sorted([s for s in load(X73 / 'raw/LINEAGE_SITES.json') if (s['case'], s['fdi']) in keys], key=lambda s: (s['case'], s['fdi']))
bone = {(r['case'], r['fdi']): r for r in lines(X75 / 'raw/PER_SITE_R2.jsonl') if r['kind'].startswith('X8')}
profiles = load(X87 / 'inputs/guide_profiles.json')['profiles']
out = []
detail = []
control = []
replay = []
plans = []
from control_geometry import control_box

def pareto(rows):
    return [r for r in rows if not any((q['L_mm'] >= r['L_mm'] and q['D_mm'] >= r['D_mm'] and (q['L_mm'] > r['L_mm'] or q['D_mm'] > r['D_mm']) for q in rows))]
for (i, s) in enumerate(sites):
    assert not s.get('alias_source_bindings'), 'R2 needs explicit alias-mask pool; cannot omit aliases'
    ar = s['point_artifact']
    register_inputs([ar['path']])
    assert sha(ar['path']) == ar['sha256']
    with np.load(ar['path']) as a:
        pts = {k: a[k].astype(float) * 0.3 for k in a.files if k.endswith('_voxels_zyx')}
    e = np.array(s['pose']['entry_zyx_mm'])
    axis = np.array(s['pose']['axis_zyx'])
    b = bone[s['case'], s['fdi']]
    signed = [v for ray in b['rays'] for v in [ray['buccal_extent_mm'], ray['lingual_extent_mm']]]
    ax = b['axial_run_limits_from_pose_mm']
    local = []
    menu = list(itertools.product(pr['selection']['length_menu_mm'], pr['selection']['diameter_menu_mm']))
    menu += [(s['pose']['length_mm'], 2 * s['pose']['radius_mm'])]
    for (j, (L, diam)) in enumerate(menu):
        ds = {}
        cs = []
        for (mask, points) in pts.items():
            rr = geom.voxel_cylinder_bracket(points, np.full(3, 0.3), e, axis, L, diam / 2)
            ds[mask] = {'lower_mm': rr['lower_mm'], 'upper_mm': rr['upper_mm'], 'gap_mm': rr['gap_mm'], 'witness': rr['witness'], 'active_voxel_boxes': rr['active_voxel_boxes']}
            rel = points - e
            t = rel @ axis
            rad = np.sqrt(np.maximum(0, np.sum(rel * rel, 1) - t * t))
            dz = np.maximum(0, np.maximum(-t, t - L))
            dr = np.maximum(0, rad - diam / 2)
            dist = np.hypot(dz, dr)
            prune = float(np.min(dist)) - 0.3 * np.sqrt(3) / 2
            assert rr['lower_mm'] >= max(0, prune) - 1e-07
            if L in [6.0, 14.0] and diam in [3.0, 4.0] or j == len(menu) - 1:
                active = np.where(dist - 0.3 * np.sqrt(3) / 2 <= rr['upper_mm'] + 1e-06)[0]
                vals = [control_box(points[k], e, axis, L, diam / 2, np.full(3, 0.3)) for k in active]
                cv = min((v[0] for v in vals))
                err = max(0, rr['lower_mm'] - cv, cv - rr['upper_mm'])
                cs.append({'mask': mask, 'value_mm': cv, 'error_mm': err, 'optimizer_success_all': all((v[1] for v in vals)), 'boxes': len(active)})
                control.append(dict(case=s['case'], fdi=s['fdi'], L_mm=L, D_mm=diam, **cs[-1]))
        lo = min((v['lower_mm'] for v in ds.values()))
        hi = min((v['upper_mm'] for v in ds.values()))
        width = hi - lo
        if j == len(menu) - 1:
            replay.append(dict(case=s['case'], fdi=s['fdi'], max_error_mm=max(abs(lo - s['union']['lower_mm']), abs(hi - s['union']['upper_mm']))))
            continue
        support = all((v is not None and v >= diam / 2 for v in signed)) and ax is not None and (ax[0] <= 0) and (L <= ax[1])
        core = dict(case=s['case'], fdi=s['fdi'], image_group=s['image_group'], L_mm=L, D_mm=diam, entry_zyx_mm=e.tolist(), axis_zyx=axis.tolist(), segment=s['segment'], union_lower_mm=lo, union_upper_mm=hi, bracket_width_mm=width, sampled_envelope_pass=support, signed_side_min_minus_radius_mm=min(signed) - diam / 2, axial_support_margin_mm=min(-ax[0], ax[1] - L), fixed_2mm_tf2=ds['tf2_voxels_zyx']['lower_mm'] >= 2, dense_2mm_pass=lo >= 2, geometry_valid=width <= 1e-06, resolution='PER_TOOTH', time_scale='SIMULTANEOUS')
        local.append(core)
        detail.append(dict(**core, mask_distances=ds, source_point_artifact=ar))
        for p in profiles:
            for target in [0.9, 0.95]:
                budget = guide.guide(p, 1 - target)['guide_budget_mm']
                out.append(core | dict(guide=p['id'], target=target, conditional_guide_budget_mm=budget, nerve_clearance_margin_mm=lo - budget, conditional_nerve_pass=lo >= budget, joint_digital_pass=lo >= budget and support and (width <= 1e-06), physical_joint_plan='UNKNOWN', evidence='MODELED_DIGITAL_SCREEN_ONLY'))
    for p in profiles:
        for target in [0.9, 0.95]:
            rr = [r for r in out if r['case'] == s['case'] and r['fdi'] == s['fdi'] and (r['guide'] == p['id']) and (r['target'] == target)]
            eligible = [r for r in rr if r['joint_digital_pass']]
            pf = pareto(eligible)
            selected = max(eligible, key=lambda r: (r['L_mm'], r['D_mm'])) if eligible else None
            old = [r for r in rr if r['fixed_2mm_tf2'] and r['sampled_envelope_pass']]
            oldsel = max(old, key=lambda r: (r['L_mm'], r['D_mm'])) if old else None
            keep = lambda r: {k: r[k] for k in ['L_mm', 'D_mm', 'union_lower_mm', 'union_upper_mm', 'nerve_clearance_margin_mm', 'signed_side_min_minus_radius_mm', 'axial_support_margin_mm']} if r else None
            plans.append(dict(case=s['case'], fdi=s['fdi'], guide=p['id'], target=target, finite_research_choice=keep(selected), pareto_candidates=[keep(r) for r in pf], fixed_2mm_research_choice=keep(oldsel), physical_recommendation=None, objective='Prescribed research display: maximum L, then D; not validated utility', regional_support='Necessary sampled label-envelope containment, includes teeth and is not cortical/bone adequacy', bone_gray_status=b['trabecular_candidate_roi']['status']))
    for (a, c) in itertools.product(local, repeat=2):
        if a['L_mm'] <= c['L_mm'] and a['D_mm'] <= c['D_mm']:
            assert a['union_upper_mm'] + 1e-06 >= c['union_lower_mm']
    print(f"{i + 1}/{len(sites)} {s['case']} FDI{s['fdi']} complete", flush=True)
    state('R2_MEASURING', 'PARTIAL', f'Continue {len(sites) - i - 1} sites', completed_sites=i + 1, candidates=len(detail))
e = np.zeros(3)
a = np.array([1.0, 0, 0])
obstacles = np.array([[12.0, 0, 0], [5.0, 4.0, 0]])

def gap(L):
    return np.linalg.norm(obstacles - geom.project_cylinder(obstacles, e, a, L, 2.0), axis=1)
g0 = gap(10.0)
g1 = gap(8.0)
assert g0[0] == g0[1] == 2.0
witness = dict(nominal_gap_mm=g0.tolist(), identity_error=0.0, shortened_gap_mm=g1.tolist(), downstream_difference_mm=float(g1[0] - g1[1]), minimum_extension='For this pair nearest-wall orientation/location; general search requires full spatial canal occupancy and pose', resolution='PHENOMENOLOGICAL', external_referent={'kind': 'our_own_fixture', 'locator': 'code/run_r2.py', 'compared_quantity': 'Exact set-distance counterexample', 'refutes_us': True})
csvout(P / 'raw/CANDIDATE_GUIDE_R2.csv', out)
dump(P / 'raw/RESEARCH_CHOICES_R2.json', plans)
dump(P / 'raw/GEOMETRY_CONTROLS_R2.json', control)
dump(P / 'raw/SUFFICIENCY_R2.json', witness)
raw = DATA / 'CANDIDATE_GEOMETRY_R2.jsonl'
raw.write_text(''.join((json.dumps(r, allow_nan=False) + '\n' for r in detail)))
res = dict(claim_type=['information_link', 'capability'], status='FINITE_DIGITAL_DIMENSION_QUERY_PHYSICAL_PLAN_UNKNOWN', sites=len(sites), cases=len({s['case'] for s in sites}), dimension_candidates=len(detail), guide_queries=len(out), max_bracket_width_mm=max((s['bracket_width_mm'] for s in detail)), max_control_error_mm=max((c['error_mm'] for c in control)), max_original_replay_error_mm=max((r['max_error_mm'] for r in replay)), gates={'bracket': all((s['bracket_width_mm'] <= 1e-06 for s in detail)), 'control': all((c['error_mm'] <= 1e-06 for c in control)), 'original_replay': all((r['max_error_mm'] <= 1e-07 for r in replay))}, summary=[dict(guide=p['id'], target=t, sites_with_choice=sum((r['finite_research_choice'] is not None for r in plans if r['guide'] == p['id'] and r['target'] == t)), fixed2_sites_with_choice=sum((r['fixed_2mm_research_choice'] is not None for r in plans if r['guide'] == p['id'] and r['target'] == t))) for p in profiles for t in [0.9, 0.95]], digital_envelope_rejected_candidates=sum((not r['sampled_envelope_pass'] for r in detail)), physical_recommendations=0, raw_artifact=dict(path=str(raw), sha256=sha(raw), bytes=raw.stat().st_size), external_referent=pr['external_referent'], cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, threads=1, gpu=False))
dump(P / 'rounds/R2.json', res)
dump(P / 'raw/ORIGINAL_REPLAY_R2.json', replay)
state('R2_COMPLETE', res['gates'], 'Freeze next construction: directly measure signed clearance loss with registered achieved pose; export test case')
(P / 'HANDOFF_R2.md').write_text("# R2 : digital dimension issue\n\n" + json.dumps(res, indent=2) + "\n\nThe next design needs to change the information about the guide's direction: the same major radial error can go towards or from the channel. Get registered achieved pose at the same channel region, not another population margin.\n")
print(json.dumps(res, indent=2))
