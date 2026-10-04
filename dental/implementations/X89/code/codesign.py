"""Conditional local milling requirement. N, physical tool mm, final crown mm.

No clinical threshold, actual CAM or physical manufacturing admission.
"""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import sys, json, hashlib, math, time, datetime, resource
from pathlib import Path
import numpy as np
from scipy.integrate import quad
ROOT = Path(__file__).resolve().parents[1]
DENTAL = ROOT.parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X89-milling-codesign'))
PARENT = DENTAL / 'results/LANE_X53_MILLING_TOOLS'
sys.dont_write_bytecode = True
sys.path.insert(0, str(PARENT / 'code'))
from collision import Scene, directions, tool_capsules, EPS

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def read(p):
    return json.loads(Path(p).read_text())

def clean(x):
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {k: clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not math.isfinite(x)):
        return None
    return x

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    q.replace(p)

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def state(phase, gate, next_op):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X89-milling-codesign', claim_type='capability', phase=phase, latest_gate=gate, next_operation=next_op, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))

def lock_prereg(name='R1'):
    p = ROOT / f'PREREG_{name}.json'
    if sha(p) != p.with_suffix('.sha256').read_text().strip():
        raise ValueError('PREREG hash drift')
    d = read(p)
    for (p, h) in d.get('input_sha256', {}).items():
        if sha(p) != h:
            raise ValueError('Input drift ' + p)
    return d

def compliance(t, E):
    """Physical-tip load, clamp at gauge+radius; neck card is centre referenced.
    Cylindrical solid extending to tip is an explicit scenario, not flute inertia.
    """
    d = t['diameter_mm']
    n = t['neck_reach_mm'] + d / 2
    L = t['gauge_mm'] + d / 2
    if min(E, d, t['shank_mm'], n, L) <= 0 or n > L:
        raise ValueError('Invalid profile')
    return 64 / (3 * math.pi * E) * (n ** 3 / d ** 4 + (L ** 3 - n ** 3) / t['shank_mm'] ** 4)

def mechanical(t, m):
    c = compliance(t, m['modulus_N_mm2'])
    s = m['scale_green_to_final']
    F = m['force_N']
    e = m['other_geometry_error_final_mm']
    tol = m['film_error_tolerance_final_mm']
    eps = m['roundoff_allowance_mm']
    room = (tol - 2 * e - eps) / 2
    h = m['holder_compliance_mm_N']
    Fmax = max(0.0, room / (s * (c + h))) if c + h > 0 else float('inf')
    hmax = room / (s * F) - c if F > 0 else float('inf')
    neck = t['neck_reach_mm'] + t['diameter_mm'] / 2
    L = t['gauge_mm'] + t['diameter_mm'] / 2
    shank = 64 * (L ** 3 - neck ** 3) / (3 * math.pi * m['modulus_N_mm2'] * t['shank_mm'] ** 4)
    remaining = room / (s * F) - h - shank if F > 0 else float('inf')
    dmin = (64 * neck ** 3 / (3 * math.pi * m['modulus_N_mm2'] * remaining)) ** 0.25 if remaining > 0 else None
    u = s * F * (c + h)
    return dict(tool_compliance_mm_N=c, total_compliance_mm_N=c + h, force_max_N=Fmax, holder_compliance_max_mm_N=hmax, holder_stiffness_min_N_mm=1 / hmax if hmax > 0 else None, neck_diameter_min_mm=dmin, neck_length_physical_tip_mm=neck, tip_to_clamp_mm=L, bending_final_mm=u, mechanical_status='SCENARIO_FEASIBLE' if 2 * (u + e) + eps <= tol else 'SCENARIO_INFEASIBLE', error_budget_resolution='PER_POINT', inertia_status='UNKNOWN_REAL_FLUTES', modulus_status='TYPICAL_SCENARIO_NOT_CERTIFIED_LOWER_BOUND')

def robust_search(scene, p, n, t, m, installed_ports=None):
    """One displacement ball bounds all points of the declared translation sweep.
    Offset compensates any displacement sign, hence factor 2 in stock budget.
    """
    mech = mechanical(t, m)
    scale = m['scale_green_to_final']
    out = dict(mech)
    scaled = {k: v * scale if k in ('diameter_mm', 'neck_reach_mm', 'shank_mm', 'gauge_mm', 'holder_diameter_mm', 'holder_length_mm') else v for (k, v) in t.items()}
    r = scaled['diameter_mm'] / 2
    ds = directions(5)
    order = np.argsort(-(ds @ n), kind='stable')
    approach = float(np.linalg.norm(scene.bounds[1] - scene.bounds[0]) + 2 * scaled['holder_diameter_mm'] + 2)
    ball = scene.ball_clearance(p + (r + 5 * EPS) * n, r) >= -EPS
    out.update(ball_only=bool(ball), tool_id=t['id'], nominal_tip_diameter_mm=t['diameter_mm'], geometry_source=t['model_status'], pose=None, stock_interval_final_mm=None)
    if mech['mechanical_status'] == 'SCENARIO_INFEASIBLE' and (not installed_ports):
        out['status'] = 'MECHANICS_REJECTED'
        return out
    u = mech['bending_final_mm'] + m['other_geometry_error_final_mm']
    off = u + m['roundoff_allowance_mm']
    reasons = {}
    for k in order:
        if np.dot(n, ds[k]) < 0.001:
            continue
        selected_port = None
        if installed_ports:
            from calibration import displacement_bound
            pose_id = f"{t['id']}:{int(k)}"
            port = installed_ports.get(pose_id)
            if port is not None:
                bound = displacement_bound(port, m['force_N'], pose_id)
                if bound['status'] == 'BOUNDED_WITHIN_DECLARED_CALIBRATION':
                    selected_port = bound
                    u = scale * bound['upper_mm'] + m['other_geometry_error_final_mm']
                    off = u + m['roundoff_allowance_mm']
        if selected_port is None:
            if mech['mechanical_status'] == 'SCENARIO_INFEASIBLE':
                continue
            u = mech['bending_final_mm'] + m['other_geometry_error_final_mm']
            off = u + m['roundoff_allowance_mm']
        if 2 * u + m['roundoff_allowance_mm'] > m['film_error_tolerance_final_mm']:
            reasons['measured_response_exceeds_budget'] = reasons.get('measured_response_exceeds_budget', 0) + 1
            continue
        if u >= r:
            reasons['normal_ray_ball_support_unknown'] = reasons.get('normal_ray_ball_support_unknown', 0) + 1
            continue
        c = p + (r + off) * n
        ok = True
        for (name, a, b, rr) in tool_capsules(c, ds[k], scaled, approach):
            cl = scene.clearance(a, b, rr + u)
            if cl < -EPS:
                reasons[name] = reasons.get(name, 0) + 1
                ok = False
                break
        if ok:
            if selected_port is not None:
                out['installed_response'] = selected_port
                out['mechanical_status'] = 'BOUNDED_BY_INSTALLED_RESPONSE_PORT'
                out['nominal_inverse_requirements_status'] = 'Superseded for this pose; nominal force/holder decomposition remains scenario-only'
            out.update(status='CONDITIONAL_LOCAL_WITNESS', pose=dict(index=int(k), centre=c, direction=ds[k], offset_mm=off), stock_interval_final_mm=[off - u, off + u], rejections=reasons)
            return out
    out.update(status='NOT_FOUND_ON_DECLARED_GRID', rejections=reasons)
    return out

def intake():
    r = []
    manifest = read(PARENT / 'INPUT_MANIFEST.json')
    for x in manifest['rows']:
        if x['family'] != 'crown_loop':
            continue
        if sha(x['file']) != x['sha256']:
            raise ValueError('mesh drift')
        a = np.load(x['file'], allow_pickle=False)
        v = a['vertices']
        f = a['faces']
        inner = a['inner_faces']
        z = v[f].mean(1)[:, 2]
        rim = np.flatnonzero(z <= z.min() + 0.25)
        r.append(dict(id=x['id'], family='crown_loop', file=x['file'], sha256=x['sha256'], source=x['source'], source_sha256=x['source_sha256'], dataset=x['dataset'], license=x['license'], locator=x['locator'], regions={'intaglio': inner, 'margin_band': rim}, geometry_error_source_mm=x['geometry_error_mm'], region_status='X53 inherited heuristic intaglio; lower-z 0.25 mm margin band'))
    exports = read(DENTAL / 'results/PROOF_LANE_FULL_CROWN_R4/FROZEN_EXPORTS.json')
    for x in exports['exports']:
        if sha(x['mesh_path']) != x['mesh_sha256'] or sha(x['stl_path']) != x['stl_sha256']:
            raise ValueError('full crown drift')
        a = np.load(x['mesh_path'], allow_pickle=False)
        r.append(dict(id='full_' + x['family'], family='full_crown_R4', file=x['mesh_path'], sha256=x['mesh_sha256'], source=x['stl_path'], source_sha256=x['stl_sha256'], dataset='Bits2Bites derived designs', license='CC BY-NC-SA 4.0; confirm parent SOURCE_SELECTION', locator='https://doi.org/10.1038/s41597-025-05580-6', regions={'intaglio': np.flatnonzero(a['roles'] == 1), 'margin_band': np.flatnonzero(a['roles'] == 2)}, geometry_error_source_mm=None, region_status='Frozen parent roles: inner=1; annulus=2; original design NOT_FABRICATION_QUALIFIED'))
    return r

def samples(scene, regions, seed, n):
    rng = np.random.default_rng(seed)
    out = []
    for (region, ids) in regions.items():
        ids = np.asarray(ids, dtype=int)
        if not len(ids):
            continue
        chosen = rng.choice(ids, n, p=scene.area[ids] / scene.area[ids].sum(), replace=True)
        bary = rng.dirichlet([1, 1, 1], n)
        for (j, (i, w)) in enumerate(zip(chosen, bary)):
            tri = scene.tri[i]
            nn = np.cross(tri[1] - tri[0], tri[2] - tri[0])
            nn /= np.linalg.norm(nn)
            out.append(dict(point_index=j, face_id=int(i), region=region, point=w @ tri, normal=nn, area_weight_mm2=scene.area[ids].sum() / n))
    return out

def controls(card, m):
    q = []
    for t in card:
        n = t['neck_reach_mm'] + t['diameter_mm'] / 2
        L = t['gauge_mm'] + t['diameter_mm'] / 2
        E = m['modulus_N_mm2']
        val = quad(lambda s: s * s / (E * math.pi * t['diameter_mm'] ** 4 / 64), 0, n, epsabs=1e-13)[0] + quad(lambda s: s * s / (E * math.pi * t['shank_mm'] ** 4 / 64), n, L, epsabs=1e-13)[0]
        c = compliance(t, E)
        rel = abs(val - c) / c
        mech = mechanical(t, m)
        mut = dict(mech, force_max_N=2 * mech['force_max_N'])

        def inverse_check(x):
            return 2 * (m['scale_green_to_final'] * x['force_max_N'] * (c + m['holder_compliance_mm_N']) + m['other_geometry_error_final_mm']) + m['roundoff_allowance_mm'] <= m['film_error_tolerance_final_mm'] + 1e-10
        q.append(dict(id=t['id'], quadrature_rel_error=rel, baseline_pass=rel <= 1e-08, injected_compliance_2x_rejected=abs(val - 2 * c) / val > 1e-08, inverse_baseline_pass=inverse_check(mech), injected_force_2x_rejected=not inverse_check(mut)))
    a = dict(card[0], neck_reach_mm=3 - 0.3)
    b = dict(a, neck_reach_mm=8 - 0.3)
    summary = lambda t: [t['diameter_mm'], t['shank_mm'], t['total_length_mm'], t['gauge_mm']]
    (ca, cb) = (compliance(a, m['modulus_N_mm2']), compliance(b, m['modulus_N_mm2']))
    ident = max((abs(x - y) for (x, y) in zip(summary(a), summary(b))))
    suff = dict(summary='tip diameter, shank diameter, total length, installed gauge', summary_resolution='PER_TOOTH tool selection', state_A=a, state_B=b, identity_error_mm=ident, downstream_compliance_A_mm_N=ca, downstream_compliance_B_mm_N=cb, downstream_error_difference_final_mm=m['scale_green_to_final'] * m['force_N'] * abs(cb - ca), status='SUMMARY_INSUFFICIENT', minimum_extension='Bending compliance integral (or measured installed tip compliance), locally linked to force and pose', rigorous_floating_enclosure=False)
    impossible = [dict(m, force_N=100), dict(m, holder_compliance_mm_N=1), dict(m, other_geometry_error_final_mm=0.02)]
    bad = [all((mechanical(t, mm)['mechanical_status'] == 'SCENARIO_INFEASIBLE' for t in card)) for mm in impossible]
    return dict(quadrature_and_inverse=q, sufficiency=suff, injected_impossible=dict(high_force_rejected=bad[0], soft_holder_rejected=bad[1], exhausted_error_budget_rejected=bad[2]))

def run(replay=False):
    pre = lock_prereg()
    m = pre['metrics']
    start = time.perf_counter()
    t0 = time.perf_counter()
    rows = intake()
    prep_time = time.perf_counter() - t0
    card = read(PARENT / 'NOMINAL_TIP_REFERENCE_CARD.json')['tools']
    summary = []
    hashes = {}
    for (ix, row) in enumerate(rows):
        state('R1_RUNNING', f'{len(summary)}/{len(rows)} designs processed', 'Pointwise coupled access and inverse stiffness')
        a = np.load(row['file'], allow_pickle=False)
        try:
            scene = Scene(a['vertices'], a['faces'])
        except ValueError as ex:
            summary.append(dict(id=row['id'], family=row['family'], points=0, found=0, ball_only=0, regions={}, status='INPUT_REJECTED', reason=str(ex), physical_status='UNKNOWN_INPUT_REJECTED'))
            print(row['id'], 'INPUT_REJECTED', str(ex), flush=True)
            continue
        pts = samples(scene, row['regions'], m['seed'] + ix, m['n_per_region'])
        rr = []
        for p in pts:
            tools = [robust_search(scene, p['point'], p['normal'], t, m) for t in card]
            viable = [t for t in tools if t['status'] == 'CONDITIONAL_LOCAL_WITNESS']
            rr.append(dict(**p, tools=tools, selected_tool=min(viable, key=lambda x: x['nominal_tip_diameter_mm'])['tool_id'] if viable else None, status='CONDITIONAL_LOCAL_WITNESS' if viable else 'NO_LIBRARY_WITNESS_AT_DECLARED_LOAD', resolution='PER_POINT', time_scale='HANDOVER'))
        out = DATA / 'R1' / f"{row['id']}.json"
        dump(out, dict(input={k: v for (k, v) in row.items() if k != 'regions'}, prereg_sha256=sha(ROOT / 'PREREG_R1.json'), rows=rr))
        hashes[str(out)] = sha(out)
        ss = dict(id=row['id'], family=row['family'], points=len(rr), found=sum((x['selected_tool'] is not None for x in rr)), ball_only=sum((any((t['ball_only'] for t in x['tools'])) for x in rr)), physical_status='UNKNOWN_MISSING_ASSEMBLY_FORCE_CAM_PROCESS_AND_WHOLE_SURFACE_BOUND', raw_file=out, regions={reg: dict(n=sum((x['region'] == reg for x in rr)), found=sum((x['region'] == reg and x['selected_tool'] is not None for x in rr))) for reg in row['regions']})
        summary.append(ss)
        print(row['id'], ss['found'], '/', ss['points'], flush=True)
        dump(ROOT / 'raw/R1_PARTIAL.json', summary)
    v0 = time.perf_counter()
    checks = controls(card, m)
    validation_s = time.perf_counter() - v0
    dump(ROOT / 'raw/CONTROLS_R1.json', checks)
    assert all((all((x[k] for k in ('baseline_pass', 'injected_compliance_2x_rejected', 'inverse_baseline_pass', 'injected_force_2x_rejected'))) for x in checks['quadrature_and_inverse']))
    assert all(checks['injected_impossible'].values())
    assert checks['sufficiency']['identity_error_mm'] == 0 and checks['sufficiency']['downstream_error_difference_final_mm'] >= 0.001
    result = dict(claim_type='capability', round='R1', status='CONDITIONAL_INVERSE_REQUIREMENTS; PHYSICAL_UNKNOWN', prereg_sha256=sha(ROOT / 'PREREG_R1.json'), metrics=m, designs=summary, tool_mechanics=[dict(id=t['id'], **mechanical(t, m)) for t in card], controls=checks, artifact_sha256=hashes, external_referent={'kind': 'closed_form', 'locator': 'Euler-Bernoulli cantilever, XBREAK_HUNT_3/code/beam.py and independent SciPy quad', 'compared_quantity': 'Declared circular static beam tip compliance [mm/N]', 'refutes_us': False}, published_measurement_validation='PENDING_SOURCE_SCOPE_R2', resolution='PER_POINT', time_scale='HANDOVER', dropout=dict(design_rejected=sum((x.get('status') == 'INPUT_REJECTED' for x in summary)), requested=8, physical_qualifications_missing=8), cost=dict(preparation_s=prep_time, query_and_discovery_s=time.perf_counter() - start - validation_s - prep_time, validation_s=validation_s, total_wall_s=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, fit_s=0, questions_s=0, fallback_s=None, prior_search_preparation_s=None, gpu=0), rigorous_analytic_enclosure='Triangle-inequality displacement tube under declared quasistatic solid-beam assumptions', rigorous_floating_enclosure=False, physical_uncertainty='UNBOUNDED because required measurement leaves are missing')
    if replay:
        previous = read(ROOT / 'RESULTS_R1.json')
        assert clean(result['designs']) == previous['designs'], 'Frozen R1 outcome drift'
        assert clean(result['controls']) == previous['controls'], 'Frozen R1 controls drift'
        dump(ROOT / 'raw/R1_REPLAY.json', result)
    else:
        dump(ROOT / 'RESULTS_R1.json', result)
        dump(ROOT / 'results.json', result)
    state('R1_COMPLETE', 'Inverse algebra, quadrature, exact sufficiency pair and impossible mechanics injections PASS', 'R2: external milling error scope; measured assembly port and X49 CLI integration')
    return result
if __name__ == '__main__':
    run(replay='--replay' in sys.argv)
