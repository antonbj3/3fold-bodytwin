"""Executed independent controls and falsifiable contract tests, never cached booleans."""
import copy
import math
import numpy as np
from implant_safety import ContractError
from implant_safety.module import classify, guide_components
from implant_safety.vendor.x8_geometry import project_cylinder, voxel_cylinder_bracket
from implant_safety.vendor.proof_lane_control import control_box
from implant_safety.vendor.x87_science import controls as x87_controls, guide as x87_guide
from implant_safety.vendor.x96_guide import scalar_control
from scipy.optimize import minimize

def check_distance(interval, control_mm, tolerance=1e-06):
    if not (math.isfinite(control_mm) and interval[0] - tolerance <= control_mm <= interval[1] + tolerance):
        raise AssertionError('independent distance outside candidate bracket')

def independent_geometry(points, entry, axis, length, radius, spacing):
    d = np.linalg.norm(points - project_cylinder(points, entry, axis, length, radius), axis=1)
    indices = np.flatnonzero(d - np.linalg.norm(spacing / 2) <= float(d.min()) + 1e-08)
    values = [control_box(points[i], entry, axis, length, radius, spacing) for i in indices]
    return (min((v[0] for v in values)), all((v[1] for v in values)), len(indices))

def control_box_enclosure(center, entry, axis, L, r, sp):
    """PROOF_LANE's same independent primal variables, with X8 convex support bounds.

    This records real-arithmetic bounds, not rigorous IEEE enclosures. It does
    not require success of irrelevant optimizer subproblems to certify a gap.
    """
    lo = center - sp / 2
    hi = center + sp / 2
    v = np.cross(axis, [1.0, 0, 0] if abs(axis[0]) < 0.9 else [0, 1.0, 0])
    v /= np.linalg.norm(v)
    w = np.cross(axis, v)

    def xyz(z):
        return entry + z[3] * axis + z[4] * v + z[5] * w

    def fun(z):
        dv = z[:3] - xyz(z)
        return dv @ dv

    def jac(z):
        dv = z[:3] - xyz(z)
        return np.r_[2 * dv, -2 * dv @ axis, -2 * dv @ v, -2 * dv @ w]
    t = np.clip((center - entry) @ axis, 0, L)
    uv = np.array([(center - entry) @ v, (center - entry) @ w])
    uv *= min(1, r / max(np.linalg.norm(uv), 1e-30))
    res = minimize(fun, np.r_[center, t, uv], jac=jac, bounds=list(zip(lo, hi)) + [(0, L), (-r, r), (-r, r)], constraints=[{'type': 'ineq', 'fun': lambda z: r * r - z[4] ** 2 - z[5] ** 2, 'jac': lambda z: np.array([0.0, 0, 0, 0, -2 * z[4], -2 * z[5]])}], method='SLSQP', options={'ftol': 1e-12, 'maxiter': 200})
    z = res.x.copy()
    z[:3] = np.clip(z[:3], lo, hi)
    z[3] = np.clip(z[3], 0, L)
    z[4:6] *= min(1, r / max(np.linalg.norm(z[4:6]), 1e-30))
    distance = float(np.linalg.norm(z[:3] - xyz(z)))
    if distance == 0:
        lower = 0.0
    else:
        n = (z[:3] - xyz(z)) / distance
        box_corner = np.where(n >= 0, lo, hi)
        lower = max(0.0, float(n @ (box_corner - entry) - max(0.0, L * (n @ axis)) - r * math.hypot(n @ v, n @ w)))
    assert lower <= distance + 1e-10
    return dict(lower_mm=lower, upper_mm=distance, solver_success=bool(res.success), primal_feasibility_error=max(float(np.max(lo - z[:3])), float(np.max(z[:3] - hi)), -float(z[3]), float(z[3] - L), float(np.linalg.norm(z[4:6]) - r), 0.0))

def independent_enclosure(points, entry, axis, length, radius, spacing):
    d = np.linalg.norm(points - project_cylinder(points, entry, axis, length, radius), axis=1)
    lower_centers = np.maximum(0, d - np.linalg.norm(spacing / 2))
    indices = np.flatnonzero(lower_centers <= float(d.min()) + 1e-08)
    rec = [control_box_enclosure(points[i], entry, axis, length, radius, spacing) for i in indices]
    discarded = np.ones(len(points), dtype=bool)
    discarded[indices] = False
    lower = min([x['lower_mm'] for x in rec] + ([float(lower_centers[discarded].min())] if discarded.any() else []))
    upper = min((x['upper_mm'] for x in rec))
    assert upper - lower <= 1e-06 and upper - lower >= -1e-10
    assert max((x['primal_feasibility_error'] for x in rec)) <= 1e-10
    return dict(interval_mm=[lower, upper], width_mm=upper - lower, all_solver_success=all((x['solver_success'] for x in rec)), primal_feasibility_error=max((x['primal_feasibility_error'] for x in rec)), active_boxes=len(indices), rigorous_IEEE_enclosure='MISSING')

def controls(module, rows):
    seen = set()
    geom = []
    guide = []
    for r in rows:
        q = r['query']
        site = module.sites[r['site_id']]
        (nominal, union, _) = module.points(site)
        e = np.array(q['cbct_pose']['entry_zyx_mm'])
        a = np.array(q['cbct_pose']['axis_zyx'])
        L = q['length_mm']
        rad = q['radius_mm']
        h = r['terms']['drill_extra_depth']['maximum_or_example_mm']
        for (label, pts, l) in [('nominal', nominal, L), ('revision_union', union, L), ('tool_envelope', union, L + h)]:
            k = (r['site_id'], label, l, rad)
            if k in seen:
                continue
            seen.add(k)
            candidate = r['digital_geometry'][label]['interval_mm']
            enclosure = independent_enclosure(pts, e, a, l, rad, np.array(site['frame']['spacing_mm']))
            ctrl = enclosure['interval_mm'][1]
            success = enclosure['all_solver_success']
            n = enclosure['active_boxes']
            check_distance(candidate, ctrl)
            injected = False
            try:
                check_distance(candidate, ctrl + 1)
            except AssertionError:
                injected = True
            assert injected
            geom.append(dict(site_id=r['site_id'], quantity=label, length_mm=l, candidate_interval_mm=candidate, control_mm=ctrl, error_mm=max(candidate[0] - ctrl, 0, ctrl - candidate[1]), optimizer_success=success, active_boxes=n, injected_plus1mm_rejected=injected))
            geom[-1]['independent_control_enclosure'] = enclosure
    for (gid, p) in module.profiles.items():
        g = x87_guide(p, 0.05)
        c = x87_controls(p, 0.05, g)
        guide.append(c)
        for h in [0.5, 1.0, 1.3]:
            gg = guide_components(p, 0.95, 2.0, h)
            th = math.radians(gg['angle_deg'])
            r = 2 + h
            rot = math.hypot(r * math.cos(min(th, math.pi)) - r, r * math.sin(min(th, math.pi)))
            expected = max(gg['entry_mm'], gg['apex_mm']) + rot
            err = abs(expected - gg['combined_displacement_mm'])
            assert err <= 1e-10
            assert abs(expected - (gg['combined_displacement_mm'] + 1)) > 1e-10
            tail_bound = float(scalar_control(2 + gg['combined_displacement_mm'] + 1e-08, p, h))
            assert tail_bound <= 0.05 + 1e-10
            guide.append(dict(guide=gid, extra_mm=h, error_mm=err, independent_tail_bound=tail_bound, injected_plus1mm_rejected=True))
    return dict(geometry=geom, guide=guide, control_gate='PASS', geometry_max_error_mm=max((x['error_mm'] for x in geom)), guide_max_error_mm=max((x['error_mm'] for x in guide)), rigorous_IEEE_enclosure='MISSING')

def sufficiency():
    e = np.zeros(3)
    a = np.array([0.0, 0.0, 1.0])
    sp = np.full(3, 0.5)
    axial = np.array([[0.0, 0.0, 6.75]])
    lateral = np.array([[3.75, 0.0, 2.0]])
    b = [voxel_cylinder_bracket(p, sp, e, a, 4.0, 1.0) for p in [axial, lateral]]
    summary = [x['upper_mm'] for x in b]
    assert summary[0] == summary[1] == 2.5
    downstream = [voxel_cylinder_bracket(p, sp, e, a, 5.0, 1.0)['upper_mm'] for p in [axial, lateral]]
    assert downstream == [1.5, 2.5]
    ctrl = [independent_geometry(p, e, a, 5.0, 1.0, sp)[0] for p in [axial, lateral]]
    for (d, c) in zip(downstream, ctrl):
        check_distance([d, d], c)
    return dict(summary='nominal clearance', identical_summary_mm=summary, identity_error_mm=abs(summary[0] - summary[1]), downstream_swept_gap_mm=downstream, downstream_difference_mm=1.0, decision_vs2mm=[classify([v, v]) for v in downstream], outcome='SCALAR_CLEARANCE_INSUFFICIENT', minimum_extension='For this frozen sweep: signed tool-induced gap loss. For changed pose or length: wall location and direction/occupancy geometry.', external_referent=dict(kind='our_own_fixture', locator='validation.py:sufficiency', compared_quantity='logical exact-summary counterexample', refutes_us=True), resolution='PER_TOOTH', source_scope='mathematical witness; no anatomical validation')

def contract_tests(module, q):
    checks = []

    def reject(name, bad):
        try:
            module.query(**bad)
        except ContractError as e:
            checks.append(dict(name=name, result='PASS', rejection=str(e)))
            return
        raise AssertionError('injected contract error accepted: ' + name)
    for (key, val) in [('units', 'cm'), ('order', 'xyz'), ('spacing_mm', [1.0, 1.0, 1.0]), ('origin_mm', [0.0, 0.0, 1.0]), ('direction', [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])]:
        bad = copy.deepcopy(q)
        bad['cbct_pose']['frame'][key] = val
        reject('frame_' + key, bad)
    for (name, field, val) in [('zero_axis', 'axis_zyx', [0.0, 0.0, 0.0]), ('nan_entry', 'entry_zyx_mm', [0.0, 0.0, float('nan')])]:
        bad = copy.deepcopy(q)
        bad['cbct_pose'][field] = val
        reject(name, bad)
    for (key, val) in [('length_mm', -1.0), ('length_mm', float('nan')), ('radius_mm', 3.0), ('confidence', 0.99)]:
        bad = copy.deepcopy(q)
        bad[key] = val
        reject('invalid_' + key + str(val), bad)
    cases = [('unknown_system', dict(implant_system='UNREVIEWED_SYSTEM'), 'protocol'), ('unsupported_t3_length', dict(implant_system='ZimVie_T3_parallel_ACT3p85_datum_aligned'), 'protocol'), ('label_is_not_actual', dict(length_reference='manufacturer_label'), 'protocol'), ('unknown_guide', dict(guide_type='UNREVIEWED_GUIDE'), 'guide')]
    for (name, change, kind) in cases:
        bad = copy.deepcopy(q)
        bad.update(change)
        r = module.query(**bad)
        assert r['combined_conditional']['interval_mm'] is None
        assert r['physical_safety']['status'] == 'UNKNOWN'
        checks.append(dict(name=name, result='PASS', outcome='UNKNOWN retained'))
    bad = copy.deepcopy(q)
    bad['length_mm'] = 6.0
    r = module.query(**bad)
    assert r['bone_observation']['status'] == 'UNKNOWN_CHANGED_POSE_OR_LENGTH'
    assert r['query']['length_mm'] == 6.0
    checks.append(dict(name='changed_length_recomputes_and_detaches_bone', result='PASS'))
    bad = copy.deepcopy(q)
    bad['implant_system'] = 'ZimVie_T3_parallel_ACT3p85_datum_aligned'
    bad['length_mm'] = 12.6
    bad['protocol_datum'] = dict(fixture_label_mm=13.0, actual_implant_mm=12.6, actual_mark_mm=13.7, platform_depth_mm=1.0, tip_mm=1.2)
    valid = module.query(**bad)
    assert valid['terms']['drill_extra_depth']['maximum_or_example_mm'] == 1.3
    checks.append(dict(name='restricted_t3_exact_datum_supported', result='PASS'))
    bad['protocol_datum']['platform_depth_mm'] = 0.0
    assert module.query(**bad)['combined_conditional']['interval_mm'] is None
    checks.append(dict(name='wrong_t3_platform_datum', result='PASS', outcome='UNKNOWN retained'))
    r = module.query(**q)
    g = r['terms']['guide_combination']
    assert abs(g['combined_displacement_mm'] - (max(g['entry_mm'], g['apex_mm']) + g['rotation_mm'])) <= 1e-10
    assert g['combined_displacement_mm'] < g['entry_mm'] + g['apex_mm'] + g['rotation_mm']
    assert r['physical_safety']['interval_mm'] is None
    assert r['physical_safety']['injury_probability'] is None
    checks.append(dict(name='overlap_not_summed_and_physical_unknown', result='PASS'))
    return checks
