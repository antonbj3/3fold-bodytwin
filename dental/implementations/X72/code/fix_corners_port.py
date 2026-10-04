from dental_release.paths import expand as _release_expand
from common import *
import math, itertools, copy, time
import numpy as np
from field_engine.contact_port_v1 import *
from motion_engine.ncp.contact_port_adapter import solver_gap
from motion_engine.ncp.latent_contact import make_family, certify_family
from motion_engine.ncp.family_sweep import certify_sweep, verify_sweep

def qbounds(sp, box, a):
    sp = list(map(Q, sp))
    box = [list(map(Q, r)) for r in box]
    a = list(map(Q, a))
    assert min(a) > 0
    lo = [s + r[0] for (s, r) in zip(sp, box)]
    hi = [s + r[1] for (s, r) in zip(sp, box)]
    gl = [max([Q(0), lo[j]] + [lo[j] - a[j] * hi[i] / a[i] for i in range(len(a)) if i != j]) for j in range(len(a))]
    gu = [max([Q(0), hi[j]] + [hi[j] - a[j] * lo[i] / a[i] for i in range(len(a)) if i != j]) for j in range(len(a))]
    return dict(gap_lower_um=gl, gap_upper_um=gu, lift_interval_um=[max(Q(0), max((-v / w for (v, w) in zip(hi, a)))), max(Q(0), max((-v / w for (v, w) in zip(lo, a))))])

def replay(d):
    if d['input_source_sha256'] != sha(BASE / 'LANE_X10_SCANNER_GAP/PREREG_R4.json'):
        return False
    ans = qbounds(d['spacers_um'], d['box_um'], d['axis_projections'])
    return enc(ans) == d['bounds'] and d['full_contact_scope'] is True and (d['same_datum'] is True) and (d['process_bound_known'] is True)

def replay_dental_port(port):
    if not replay_witnesses(port.branch_set, {'dental_planar_seating_v1': replay}):
        return False
    w = port.branch_set.branches[0].witnesses[0]
    d = json.loads(w.payload_json)
    point = all((Q(lo) == Q(hi) for (lo, hi) in d['box_um']))
    if port.branch_set.complete != point:
        return False
    if len(port.contacts) != len(d['bounds']['gap_lower_um']):
        return False
    for (j, contact) in enumerate(port.contacts):
        gap = contact.gap.to(Unit.MM)
        if contact.frame != 'fixture_physical_v1' or gap.source != w.source or gap.lower != Q(d['bounds']['gap_lower_um'][j]) / 1000 or (gap.upper != Q(d['bounds']['gap_upper_um'][j]) / 1000):
            return False
    return port.physical_status == Status.UNKNOWN

def stl(path):
    b = record(path).read_bytes()
    n = int.from_bytes(b[80:84], 'little')
    dt = np.dtype([('normal', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
    assert len(b) == 84 + 50 * n
    return np.frombuffer(b, dtype=dt, offset=84)['v'].astype(float)

def main():
    t = time.perf_counter()
    parent = BASE / _release_expand('X10')
    sg = module('x72_signed_gap', parent / 'code/signed_gap.py')
    pr = read(parent / 'PREREG_R4.json')
    old = read(parent / 'RESULTS_R4.json')
    a = [Q(1), Q(math.sin(math.radians(pr['scenario']['half_taper_deg']))), Q(1)]
    r = Q(str(pr['scenario']['same_rms_test_um']))
    sp = [40, 40, 60]
    cases = {'same_rms_expansion': [[0, 0], [r, r], [0, 0]], 'same_rms_contraction': [[0, 0], [-r, -r], [0, 0]], 'near_threshold_reference': [[-5, 5], [-50, -40], [0, 0]], 'signed_measured_like_minus50': [[-5, 5], [-55, -45], [0, 0]], 'signed_measured_like_minus20': [[-5, 5], [-25, -15], [0, 0]], 'bridge_far_contraction': [[0, 0], [r, r], [0, 0], [0, 0], [-r, -r], [0, 0]]}
    port_rows = []
    corner_rows = []
    mutations = []
    for (name, box) in cases.items():
        aa = a * 2 if name.startswith('bridge') else a
        ss = sp * 2 if name.startswith('bridge') else sp
        bounds = qbounds(ss, box, aa)
        ctrl = sg.vertex_control(list(map(float, ss)), [[float(v) for v in r] for r in box], list(map(float, aa)))
        err = max((abs(float(v) - w) for k in bounds for (v, w) in zip(bounds[k], ctrl[k])))
        assert err <= 1e-08
        original = old['case_bounds'][name]
        olderr = max((abs(float(v) - w) for k in bounds for (v, w) in zip(bounds[k], original[k])))
        assert olderr <= 1e-08
        source = Source(name, sha(parent / 'PREREG_R4.json'), 'conditional declared planar-wall case; not physical measurement')
        payload = enc(dict(input_source_sha256=source.sha256, spacers_um=ss, box_um=box, axis_projections=aa, bounds=bounds, full_contact_scope=True, same_datum=True, process_bound_known=True))
        w = Witness('dental_planar_seating_v1', canonical(payload), digest(payload), source)
        point = all((lo == hi for (lo, hi) in box))
        bs = BranchSet((Branch('point-seating' if point else 'enclosed-seating-family', canonical(dict(lift_interval_um=bounds['lift_interval_um'])), (w,)),), point, 'point minimum lift is unique; interval-family state cardinality not enumerated', w if point else None, 'BOX_ENCLOSURE_IS_NOT_ONE_ENUMERATED_EQUILIBRIUM' if not point else '')
        contacts = tuple((ContactCandidate('wall:' + str(j), ('crown', 'die'), 'fixture_physical_v1', Gap(l / 1000, u / 1000, Unit.MM, source, sigma=None, assurance=Assurance.CONDITIONAL, condition='full-wall signed bounds, common physical datum and known manufacturing bound')) for (j, (l, u)) in enumerate(zip(bounds['gap_lower_um'], bounds['gap_upper_um']))))
        port = Port('X10/X72', 'Field+Motion contact v1', contacts, bs)
        decoded = Port.from_json(port.to_json())
        assert decoded.to_json() == port.to_json()
        assert replay_dental_port(decoded)
        from dataclasses import replace
        wrong_gap = replace(decoded.contacts[0].gap, upper=decoded.contacts[0].gap.upper + 1)
        wrong_port = replace(decoded, contacts=(replace(decoded.contacts[0], gap=wrong_gap), *decoded.contacts[1:]))
        assert not replay_dental_port(wrong_port)
        mutations.append(dict(case=name, error='wire_gap_unbound_to_witness', rejected=True))
        if not point:
            wrong_bs = replace(decoded.branch_set, complete=True, coverage=w)
            assert not replay_dental_port(replace(decoded, branch_set=wrong_bs))
            mutations.append(dict(case=name, error='false_unique_box_completeness', rejected=True))
        req = solver_gap(decoded.contacts[0].gap, dt_s='1/100', convention='CLOTH_SIGNED')
        assert req.normal_offset_m_s == (bounds['gap_lower_um'][0] / 1000000 * 100, bounds['gap_upper_um'][0] / 1000000 * 100) and req.sigma_m_s is None
        lim = Q(120)
        legacy_status = 'ACCEPT_CONDITIONAL_BOX' if bounds['gap_upper_um'][0] < lim else 'REJECT_CONDITIONAL_BOX' if bounds['gap_lower_um'][0] >= lim else 'ABSTAIN_STRADDLES_BOUNDARY'
        v1 = decoded.contacts[0].gap.classify(Q(120, 1000))
        for (key, val) in [('bounds', None), ('same_datum', False), ('process_bound_known', False), ('full_contact_scope', False)]:
            bad = copy.deepcopy(payload)
            if key == 'bounds':
                bad['bounds']['gap_upper_um'][0] = str(Q(bad['bounds']['gap_upper_um'][0]) + 1)
            else:
                bad[key] = val
            bw = Witness(w.kind, canonical(bad), digest(bad), source)
            bbs = BranchSet((Branch('tampered', canonical({}), (bw,)),), True, bs.scope, bw)
            rejected = not replay_witnesses(bbs, {w.kind: replay})
            assert rejected
            mutations.append(dict(case=name, error=key, rejected=rejected))
        write(ROOT / 'raw/ports' / (name + '.json'), json.loads(port.to_json()))
        port_rows.append(dict(case=name, bounds_um=bounds, legacy_status=legacy_status, v1_threshold_status=v1.value, branch_set_status=decoded.branch_set.status.value, interpretation='CLOSED means gap below specified120um threshold, not physical contact or clinical acceptance; interval-family branch cardinality UNKNOWN', roundtrip_identity_error=0, SI_conversion_error=0, sigma=None, physical_status='UNKNOWN', motion_normal_offset_m_s=req.normal_offset_m_s))
        corner_rows.append(dict(case=name, control_error_um=err, old_bounds_error_um=olderr, license='h decreasing in every b_i; g_j increasing in b_j, decreasing in every b_i (i!=j), for fixed positive a_i', rigorous_enclosure='EXACT_RATIONAL_BOX_WITH_REPRESENTED_PROJECTIONS', upstream_projection_rounding_enclosure='MISSING'))
    ex = next((x for x in port_rows if x['case'] == 'same_rms_expansion'))
    co = next((x for x in port_rows if x['case'] == 'same_rms_contraction'))
    suff = dict(summary='RMS of equal-magnitude signed axial-wall error', RMS_um=[r, r], identity_error_um=0, downstream_margin_difference_um=co['bounds_um']['gap_upper_um'][0] - ex['bounds_um']['gap_upper_um'][0], minimal_extension='signed error at axial constraint, tied to common datum and nominal spacer')
    corner_min = min((abs(v) for v in [Q(-1), Q(1)]))
    interior_min = abs(Q(0))
    false_claim = corner_min
    false_claim_rejected = not interior_min >= false_claim
    assert false_claim_rejected
    counter = dict(quantity='Euclidean clearance |z|, z∈[-1,1]', corner_min_mm=corner_min, interior_min_mm=interior_min, identical_endpoint_summary_error_mm=0, downstream_penetration_test_difference_mm=corner_min - interior_min, status='SAMPLING_ONLY', injected_false_corner_guarantee_rejected=false_claim_rejected)
    x49 = read(BASE / 'LANE_X49_DESIGN_GATE/rounds/R3.json')
    run = Path(x49['run_directory'])
    x49rows = []
    sweeps = []
    for name in ['valid', 'near_limit', 'axial_valid', 'axial_near', 'axial_penetration']:
        report = read(run / (name + '.json'))
        cp = Path(report['inputs']['crown']['path']).parent / 'contract.json'
        meta = read(cp)
        rules = []
        for key in ['material_wall', 'film_min', 'film_max', 'occlusal_contact']:
            rule = report['rules'][key]
            licensed = 'Hausdorff set triangle inequality, not corner test' if 'supplied_errors_mm' in rule else 'separately monotone barycentric signed gap in fixed-XY upper/lower vertex heights' if 'axial_error_mm' in rule else 'NONE; original UNKNOWN retained'
            rules.append(dict(quantity=key, status_before=rule['status'], status_after=rule['status'], license=licensed, proof_level='ANALYTIC_PERTURBATION_BOUND; nominal FP enclosure MISSING', source_of_box=meta.get('axial_error_source') if 'axial_error_mm' in rule else meta.get('error_source')))
        x49rows.append(dict(case=name, rules=rules, physical_calibrated_sigma=None))
        if not name.startswith('axial'):
            continue
        U = stl(report['inputs']['antagonist']['path'])[meta['regions']['antagonist']]
        L = stl(report['inputs']['crown']['path'])[meta['regions']['occlusal']]
        e = Q(str(meta['axial_error_mm']['crown'])) + Q(str(meta['axial_error_mm']['antagonist']))
        pb = make_family([[1]], [1], [0], [-1, 1], step=1, source_id='X49:' + name, source_sha256=sha(cp), unit_system='EXTERNAL_UNSPECIFIED')
        fam = certify_family(pb)
        geometry = []
        for (i, point) in enumerate(np.unique(U.reshape(-1, 3), axis=0)):
            for (j, tri) in enumerate(L):
                pts = [[Q(float(v)) for v in p] for p in [point, *tri]]
                base = copy.deepcopy([pts, pts])
                slope = [[[0, 0, e], [0, 0, 0], [0, 0, 0], [0, 0, 0]]] * 2
                geometry.append(dict(id=f'uppervertex{i}:lowerface{j}', kind='PT', skin=0, cells=[dict(base=base, slope=slope) for _ in fam['cells']]))
        ans = certify_sweep(pb, fam, geometry, complete=True, directions=[(0, 0, 1)])
        passed = verify_sweep(pb, fam, geometry, ans, complete=True) if ans['status'].startswith('SAFE') else False
        injection = None
        if passed:
            bad = copy.deepcopy(ans)
            bad['certificates'][0]['margin'] += 1
            injection = not verify_sweep(pb, fam, geometry, bad, complete=True)
            assert injection
        sweeps.append(dict(case=name, status=ans['status'], verified=passed, mutated_margin_rejected=injection, declared_pairs=len(geometry), scope='Fixed axial positive separation of all upper vertices vs all lower facets; conservative scalar row-error relaxation ±e; no unsigned norm or full insertion path claim'))
        write(ROOT / 'raw/sweeps' / (name + '.json'), dict(problem=pb, family=fam, geometry=geometry, answer=ans))
    assert len(port_rows) == 6 and len(corner_rows) == 6
    touching = Gap(Q(120, 1000), Q(120, 1000), Unit.MM, Source('threshold', sha(parent / 'PREREG_R4.json'), 'exact fixture')).classify(Q(120, 1000))
    assert touching == Status.UNKNOWN
    extra_licenses = [dict(quantity='X10R1:40+global_scanner_mean', license='affine, but transferred scalar is not a local error bound', status='DESCRIPTIVE_ONLY'), dict(quantity='X10R2:D+2000r sin(theta/2)', license='monotone D,r,theta for r>=0 and theta in[0,180deg]', status='MEDIAN_COMPOSITION_NOT_PROBABILITY_BOUND'), dict(quantity='X10R2:max(0,(80-D)/(2000sin(theta/2)))', license='nonincreasing D and theta for theta in(0,180deg]', status='DESCRIPTIVE_THRESHOLD_NOT_CALIBRATED_INTERVAL'), dict(quantity='X49:mesh/scale/insertion/milling/export', license='discrete checks, exact recession/containment or convex bur certificates; no corner-box shortcut', status='original scopes retained; normals/pose/enclosure/calibration debts remain'), dict(quantity='X49:arbitrary unsigned distance under vertex offsets', license='NONE for corner extrema; Hausdorff inequality is a separate whole-set guarantee', status='SAMPLING_ONLY for original perturbation checks; analytic envelope retains conditional scope')]
    out = dict(claim_type='capability', resolution='PER_SURFACE_REGION', corners_X10=corner_rows, corners_X49=x49rows, other_quantity_licenses=extra_licenses, family_sweep=sweeps, counterexample=counter, ports=port_rows, sufficiency=suff, injections=mutations, threshold_equality_difference=dict(X10='REJECT_CONDITIONAL_BOX', v1='UNKNOWN'), changed_conclusion=False, changed_certification=True, box_provenance=dict(X10_R1_R2='Published global means/SD or unsigned position/angular medians: descriptive only; no signed deterministic box and no calibrated sigma', X10_R3_R4='Declared signed numerical cases motivated by published regional bins, not measured full-contact bounds; published table counts are independent metrology, not seating facit', X49='Synthetic supplied surface/axial bounds in contract.json; no scanner calibration'), calibrated_sigma_requirements='same object and region, physical datum, independent reference with uncertainty, repeated signed displacement vector fields, joint spatial covariance including registration and manufacturing, heldout coverage/tails and observation law; sigma_gap^2=n^T(Cu+Cl-Cul-Clu)n; sigma is not a deterministic radius', v1_gaps=['No um enum; exact external conversion supplied', 'No built-in dental shared-lift seating solver or datum/process/calibration validator; explicit replayer supplied', 'No exact unit normal from approximate sin/cos: normal omitted, projection approximation disclosed', 'Piecewise max-affine seating cannot be put in one global AffineGap law; opaque bound evidence preserves it', 'Completeness and physical status remain distinct; no physical certification'], external_referent=dict(kind='published_code', locator='Field1216032 contact_port_v1; Motion6e2defd adapter; Motion95ec803 family_sweep; PMC10265779 tables retained', compared_quantity='exact wire/bounds and fixed signed-direction path certificates; published datum-sensitive deviations', refutes_us=True), physical_validation='UNKNOWN; own seating and X49 fixtures not external empirical facit', cost_wall_s=time.perf_counter() - t)
    write(ROOT / 'raw/CORNERS_PORT.json', out)
    state('CORNERS_PORT_DONE', 'Analytic licenses and exact codec/sweep replay recorded', 'Local force budgets and inclined frictional X18 FE')
    print('CORNERS_PORT', [(r['case'], r['status']) for r in sweeps])
if __name__ == '__main__':
    main()
