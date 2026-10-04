"""Same-patient adapter; all parent operations remain read-only."""
from dental_release.paths import expand as _release_expand
import ast, csv, importlib.util, json, math, os, sys
from pathlib import Path
import numpy as np
import trimesh
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
R = Path(__file__).resolve().parents[1]
D = R.parent.parent
RESULTS = D / 'results'
P = RESULTS / 'DEMO48_PACKAGE'
V4 = RESULTS / _release_expand('GENCAD_V4')
X69 = RESULTS / _release_expand('X69')
sys.dont_write_bytecode = True
sys.path[:0] = [str(V4 / 'payload/participant_code'), str(V4 / 'payload/participant_code/legacy/vendor'), str(X69 / 'code')]

def module(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sys.modules[name] = m
    sp.loader.exec_module(m)
    return m
x70 = module('packaged_x70_worker', RESULTS / 'LANE_X70_END_TO_END/worker.py')
write = x70.write
read = x70.read
sha = x70.sha
npz = x70.npz

def pure_functions(path, names, globals_):
    """Execute unchanged pure package functions without unrelated volume-I/O imports."""
    tree = ast.parse(Path(path).read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    if {n.name for n in nodes} != set(names):
        raise ValueError('Missing package pure function')
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), globals_)
    return [globals_[n] for n in names]
(stl, transform) = pure_functions(X69 / 'code/registration.py', ['stl', 'transform'], {'np': np, 'pathlib': __import__('pathlib')})
from legacy.geometry import grid, height, crop, shell, area_weights
SEL = read(R / 'PATIENT_SELECTION.json')
PR = read(R / 'PREREG_R1_SAME_PATIENT.json')
BASE = Path(SEL['source_base'])
CASE = SEL['patient_id']

def task(out):
    t = read(out / 'TASK.json')
    a = np.load(out / 'TASK_ARRAYS.npz')
    for k in ['xy', 'faces', 'uv', 'prior', 'ceiling', 'weights', 'preparation_z']:
        t[k] = a[k]
    t['A'] = csr_matrix((a['a_data'], a['a_indices'], a['a_indptr']), shape=(int(a['a_rows']), len(a['xy'])))
    t['obstacle_b'] = a['b']
    return t

def source(out):
    w = np.load(SEL['witness_npz'])
    T = w['lower_T']
    U = w['upper_T']
    tri = transform(stl(BASE / 'IOS/teeth_seg/L/36.stl').reshape(-1, 3), T).reshape(-1, 3, 3)
    base = (tri.reshape(-1, 3).min(0) + tri.reshape(-1, 3).max(0)) / 2
    local = tri - base
    half = np.ptp(local.reshape(-1, 3)[:, :2], axis=0) / 2 * PR['construction']['footprint_fraction']
    (xy, f, uv) = grid(half, 'molar_crown', PR['construction']['grid_n'])
    ref = height(local, xy)
    missing = ~np.isfinite(ref)
    if missing.all():
        raise ValueError('No measured roof support at frozen footprint')
    if missing.any():
        (_, ii) = cKDTree(xy[~missing]).query(xy[missing])
        ref[missing] = ref[~missing][ii]
    prep = float(np.median(ref) - PR['construction']['preparation_below_median_roof_mm'])
    upp = []
    origins = []
    for p in sorted((BASE / 'IOS/teeth_seg/U').glob('*.stl')):
        tt = transform(stl(p).reshape(-1, 3), U).reshape(-1, 3, 3) - base
        kept = crop(tt, xy)
        upp.append(kept)
        origins.append(dict(fdi=int(p.stem), sha256=sha(p), total_facets=len(tt), retained_projected_facets=len(kept)))
    upper = np.concatenate(upp)
    if not len(upper):
        raise ValueError('No same-patient antagonist triangles over selected roof')
    from obstacle import constraints
    (A, b) = constraints(upper, xy, f, PR['construction']['requirements']['clearance_mm'])
    ceiling = height(upper, xy, upper=True)
    (v, ff) = shell(xy, ref, np.full(len(ref), prep - 2), f)
    roof = trimesh.Trimesh(v, ff, process=False)
    if roof.volume < 0:
        roof.invert()
    if not roof.is_watertight:
        raise ValueError('Declared roof cap is not watertight')
    roof.export(out / 'measured_capped_roof.stl')
    triangles = upper.astype(np.float32).astype(float)
    antagonist = trimesh.Trimesh(triangles.reshape(-1, 3), np.arange(len(triangles) * 3).reshape(-1, 3), process=False)
    antagonist.export(out / 'antagonist.stl')
    trimesh.Trimesh(local.reshape(-1, 3), np.arange(len(local) * 3).reshape(-1, 3), process=False).export(out / 'segmented_tooth.stl')
    t = dict(task_id='X77_Demo1_FDI36_R1', case_key=CASE, family='molar_crown', level='frozen_scenario', units='mm', frame='Demo1_local_translated_CBCT', status='READY', geometry_file='TASK_ARRAYS.npz', requirements=PR['construction']['requirements'], preparation_height_mm=prep, material=PR['construction']['material'], connector_x=[-half[0] / 3, half[0] / 3], implant_axis=[0, 0, 1], minimum_strut_mm=1.0, scope='Restricted roof; complete axial crown UNKNOWN; independent metric scale/bite UNKNOWN')
    write(out / 'TASK.json', t)
    npz(out / 'TASK_ARRAYS.npz', xy=xy, faces=f, uv=uv, prior=ref, ceiling=ceiling, weights=area_weights(xy, f), preparation_z=np.full(len(xy), prep), a_data=A.data, a_indices=A.indices, a_indptr=A.indptr, a_rows=np.array(A.shape[0]), b=b)
    npz(out / 'scene.npz', xy=xy, faces=f, reference_z_mm=ref, base_mm=base, local_to_cbct_rotation=np.eye(3), lower_IOS_to_CBCT=T, upper_IOS_to_CBCT=U, source_tooth_triangles_mm=local)
    write(out / 'SCAN.json', dict(patient_id=CASE, fdi=36, frame=t['frame'], base_CBCT_mm=base, local_to_CBCT_rotation=np.eye(3), source_lower_IOS_to_CBCT=T, source_upper_IOS_to_CBCT=U, units=SEL['units'], resolution='PER_POINT', time_scale='SIMULTANEOUS', source_tooth_sha256=sha(BASE / 'IOS/teeth_seg/L/36.stl'), source_triangles=len(tri), roof_nodes=len(xy), roof_missing_nodes=int(missing.sum()), roof_missing_fraction=float(missing.mean()), completion='Nearest observed projected roof only at missing nodes; virtual flat basal cap', antagonist_sources=origins, antagonist_facets=len(upper), projected_constraint_rows=A.shape[0], source_status='EXTERNALLY_MEASURED_COORDINATES_WITH_UPSTREAM_SEGMENTATION_CLOSURE', interjaw_acquisition_bite='UNKNOWN', anatomical_regions='Publisher FDI36; no clinical finish line', preparation='VIRTUAL_PLANE_NO_DENTIN_OR_PULP_LABELS', external_referent=PR['external_referent']))
    xr = read(X69 / 'raw/R2_RESULTS.json')
    rows = xr['regions']
    with (out / 'REGISTRATION_PER_REGION.csv').open('w') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    checks = []
    for row in rows:
        j = row['jaw']
        mask = (w[j + '_fdi'] == row['fdi']) & (w[j + '_region'] == row['spatial_z_half'])
        ds = w[j + '_residual_mm'][mask]
        measured = float(np.percentile(ds, 95))
        checks.append(dict(jaw=j, fdi=row['fdi'], half=row['spatial_z_half'], n=len(ds), p95_mm=measured, source_p95_mm=row['p95_mm'], agreement_error_mm=abs(measured - row['p95_mm']), pass_gate=measured <= 0.5, heldout=row['heldout_from_fit']))
    cb = stl(BASE / 'CBCT/reconstruction/out_smoothed.stl')
    cb = cb[cb.mean(1)[:, 2] < 39]
    m = trimesh.Trimesh(cb.reshape(-1, 3), np.arange(len(cb) * 3).reshape(-1, 3), process=False)
    ids = np.flatnonzero(w['lower_fdi'] == 36)[::38]
    q = w['lower_transformed_points'][ids]
    d = trimesh.proximity.closest_point(m, q)[1]
    d0 = trimesh.proximity.closest_point_naive(m, q)[1]
    error = float(np.max(abs(d - d0)))
    write(out / 'REGISTRATION.json', dict(patient_id=CASE, source_result_sha256=sha(X69 / 'raw/R2_RESULTS.json'), immutable_witness_sha256=sha(SEL['witness_npz']), regions=checks, regional_max_statistic_discrepancy_mm=max((x['agreement_error_mm'] for x in checks)), regions_failed=sum((not x['pass_gate'] for x in checks)), regions_count=len(checks), heldout_regions_pass=sum((x['pass_gate'] for x in checks if x['heldout'])), heldout_regions_count=sum((x['heldout'] for x in checks)), selected_tooth_regions=[x for x in checks if x['fdi'] == 36], metric='Nearest CBCT-only reconstructed tooth surface consistency; NOT independent TRE', independent_TRE_mm=None, uniform_registration_error_mm=None, scope='PER_SURFACE_REGION samples; sampling maximum is not a uniform anatomical bound', published_triangle_control_error_mm=error, published_triangle_control_pass=error <= 1e-09, injected_distance_plus_1mm_rejected=float(np.max(abs(d + 1 - d0))) > 1e-09, external_referent=PR['external_referent']))

def packaged(stage, out):
    x70.task = lambda : task(out)
    x70.R = R
    getattr(x70, stage)(out)
    if stage == 'generate':
        p = read(out / 'GENERATION.json')
        p['patient_case'] = CASE
        p['patient_id'] = CASE
        p['fdi'] = 36
        p['metadata_adapter'] = 'X70 original fixed case27 label replaced by actual TASK identity; no numerical change'
        p['physical_force'] = 'UNKNOWN_UNMATCHED_PATIENT_SPECIMEN'
        write(out / 'GENERATION.json', p)
        from legacy.checks import all_checks
        from legacy.generators import submit
        t = task(out)
        f = np.load(out / 'DESIGN_FIELDS.npz')
        z = f['outer_z_mm']
        inner = f['inner_z_mm']
        delta = 1.0 + max(0.0, float(np.max(t['obstacle_b'] - t['A'] @ z)))
        injected = all_checks(t, submit(t, z + delta, inner))
        direct = t['A'] @ (z + delta) - t['obstacle_b']
        write(out / 'GENERATOR_FAULT_REPAIR.json', dict(round='R1B', inherited_plus_1mm_instrument=read(out / 'GENERATOR_FAULTS.json'), reason='A +1mm move does not collide with this patients 5mm gap; original failure is preserved.', new_shift_mm=delta, minimum_direct_collision_mm=float(direct.min()), injected_check=injected, actual_collision_rejected=injected['validity'] == 'FAIL' and injected['checks'].get('antagonist', {}).get('status') == 'FAIL', unchanged_tolerance_mm=1e-06))
    if stage == 'gate':
        write(out / 'EXPORT_FRAME.json', dict(patient_id=CASE, frame=task(out)['frame'], local_to_CBCT=dict(rotation=np.eye(3), translation_mm=np.load(out / 'scene.npz')['base_mm']), source_scale='conditional_mm', generation_hash=sha(out / 'crown.stl'), physical_release='NOT_SUPPORTED', rigid_same_jaw_design_checks='Euclidean thickness/film/milling within the same IOS frame invariant to common rigid pose; cross-jaw contact depends on both maps and unknown bite state'))

def physical(out):
    tree = ast.parse((V4 / 'code/literature.py').read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'calibrated_force' or (isinstance(n, ast.Assign) and any((isinstance(a, ast.Name) and a.id == 'PROTOCOL_FIELDS' for a in n.targets)))]
    ns = {'math': math}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'packaged_calibrated_force', 'exec'), ns)
    t = task(out)
    f = np.load(out / 'DESIGN_FIELDS.npz')
    query = dict(thickness_mm=float(np.median(f['local_thickness_mm'])), product=t['material'], specimen_geometry='X77_Demo1_FDI36_restricted_roof', support='X60 rigid basal nodes', cement='UNKNOWN_SEATED_FILM', angle_deg=0, indenter_diameter_mm=None, crosshead_mm_min=None, contact_layer='X60 normalized unit nodal load', aging=None)
    answer = ns['calibrated_force']([], query)
    lit = read(V4 / 'FROZEN_LITERATURE_PREDICTIONS.json')['rows']
    rows = []
    for q in lit:
        force = q['kind'] == 'crown'
        error = abs(math.log(q['observed'] / q['predicted'])) if force else abs(q['observed'] - q['predicted'])
        tol = 0.2 if force else 20.0
        wrong = q['observed'] * math.exp(1) if force else q['observed'] + 100
        bad = abs(math.log(wrong / q['predicted'])) if force else abs(wrong - q['predicted'])
        rows.append({**q, **dict(error=error, tolerance=tol, gate=error <= tol, injected_observation=wrong, injected_error=bad, injected_wrong_value_rejected=bad > tol, resolution='POPULATION')})
    curated = read(V4 / 'payload/literature/CROWN_CURATED.json')
    groups = {}
    for q in curated:
        if q['eligible'] and q['protocol'].startswith('crown'):
            groups.setdefault((q['study'], q['protocol'], q['product']), []).append(q)
    ctrl = []
    for (key, gg) in groups.items():
        gg = sorted(gg, key=lambda q: q['t'])
        if len(set((q['t'] for q in gg))) < 3:
            continue
        (a, b) = (gg[0], gg[-1])
        for q in gg[1:-1]:
            alpha = math.log(q['t'] / a['t']) / math.log(b['t'] / a['t'])
            pred = math.exp((1 - alpha) * math.log(a['mean']) + alpha * math.log(b['mean']))
            src = next((x for x in lit if x['kind'] == 'crown' and x['study'] == q['study'] and (x['protocol'] == q['protocol']) and (x['x'] == q['t'])))
            ctrl.append(dict(locator='doi:' + q['doi'] + ' ' + q['locator'], predicted_N=pred, observed_N=q['mean'], source_prediction_N=src['predicted'], difference_N=abs(pred - src['predicted'])))
    force = [q for q in rows if q['kind'] == 'crown']
    gap = [q for q in rows if q['kind'] == 'cement']
    write(out / 'PHYSICAL.json', dict(patient_id=CASE, patient_force_query=query, calibrated_patient_force=answer, calibrated_force05_N=None, seated_spatial_cement_film_um=None, nominal_film_um=(f['inner_z_mm'] - f['virtual_prep_z_mm']) * 1000, nominal_film_resolution='PER_POINT', source_controls=rows, independent_endpoint_control=ctrl, source_force_pass_count=sum((q['gate'] for q in force)), source_force_count=len(force), source_gap_pass_count=sum((q['gate'] for q in gap)), source_gap_count=len(gap), source_support_gate=False, force_likelihood_transfer='UNKNOWN_NO_MATCHED_PROTOCOL', same_jaw_registration_influence='Exactly zero for Euclidean CAD self-checks under a shared rigid transform; physical support/load/material still UNKNOWN', cross_jaw_force_influence='UNKNOWN_ACTUAL_CONTACT_AND_BITE_STATE', affine_design_sensitivity_enclosure='MISSING; no robust mechanics claim', external_referent=dict(kind='independent_measurement', locator='doi:10.4047/jap.2021.13.5.269 Table1; doi:10.3390/ma17020365 Table2', compared_quantity='Published protocol-specific crown group fracture force; not generated specimen', refutes_us=True)))

def rational_decisions(out):
    from fractions import Fraction as F
    from scipy.linalg import solve
    (rational_inverse, apply_rational, greedy) = pure_functions(X69 / 'code/round3.py', ['rational_inverse', 'apply_rational', 'greedy'], {'F': F, 'np': np})
    w = np.load(SEL['witness_npz'])
    (anchors, ids) = greedy(w['lower_source_points'])
    M = np.vstack((anchors.T, np.ones(4)))
    inv = rational_inverse(M)
    s = np.load(out / 'scene.npz')
    base = s['base_mm']
    T = s['lower_IOS_to_CBCT']
    it = np.linalg.inv(T)
    f = np.load(out / 'DESIGN_FIELDS.npz')
    pts = np.c_[s['xy'], f['virtual_prep_z_mm']] + base
    cb = stl(BASE / 'CBCT/reconstruction/out_smoothed.stl')
    root = cb[:, 0][cb[:, 0, 2] < 30]
    root = root[np.linspace(0, len(root) - 1, 12, dtype=int)]
    queries = np.concatenate([pts, root])
    q = transform(queries, it)
    records = []
    weights = []
    maxctrl = 0.0
    exact = True
    for (i, x) in enumerate(q):
        h = np.r_[x, 1.0]
        lam = apply_rational(inv, h)
        weights.append(lam)
        identity = all((sum((lam[k] * F(float(M[j, k])) for k in range(4))) == F(float(h[j])) for j in range(4)))
        exact &= identity
        L = sum((abs(z) for z in lam))
        err = float(np.max(abs(solve(M, h) - np.array([float(z) for z in lam]))))
        maxctrl = max(maxctrl, err)
        extreme = sum((z * F(1 if z >= 0 else -1) for z in lam))
        records.append(dict(query_id=i, role='preparation_to_pulp_input' if i < len(pts) else 'CBCT_tooth_root_witness_not_canal', CBCT_point_mm=queries[i], IOS_point_mm=x, lambda_exact=[str(z) for z in lam], absolute_error_amplification_exact=str(L), absolute_error_amplification_upper=math.nextafter(float(L), math.inf), resolution='PER_POINT', exact_reconstruction=identity, independent_scipy_error=err, actual_affine_extremal_displacement_per_unit_epsilon=str(extreme), extremal_anchor_sign_fault_half_bound_rejected=abs(extreme) > L / 2, epsilon_anchor_mm=None, physical_distance_mm=None))
    pairs = []
    for i in range(0, len(pts) - 1, 5):
        j = i + 1
        dx = [a - b for (a, b) in zip(weights[i], weights[j])]
        relative = sum((abs(z) for z in dx))
        naive = sum((abs(z) for z in weights[i])) + sum((abs(z) for z in weights[j]))
        pairs.append(dict(i=i, j=j, relative_coeff_exact=str(relative), relative_coeff_upper=math.nextafter(float(relative), math.inf), independent_absolute_sum=float(naive), gain_in_conditional_coefficient=float(naive / relative) if relative else None, signed_difference_weights=[str(z) for z in dx], same_pose_cancellation_exact=sum(dx) == 0, physical_calibration='UNKNOWN; these are conditional coefficients, not measured errors'))
    dx = F(1, 4)
    xp = F(0)
    tissue = F(2)
    distance_plus = abs(xp + dx - tissue)
    distance_minus = abs(xp - dx - tissue)
    sufficient = dict(summary='Each of four anchor error magnitudes, RMS and p95=.25mm', state_A_translation=[0.25, 0, 0], state_B_translation=[-0.25, 0, 0], summary_identity_error=0.0, same_summary_machine_exact=True, downstream_fixed_tissue_distance_A_mm=float(distance_plus), downstream_fixed_tissue_distance_B_mm=float(distance_minus), downstream_difference_mm=float(abs(distance_plus - distance_minus)), same_pose_relative_distance_A_mm=2.0, same_pose_relative_distance_B_mm=2.0, minimal_extension='For the tested query keep signed projection of the shared pose onto the tissue direction; for arbitrary queries keep directed shared map/error state. Magnitude summaries alone cannot decide.', external_referent=dict(kind='our_own_fixture', locator='PREREG_R2_CORRELATED_DECISIONS.json; exact dyadic rigid-translation construction', compared_quantity='Summary sufficiency counterexample, not patient anatomical truth', refutes_us=True))
    support = np.array([[-1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0]])
    aa = np.eye(4)
    bb = np.eye(4)
    bb[:3, :3] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    da = cKDTree(support).query(transform(support, aa))[0]
    db = cKDTree(support).query(transform(support, bb))[0]
    rp = np.array([[1.0, 0.0, -10.0]])
    nearest_summary = dict(identity_error_mm=float(np.max(abs(da - db))), machine_exact=bool(np.array_equal(da, db)), same_RMS=True, same_p95=True, downstream_root_displacement_mm=float(np.linalg.norm(transform(rp, aa) - transform(rp, bb))), minimal_extension='Directed asymmetric independent correspondences plus jaw state', kind='our_own_fixture', locator=str(X69 / 'code/contracts.py'))
    from guards import registration_bound, admit_tissue, rejects
    bound_bad = dict(kind='p95_surface_residual', epsilon_mm=0.3, independently_measured=True)
    guard_controls = dict(source_residual_as_uniform_bound_rejected=rejects(lambda : registration_bound(bound_bad)), missing_pulp_geometry_rejected=rejects(lambda : admit_tissue(None, 'pulp', CASE)), foreign_patient_mask_rejected=rejects(lambda : admit_tissue(dict(patient_id='ToothFairy_P001', tissue='pulp'), 'pulp', CASE)), tooth_exterior_as_pulp_rejected=rejects(lambda : admit_tissue(dict(patient_id=CASE, tissue='tooth_exterior'), 'pulp', CASE)), unlocated_bone_rejected=rejects(lambda : admit_tissue(dict(patient_id=CASE, tissue='bone'), 'bone', CASE)))
    faults = dict(exact_weights=exact, independent_control_pass=maxctrl <= 1e-08, half_bound_faults_rejected=all((r['extremal_anchor_sign_fault_half_bound_rejected'] for r in records)), scalar_summary_injection_rejected=sufficient['downstream_difference_mm'] >= 0.01, **guard_controls)
    write(out / 'UNCERTAINTY.json', dict(patient_id=CASE, anchor_points_IOS_mm=anchors, anchor_point_indices=ids, records=records, shared_pairs=pairs, max_scipy_weight_error=maxctrl, proof='For affine D and partition-of-unity lambda, D(x)=sum lambda_i D(a_i), D(x)-D(y)=sum (lambda_i(x)-lambda_i(y))D(a_i). Triangle inequality; distance-to-fixed-set is 1-Lipschitz. All weights exactly rational on stored coordinates. Coefficients outward rounded. epsilon and tissue Hausdorff error UNKNOWN.', enclosure_status='RIGOROUS_CONDITIONAL_AFFINE_ENCLOSURE_NOT_MEASURED_ERROR', registration_epsilon_mm=None, tissue_error_mm=None, sufficiency=sufficient, nearest_surface_summary_counterexample=nearest_summary, controls=faults, external_referent=read(R / 'PREREG_R2_CORRELATED_DECISIONS.json')['external_closed_form']))
    event = module('packaged_postcut_event', RESULTS / 'LANE_X31_CLINICAL_ANSWERS/minimum_observation.py')
    absent = event.event_interval(0, 0, 1.5, 1.5, False)
    decisions = []
    for (name, tissue) in [('pulp_margin', 'pulp'), ('residual_dentin', 'DEJ_and_pulp'), ('implant_bone_width_height', 'bone'), ('implant_canal_clearance', 'canal')]:
        decisions.append(dict(patient_id=CASE, quantity=name, status='UNKNOWN', tissue=tissue, geometry_present=False, nominal_value_mm=None, registration_contribution_mm=None, segmentation_contribution_mm=None, decision=None, resolution='PER_POINT port; PER_TOOTH or PER_SURFACE_REGION aggregation blocked', time_scale='HANDOVER' if 'pulp' in name or 'dentin' in name else 'SIMULTANEOUS', reason='No same-patient ' + tissue + ' geometry; no independent uniform directed registration error', minimum_observation='Same patient calibrated CBCT volume/semantic tissue, independent region error bounds; independent directed registration anchors. Residual dentin also needs DEJ or same-ray postcut material event.', rigorous_formula='[max(0,d_nom - epsilon*L - eta_tissue - eta_prep), d_nom + epsilon*L + eta_tissue + eta_prep] for measured finite d_nom; no values inserted'))
    write(out / 'ANATOMICAL_DECISIONS.json', dict(patient_id=CASE, decisions=decisions, packaged_X31_absent_event=absent, borrowed_masks_admitted=0, guard_controls=guard_controls, physical_complete_flow=False, candidate_implant_site='UNKNOWN_NO_ANNOTATED_ALVEOLAR_REGION; published crown gaps alone do not define a surgical site', external_referent=PR['external_referent']))
    ports = [dict(consumer='SDF / virtual preparation location', registration_effect='absolute pose error when placed in CBCT; local shape invariant', rigorous_term='epsilon_lower * L_prepared_point', numeric_mm=None, resolution='PER_POINT', time_scale='HANDOVER'), dict(consumer='CAD wall / nominal film / ideal insertion', registration_effect='0 under shared proper rigid map of crown + die', rigorous_term='0 common-pose contribution; source geometry and manufacture errors separate', numeric_common_rigid_mm=0.0, resolution='PER_POINT', time_scale='SIMULTANEOUS'), dict(consumer='opposing-jaw contact', registration_effect='independent jaw maps and bite acquisition state', rigorous_term='epsilon_lower*L_lower + epsilon_upper*L_upper + eta_lower_surface + eta_upper_surface + unknown_bite_state', numeric_mm=None, resolution='PER_POINT', time_scale='SIMULTANEOUS'), dict(consumer='modeled fracture force ratio', registration_effect='0 if entire modeled load/support/material scene rotates together; actual contact force UNKNOWN', rigorous_term='Rigid rotation invariance of isotropic model, NOT a physical force uncertainty enclosure', numeric_common_rigid_mm=0.0, resolution='PER_POINT -> PER_TOOTH conditional model', time_scale='HANDOVER'), dict(consumer='milling tool envelope', registration_effect='0 if tool and crown poses share the transformed scene; machine fixture frame error remains UNKNOWN', rigorous_term='0 common pose; no calibrated machine/fixture envelope', numeric_common_rigid_mm=0.0, resolution='PER_POINT', time_scale='HANDOVER'), dict(consumer='export / metrology datum', registration_effect='local mesh exact transport; CBCT placement carries directed lower pose; actual lab datum UNKNOWN', rigorous_term='local transport error + epsilon_lower*L_point when referencing CBCT', numeric_mm=None, resolution='PER_POINT', time_scale='HANDOVER'), dict(consumer='preparation -> CBCT pulp margin', registration_effect='IOS-defined preparation vs fixed CBCT pulp', rigorous_term='epsilon_lower*L_prepared_point + eta_pulp + eta_preparation; pulp missing', numeric_mm=None, resolution='PER_POINT', time_scale='HANDOVER'), dict(consumer='remaining dentin', registration_effect='same local spatial enclosure plus DEJ or observed material event', rigorous_term='pulp distance interval does not identify dentin without DEJ/postcut event', numeric_mm=None, resolution='PER_POINT', time_scale='HANDOVER'), dict(consumer='bone width/height in CBCT-native anatomy', registration_effect='0 IOS registration effect for a site/axis defined wholly in CBCT; tissue and site currently missing', rigorous_term='0 if native; epsilon_lower*L for IOS-defined site/orientation. Origin/axis contract required.', numeric_common_IOS_mm=0.0, resolution='PER_POINT / PER_SURFACE_REGION', time_scale='SIMULTANEOUS'), dict(consumer='implant -> canal clearance', registration_effect='0 common IOS error if implant and canal defined CBCT-native; IOS-planned tool against fixed CBCT canal needs mapped whole-tool points', rigorous_term='sup_tool_points epsilon_lower*L_tool_point + eta_canal + eta_tool; no canal or qualified site', numeric_mm=None, resolution='PER_POINT before region minimum', time_scale='SIMULTANEOUS')]
    write(out / 'DECISION_ERROR_MAP.json', dict(patient_id=CASE, epsilon_lower_mm=None, epsilon_upper_mm=None, ports=ports, scope='Each term is rigorously conditional in affine point geometry; full physical, mechanical and scanner uncertainty not calibrated. No surface p95 used as epsilon.'))

def main():
    stage = sys.argv[1]
    out = Path(sys.argv[2])
    out.mkdir(exist_ok=True, parents=True)
    if stage == 'source':
        source(out)
    elif stage == 'physical':
        physical(out)
    elif stage == 'uncertainty':
        rational_decisions(out)
    else:
        packaged(stage, out)
if __name__ == '__main__':
    main()
