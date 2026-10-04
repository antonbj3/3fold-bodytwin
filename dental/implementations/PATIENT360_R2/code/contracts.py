"""Decision-relative X26 classes and executable shared-theta prerequisites."""
import copy, time, sys
from common import *
from region_field import RegionField
LEVELS = {'PER_POINT': 0, 'PER_SURFACE_REGION': 1, 'PER_TOOTH': 2, 'PER_ARCH': 3, 'POPULATION': 4}

def validate_decision(row, registry):
    from decision_binding import validate_identity
    validate_identity(row, registry)
    for key in ['id', 'question', 'status', 'resolution', 'time_scale', 'decidability', 'binding_quantity', 'evidence', 'theta_sha256', 'physical_status', 'support']:
        if key not in row:
            raise ValueError('DECISION_MISSING:' + key)
    if row['resolution'] not in LEVELS:
        raise ValueError('RESOLUTION_REQUIRED')
    if row['time_scale'] not in ['SIMULTANEOUS', 'HANDOVER']:
        raise ValueError('TIME_SCALE_REQUIRED')
    d = row['decidability']
    if d.get('class') not in [1, 2, 3] or not d.get('refinement_variable') or (not d.get('justification')):
        raise ValueError('DECIDABILITY_QUERY_REQUIRED')
    if d['class'] == 3 and (not d.get('divergence_proof')):
        raise ValueError('NO_DIVERGENCE_PROOF')
    if d['class'] == 2 and (not row['binding_quantity']):
        raise ValueError('BINDING_MEASUREMENT_REQUIRED')
    if row['theta_sha256'] != digest(registry):
        raise ValueError('SHARED_THETA_MISMATCH')
    if not row['support'].get('patient_id') or not row['support'].get('frame_id'):
        raise ValueError('SUPPORT_IDENTITY_REQUIRED')
    if not row['evidence']:
        raise ValueError('EVIDENCE_REQUIRED')
    if row['physical_status'] == 'CERTIFIED' and (not row.get('physical_validation_locator')):
        raise ValueError('PHYSICAL_PROMOTION_WITHOUT_FACIT')
    return True

def edge(producer, consumer, producer_resolution, consumer_resolution, *, time_scale, patient_id, source_support):
    fine = min([producer_resolution, consumer_resolution], key=LEVELS.get)
    if time_scale not in ['SIMULTANEOUS', 'HANDOVER']:
        raise ValueError('EDGE_TIME_SCALE')
    return dict(producer=producer, consumer=consumer, resolution=fine, producer_resolution=producer_resolution, consumer_resolution=consumer_resolution, time_scale=time_scale, patient_id=patient_id, source_support=source_support)

def validate_edge(e):
    required = min([e['producer_resolution'], e['consumer_resolution']], key=LEVELS.get)
    if e['resolution'] != required:
        raise ValueError('EDGE_NOT_FINEST_COMMON_SUPPORT')
    if e['time_scale'] not in ['SIMULTANEOUS', 'HANDOVER'] or not e.get('source_support'):
        raise ValueError('EDGE_SUPPORT_OR_TIME')
    return True
IOS_SPECS = {'X11': ('Are the tooth/ FDI labels anatomically correct?', 'PER_POINT', 'SIMULTANEOUS', 2, 'Blind FDI/facet annotation on this source scan'), 'X21': ('What surfaces do physical carry bite power?', 'PER_POINT', 'SIMULTANEOUS', 2, 'Loaded relative pose and source-addressed regional force'), 'X7': ('Is the reporting category valid on the exact region under consideration?', 'PER_SURFACE_REGION', 'SIMULTANEOUS', 2, 'Scope-matched independent clinical report; no global conclusion from lateral absence'), 'X18/X18b': ('Is the generated surface a complete possible preparation/crown ?', 'PER_POINT', 'HANDOVER', 2, 'Actual preparation, cervical margin and intaglio'), 'X34': ('Can this design be released as crown ?', 'PER_POINT', 'HANDOVER', 2, 'Measured preparation, valid wall/film and applicable material design limits; original FAIL retained'), 'X1b': ('Does the frozen break load forecast apply to this individual?', 'PER_TOOTH', 'SIMULTANEOUS', 2, 'Same-specimen geometry/support/material; other-patient X1b cannot bind'), 'X18_FE': ('What physical top voltage does this tooth ?', 'PER_POINT', 'SIMULTANEOUS', 2, 'Force allocation, loaded pose, support, material and mesh-convergence bound'), 'X38': ('Is source identity, coordinates and contracts retained during data transport?', 'PER_POINT', 'SIMULTANEOUS', 1, 'None for byte/typed metadata transport; manufacture remains blocked')}
CBCT_SPECS = {'X12': ('Is the pulp measure actually remaining dentin after preparation ?', 'PER_POINT', 'HANDOVER', 2, 'DEJ and actual prepared surface on this tooth; total hard tissue is not dentin'), 'X8/X31': ('Which individual nerve injury risk has the guide placement?', 'PER_SURFACE_REGION', 'HANDOVER', 2, 'Same-case achieved pose and anatomical canal boundary/error law'), 'X15': ('How can tooth be moved biologically within the bone?', 'PER_SURFACE_REGION', 'HANDOVER', 2, 'CEJ/full root and longitudinal biological response; digital intervals do not identify movement'), 'X20': ('Is there a background for the implant /sine query here?', 'PER_TOOTH', 'HANDOVER', 2, 'Actual edentulous site, bone/sinus boundary and intended axis; occupied sites abstain'), 'X16': ('Does the published heat/output model apply to this individual?', 'PER_POINT', 'HANDOVER', 2, 'Same-specimen heat/cooling protocol and measured thermal response'), 'X33': ('What physical powder temperature gives the preparation?', 'PER_POINT', 'HANDOVER', 2, 'Registered tool path, heat input, coolant and pulp thermometry'), 'X26': ('Can Numeric refinement determine the physical query ?', 'PER_SURFACE_REGION', 'SIMULTANEOUS', 2, 'Question-specific physical error floor; scalar scenario is not a calibrated anatomical bound')}

def run(out, field_result, fields, contact, design):
    t0 = time.perf_counter()
    out.mkdir(parents=True, exist_ok=True)
    checks = []
    patients = []
    edges = []
    allrows = []
    shared_execution = []
    freeze(out / 'FROZEN_PREDICTIONS.json', dict(claim_type='capability', prereg=source(ROOT / 'PREREG_C3_CONTRACT.json'), expected_release='ABSTAIN for all cases; absent physical observations cannot be replaced by numerical refinement', measurement_status='Physical measurement NOT_RUN', expected_untyped_decisions=0))
    for pid in field_result['patients']:
        own = [f for f in fields if f.meta['patient_id'] == pid]
        frameid = own[0].meta['frame_id']
        ios = pid.startswith('Bite2Text')
        oldpath = PARENT / 'patients' / pid / 'DECISION_CONTRACT.json'
        if oldpath.exists():
            old = load(oldpath)
            source(oldpath)
        else:
            old = dict(decisions=[dict(tool='X11', status='REUSED_SOURCE_LABELS', physical='UNKNOWN_FDI'), dict(tool='X21', status='SOURCE_GEOMETRY_QUERIED', physical='UNKNOWN_CONTACT')], release='ABSTAIN')
        if ios:
            reg = dict(patient_id=pid, frame_id=frameid, parameters={'pose_delta_mm': dict(support=[-0.05, 0.05], unit='mm', resolution='PHENOMENOLOGICAL', replacement='Same-patient registered loaded pose'), 'force_N': dict(value=100.0, unit='N', resolution='PHENOMENOLOGICAL', replacement='Per-tooth force and per-region allocation'), 'E_MPa': dict(value=210000.0, unit='MPa', resolution='PHENOMENOLOGICAL', replacement='Same material batch stiffness'), 'support': dict(value=None, resolution='PHENOMENOLOGICAL', replacement='Actual preparation/cement/tooth or die support')}, rows=[dict(pose_delta_mm=d, force_N=100.0, E_MPa=210000.0, nu=0.3) for d in [-0.05, 0.0, 0.05]], probability_law='NONE; deterministic scenarios, not independent marginals')
            specs = IOS_SPECS
        else:
            reg = dict(patient_id=pid, frame_id=frameid, parameters={'boundary_error_mm': dict(support=None, unit='mm', resolution='PER_POINT', replacement='Independent anatomical boundaries'), 'retained_heat_fraction': dict(value=0.02, resolution='PHENOMENOLOGICAL', replacement='Same-tooth retained heat'), 'cooling_h_W_m2K': dict(support=[0.0, 2000.0], unit='W/m2/K', resolution='PHENOMENOLOGICAL', replacement='Measured thermal boundary flux/temperature')}, rows=[dict(retained_heat_fraction=0.02, cooling_h_W_m2K=h) for h in [0.0, 500.0, 2000.0]], probability_law='NONE')
            specs = CBCT_SPECS
        if pid == 'Bite2Text_F4775':
            upper = next((f for f in own if f.meta['jaw'] == 'upper'))
            lower = next((f for f in own if f.meta['jaw'] == 'lower'))
            roof = np.load(PARENT_DATA / pid / 'roof.npz')
            arch = np.load(PARENT_DATA / pid / 'contact_map.npz')
            xyz = np.c_[roof['xy'], roof['z']]
            ai = np.linspace(0, len(arch['xy']) - 1, 129, dtype=int)
            archxyz = np.c_[arch['xy'][ai], arch['lower_z'][ai]]
            sys.path.insert(0, str(RESULTS / 'LANE_X18_CROWN_ANTAGONIST/code'))
            from experiment import contact as contact_op
            from fe import normalize, solve_cases
            loads = {}
            for theta in reg['rows']:
                query = upper.query(xyz, patient_id=pid, frame_id=frameid, theta=theta)
                qa = upper.query(archxyz, patient_id=pid, frame_id=frameid, theta=theta)
                ql = lower.query(archxyz, patient_id=pid, frame_id=frameid, theta=theta)
                co = contact_op(query['vertical_gap_mm'], roof['z'], roof['xy'], roof['index'], roof['weights'], 0.1, 0.1)
                tid = digest(theta)
                loads[tid] = normalize(co['mask'], roof['weights'], theta['force_N'])
                shared_execution.append(dict(patient_id=pid, theta_id=tid, theta=theta, source_theta_ids=[query['theta_id'], qa['theta_id'], ql['theta_id']], crown_contact_area_mm2=co['area_mm2'], crown_active_point_ids=np.flatnonzero(co['mask']), arch_sample_point_ids=ai, arch_sample_gap_mm=qa['height_mm'] - ql['height_mm'], resolution='PER_POINT', time_scale='SIMULTANEOUS', scope='Declared query subset, no whole-arch area claim'))
            (answers, cost) = solve_cases(roof['xy'], roof['z'], roof['faces'], roof['uv'], loads, dict(thickness_mm=1.2, E_MPa=210000.0, nu=0.3))
            for rr in shared_execution:
                rr['FE'] = answers[rr['theta_id']]
                checks.append(check('shared_theta:' + rr['theta_id'], all((t == rr['theta_id'] for t in rr['source_theta_ids'])), digest({**rr['theta'], 'pose_delta_mm': 99.0}) != rr['theta_id']))
            dump(out / 'SHARED_THETA_EXECUTION.json', dict(rows=shared_execution, FE_cost=cost, physical_status='UNKNOWN', force_allocation='PHENOMENOLOGICAL equal pressure on inherited geometric patch'))
        r = []
        for prev in old['decisions']:
            tool = prev['tool']
            (question, level, ts, klass, binding) = specs[tool]
            r.append(dict(id=pid + ':' + tool, tool=tool, question=question, status=prev['status'], resolution=level, time_scale=ts, decidability=dict(class_=klass, refinement_variable='Same acquired input; numerical sampling/mesh only', justification='Finite digital transport predicates' if klass == 1 else 'Named physical/semantic input is absent; numerical refinement supplies no such observation'), binding_quantity=binding, evidence=[source(oldpath)] if oldpath.exists() else [f.meta['source'] for f in own], theta_sha256=digest(reg), physical_status='NOT_APPLICABLE_DIGITAL_TRANSPORT' if klass == 1 else 'UNKNOWN', support=dict(patient_id=pid, frame_id=frameid, region='source facets/voxels; consumer aggregation explicit'), inherited=prev, numerical_enclosure='Digital identity only' if klass == 1 else 'Physical/continuum bound MISSING'))
            r[-1]['decidability']['class'] = r[-1]['decidability'].pop('class_')
        r.append(dict(id=pid + ':RegionField', tool='RegionField', question='Does address-bound query reflect the same digital source region?', status='SOURCE_VERIFIED', resolution='PER_POINT', time_scale='SIMULTANEOUS', decidability={'class': 1, 'refinement_variable': 'Finite original source support, no anatomical refinement claim', 'justification': 'All crop centres or declared surface query points verified; unknown at ambiguous boundaries'}, binding_quantity='Anatomical accuracy still UNKNOWN', evidence=[f.meta['source'] for f in own], theta_sha256=digest(reg), physical_status='NOT_APPLICABLE_DIGITAL_QUERY', support=dict(patient_id=pid, frame_id=frameid, region='native source-addressed points'), numerical_enclosure='Exact owner at checked voxel centres; surface floating error independently checked, global floating envelope MISSING'))
        if pid == 'Bite2Text_F4775':
            for (tool, question, status, evidence) in [('PoseContact', 'Does the discrete loading surface survive each declared pose?', contact['contract_status'], str(out.parent / 'C1/CONTACT_RESULTS.json')), ('RigidOffset', 'Is there a rigid vertical adjustment that meets both geometrical requirements?', design['status'], str(out.parent / 'C4/ROBUST_DESIGN_RESULTS.json'))]:
                r.append(dict(id=pid + ':' + tool, tool=tool, question=question, status=status, resolution='PER_POINT', time_scale='SIMULTANEOUS', decidability={'class': 1, 'refinement_variable': 'Fixed digital point/facet geometry, exact rational predicates', 'justification': 'Complete pose event partition or exact necessary-inequality counterexample decides the digital question'}, binding_quantity='Loaded registration/compliance and force for a physical answer', evidence=[source(evidence)], theta_sha256=digest(reg), physical_status='UNKNOWN', support=dict(patient_id=pid, frame_id=frameid, region='generated tooth36 source support'), numerical_enclosure='Rational digital predicates; physical registration and FE bounds MISSING'))
        for row in r:
            validate_decision(row, reg)
            for key in ['resolution', 'time_scale', 'decidability', 'binding_quantity', 'theta_sha256']:
                bad = copy.deepcopy(row)
                bad.pop(key)
                checks.append(check(row['id'] + ':requires_' + key, True, rejects(validate_decision, bad, reg)))
            bad = copy.deepcopy(row)
            bad['theta_sha256'] = '0' * 64
            checks.append(check(row['id'] + ':theta_mutation', True, rejects(validate_decision, bad, reg)))
            bad = copy.deepcopy(row)
            bad['physical_status'] = 'CERTIFIED'
            checks.append(check(row['id'] + ':no_physical_promotion', True, rejects(validate_decision, bad, reg)))
        if ios:
            edges += [edge('IOS_SOURCE', 'RegionField', 'PER_POINT', 'PER_POINT', time_scale='SIMULTANEOUS', patient_id=pid, source_support=[f.meta['field_id'] for f in own]), edge('RegionField', 'FE_load', 'PER_POINT', 'PER_TOOTH', time_scale='SIMULTANEOUS', patient_id=pid, source_support='source facet/point plus theta; no physical force attached')]
        else:
            edges += [edge('CBCT_PULP_SOURCE', 'RegionField', 'PER_POINT', 'PER_POINT', time_scale='SIMULTANEOUS', patient_id=pid, source_support=own[0].meta['field_id']), edge('RegionField_preparation_state', 'X33_thermal_initial_state', 'PER_POINT', 'PER_POINT', time_scale='HANDOVER', patient_id=pid, source_support=own[0].meta['field_id'])]
        contract = dict(schema='Patient360-decision-v2', patient_id=pid, claim_type='capability', frame_id=frameid, field_descriptors=[p for p in field_result['descriptors'] if Path(p).name.startswith(pid)], shared_theta=reg, decisions=r, release='ABSTAIN', release_binding_quantities=sorted({row['binding_quantity'] for row in r if row['physical_status'] == 'UNKNOWN'}), review_state='PENDING_INDEPENDENT_REVIEW', actual_modalities=sorted({f.meta['modality'] for f in own}), missing_modalities=['CBCT' if ios else 'IOS'])
        dest = out / 'patients' / pid / 'DECISION_CONTRACT.json'
        dump(dest, contract)
        side = dict(schema='Patient360-export-sidecar-v2', patient_id=pid, frame_id=frameid, unit='mm', theta_sha256=digest(reg), theta=reg, decision_contract=source(dest), source_fields=[dict(field_id=f.meta['field_id'], source_sha256=f.meta['source_sha256'], arrays=f.meta['arrays']) for f in own], release='ABSTAIN', physical_manufacturing_status='NOT_RUN')
        sidefile = dest.parent / 'EXPORT_SIDECAR.json'
        dump(sidefile, side)
        sideback = load(sidefile)
        checks.append(check(pid + ':export_theta_roundtrip', sideback['theta_sha256'] == digest(sideback['theta']), digest({**sideback['theta'], 'patient_id': 'synthetic-wrong-identity'}) != sideback['theta_sha256']))
        patients.append(dict(patient_id=pid, contract=str(dest), sidecar=str(sidefile), decision_count=len(r)))
        allrows += r
    for e in edges:
        validate_edge(e)
        bad = dict(e)
        bad['resolution'] = 'POPULATION'
        checks.append(check(e['patient_id'] + ':finest_edge:' + e['producer'], True, rejects(validate_edge, bad)))
    a = np.array([[-1.0, 1.0], [1.0, -1.0]])
    b = np.array([[-1.0, -1.0], [1.0, 1.0]])
    suff = dict(summary_identity_error=float(np.max(np.abs(np.sort(a, axis=0) - np.sort(b, axis=0)))), identical_marginal_intervals=[[-1.0, 1.0], [-1.0, 1.0]], theta_rows_A=a, theta_rows_B=b, downstream_joint_feasible_A=bool(np.any(np.all(a >= 0, axis=1))), downstream_joint_feasible_B=bool(np.any(np.all(b >= 0, axis=1))), minimum_extension='Preserve joint theta row identity/coupling; interval marginals alone cannot answer joint existence', resolution='PER_POINT', evidence_kind='our_own_fixture')
    x26 = imported('r2_x26', RESULTS / 'LANE_X26_DECIDABILITY/decidability.py')
    margin = 0.1 - float(__import__('fractions').Fraction(design['point_gap_witness']['exact_gap']))
    x26_result = x26.ScalarPort(margin=margin, unit='mm', bias_bound=0.05, coefficient=0).evaluate(0)
    checks.append(check('X26_scalar_shared_pose_floor', x26_result['class'] == 2 and x26_result['conditional_status'] == 'ABSTAIN', x26.interval_decision(-0.01, 0.01, 0) != 'ABOVE', dict(injection='Claim ABOVE for an interval straddling zero', port=x26_result)))
    result = dict(claim_type='capability', patients=patients, decision_count=len(allrows), decisions=allrows, edges=edges, checks=checks, sufficient_summary_test=suff, X26_scalar_check=x26_result, shared_theta_execution=shared_execution, physical_releases=0, external_referent=dict(kind='our_own_fixture', locator=str(RESULTS / 'LANE_X26_DECIDABILITY/decidability.py'), compared_quantity='Decision-relative interval abstention and floor classification; local research code, not independent physical validation', refutes_us=False), runtime_seconds=time.perf_counter() - t0)
    dump(out / 'CONTRACT_RESULTS.json', result)
    return result
