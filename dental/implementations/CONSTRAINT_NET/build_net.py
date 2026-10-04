"""Build a reviewed quantity network, never admit evidence to source graphs."""
import copy, datetime, json
from collections import Counter
from common import *
from definitions import VARIABLES, CHAINS, HANDOVER, BUILT
X43 = 'results/LANE_X43_NORTHSTAR_STATUS/results.json'
REVIEW43 = 'results/LANE_XREVIEW_BATCH8/REVIEW_LANE_X43_NORTHSTAR_STATUS.json'
LEVEL_ORDER = ['PER_POINT', 'PER_SURFACE_REGION', 'PER_TOOTH', 'PER_ARCH', 'POPULATION', 'PHENOMENOLOGICAL']

def review_contract(review_name, result_name):
    r = src(review_name)
    return {'record': str(path(review_name)), 'record_sha256': SOURCES[str(path(review_name))]['sha256'], 'decision_key': '/graph_decision' if 'graph_decision' in r else '/decision', 'decision': r.get('graph_decision', r.get('decision')), 'reviewed_result': str(path(result_name)), 'reviewed_result_sha256': r.get('result_sha256'), 'scope': r.get('review_scope', ''), 'corrections': r.get('correction', []), 'scientific_admission': False}

def main():
    variables = {}
    for row in VARIABLES.strip().splitlines():
        (id_, desc, unit, level, nid) = row.split('|')
        variables[id_] = {'desc': desc, 'unit': unit, 'resolution_level': level, 'value': None, 'value_status': 'NO_GLOBAL_VALUE_USE_SCOPED_EDGE_OBSERVATIONS', 'working_bindings': [working_binding(nid)], 'chains': [], 'instance_key': ['specimen_or_population', 'region_or_point', 'protocol', 'time', 'source_lineage'], 'uncertainty': 'No patient-level distribution inferred from aggregate observations', 'knowledge_debt': {'status': 'UNKNOWN', 'replacement_measurement': 'Matched observation of ' + desc}}
    chains = []
    edges = []
    src_chains = {c['k']: (i, c) for (i, c) in enumerate(src(X43)['chains'])}
    for row in CHAINS.strip().splitlines():
        (num, ports, nid, formula) = row.split('|')
        k = 'K' + num
        (ins, outs) = [s.split() for s in ports.split(' -> ')]
        between = list(dict.fromkeys(ins + outs))
        (i, c) = src_chains[k]
        finest = min((variables[v]['resolution_level'] for v in between), key=LEVEL_ORDER.index)
        status = 'OPEN' if k in BUILT else 'UNKNOWN'
        e = {'id': 'D-E-' + k, 'between': between, 'inputs': ins, 'outputs': outs, 'constraint': formula, 'status': status, 'status_reason': 'Executed scoped suboperation documented; full matched physical relation not validated' if status == 'OPEN' else 'Required relation or physical closure not established in the admitted source subset', 'relation_role': 'CHAIN_CONTRACT', 'claim_type': 'capability', 'consumer_chains': [k], 'resolution_level': finest, 'evidence_resolution_level': 'PHENOMENOLOGICAL', 'timescale': 'HANDOVER' if k in HANDOVER else 'SIMULTANEOUS', 'time_contract': {'input': 'declared specimen state before operation' if k in HANDOVER else 'same declared state/time', 'output': 'operation terminal state passed to subsequent process' if k in HANDOVER else 'same declared state/time'}, 'evidence': [evidence(X43, f'/chains/{i}/chain', 'PHENOMENOLOGICAL', role='consumer_definition'), evidence(X43, f'/chains/{i}/remaining_requirement', 'PHENOMENOLOGICAL', role='known_missing_input'), evidence(X43, f'/chains/{i}/rationale', 'PHENOMENOLOGICAL', role='reviewed_suboperation_scope')], 'review': review_contract(REVIEW43, X43), 'working_bindings': [working_binding(nid)], 'external_referent': {'kind': 'external_review', 'locator': str(path(REVIEW43)), 'compared_quantity': 'K-chain scope and remaining requirement, not the physical map', 'refutes_us': True}, 'uncertainty': c['uncertainty'], 'rigorous_enclosure': {'status': 'MISSING_FOR_FULL_PHYSICAL_CHAIN', 'reason': c['remaining_requirement']}, 'knowledge_debt': {'status': 'UNKNOWN', 'replacement_measurement': c['remaining_requirement']}, 'scope': 'Quantity-type contract. Instantiating shared patient/specimen state requires matching instance_key; no numerical fusion.', 'no_clinical_recommendation': True}
        edges.append(e)
        chains.append({'id': k, 'desc': c['chain'], 'variables': between, 'inputs': ins, 'outputs': outs, 'edge_ids': [e['id']], 'working_ids': [nid], 'remaining_requirement': c['remaining_requirement'], 'source_evidence': e['evidence'][0], 'physical_chain_closed': False})
        for v in between:
            variables[v]['chains'].append(k)

    def add(id_, between, ks, result, review, keys, level, unit, formula, referent, scope, status='TIGHT', role='SCOPED_OBSERVATION', time='SIMULTANEOUS', unc='No transfer uncertainty inferred', enclosure='NOT_A_DERIVATIVE_OR_CONTINUUM_PREDICTION'):
        result = 'results/' + result + '/results.json'
        review = 'results/' + review
        e = {'id': id_, 'between': between, 'inputs': between[:-1], 'outputs': between[-1:], 'constraint': formula, 'status': status, 'status_reason': 'Exact reviewed external comparison within the stated scope' if status == 'TIGHT' else 'Built relationship with missing matched physical reference', 'relation_role': role, 'claim_type': 'capability', 'consumer_chains': ks, 'resolution_level': level, 'evidence_resolution_level': level, 'timescale': time, 'time_contract': {'input': 'source-defined measurement protocol', 'output': 'source-defined measurement/assay time; see scope'}, 'evidence': [evidence(result, key, level, unit) for key in keys], 'review': review_contract(review, result), 'external_referent': referent, 'scope': scope, 'uncertainty': unc, 'working_bindings': list({b['id']: b for v in between for b in variables[v]['working_bindings']}.values()), 'rigorous_enclosure': {'status': enclosure, 'reason': 'Only the declared finite source comparison is asserted; no unbounded affine extrapolation.'}, 'knowledge_debt': {'status': 'UNKNOWN', 'replacement_measurement': 'Same-specimen/protocol ' + ', '.join((variables[v]['desc'] for v in between))}, 'no_clinical_recommendation': True}
        e['external_comparison'] = {'executed': status == 'TIGHT', 'result_file': str(path(result)), 'quantity_paths': keys, 'support': 'reviewed result and explicitly scoped source comparison', 'source_result_sha256': SOURCES[str(path(result))]['sha256']}
        edges.append(e)
        for c in chains:
            if c['id'] in ks:
                c['edge_ids'].append(id_)
                for v in between:
                    if v not in c['variables']:
                        c['variables'].append(v)
                    if c['id'] not in variables[v]['chains']:
                        variables[v]['chains'].append(c['id'])

    def ext(kind, locator, qty, refutes=False):
        return dict(kind=kind, locator=locator, compared_quantity=qty, refutes_us=refutes)
    r24 = 'LANE_X24_CBCT_HU_CALIBRATION'
    rv24 = 'LANE_XREVIEW_BATCH3/REVIEW_JOB_LANE_X24_CBCT_HU_CALIBRATION.json'
    add('D-E-PHANTOM-HU', ['cbct_gray', 'hu_ref'], ['K02', 'K12'], r24, rv24, ['/key_results/R1_max_holdout_HU/value', '/key_results/R2_max_holdout_HU/value'], 'PER_SURFACE_REGION', 'HU_ref', 'Held-out max |calibrated(gray)-HU_ref| for two frozen calibration constructions', src('results/' + r24 + '/results.json')['external_referent'], 'Four scanners, Catphan insert regional responses; failure of the calibration target. No conversion from polymer HU_ref to bone E.', role='NEGATIVE_EXTERNAL_COMPARISON', unc='Observed held-out maximum, not a confidence bound')
    add('D-E-HU-ORDER-OBSTRUCTION', ['cbct_gray', 'hu_ref'], ['K02'], r24, rv24, ['/key_results/monotone_order_bound/0', '/key_results/monotone_order_bound/1'], 'PER_SURFACE_REGION', 'HU_ref', 'Reversed input order gives minimax monotone output error at least (HU_LDPE-HU_PMP)/2=58 HU_ref', src('results/' + r24 + '/results.json')['external_referent'], 'Two observed insert pairs on Elekta devices; bound refutes 40 HU target for a scalar monotone map on those pairs.', role='NEGATIVE_EXTERNAL_COMPARISON', enclosure='EXACT_FINITE_ORDER_BOUND')
    r36 = 'LANE_X36_META_REGRESSION'
    rv36 = 'LANE_XREVIEW_BATCH6/REVIEW_LANE_X36_META_REGRESSION_CROWN.json'
    contrasts = src('results/' + r36 + '/results.json')['within_study_contrasts']
    for (i, c) in enumerate(contrasts):
        cement = c['dataset'] == 'cement'
        vs = ['cement_setting', 'marginal_gap'] if cement else ['wall_thickness', 'fracture_mean']
        level = 'PER_SURFACE_REGION' if cement else 'PER_TOOTH'
        add(f'D-E-CONTRAST-{i + 1:02}', vs, ['K43', 'K44'] if cement else ['K36', 'K49'], r36, rv36, [f'/within_study_contrasts/{i}'], level, c['unit'], 'Finite two-arm contrast: ' + ('(gap_high-gap_low)/(setting_high-setting_low)*10' if cement else 'log(F_mean_high/F_mean_low)/log(t_high/t_low)'), ext('independent_measurement', 'https://doi.org/' + c['doi'] + ' ; ' + c['source_locator'], 'Study-specific marginal gap contrast' if cement else 'Study-specific mean full-system crown fracture load thickness contrast', c['local_proxy_excluded']), c['stratum'] + '; group means only, conditional independent-arm interval; not a local derivative, full law or individual fracture quantile.', unc='Conditional source interval in evidence; population sampling unit, specimen support ' + level)
    add('D-E-CEMENT-TRANSFER-FAIL', ['cement_setting', 'marginal_gap'], ['K43'], r36, rv36, ['/holdout/cement/summary'], 'PER_SURFACE_REGION', 'um', 'Held-source transfer RMSE compared with fixed-setting practice', src('results/' + r36 + '/results.json')['external_referent'], 'Retrospective external reported group means; study-balanced errors reject this transfer model, do not validate another.', role='NEGATIVE_EXTERNAL_COMPARISON')
    r28 = 'LANE_X28_HEAT_OSSEO_HANDOVER'
    rv28 = 'LANE_XREVIEW_BATCH3/REVIEW_JOB_LANE_X28_HEAT_OSSEO_HANDOVER.json'
    add('D-E-DRILL-CULTURE', ['drill_rpm', 'cell_viability', 'cell_yield'], ['K24', 'K31'], r28, rv28, ['/decisive_observed_answers/viability_500minus1000_pp', '/decisive_observed_answers/yield14_500over1000', '/assay_time_basis'], 'POPULATION', 'mixed: percentage_points; ratio; assay_time', 'Canine protocol-group contrasts: viability(500rpm)-viability(1000rpm); yield14(500rpm)/yield14(1000rpm)', src('results/' + r28 + '/results.json')['external_referent'], 'Canine bone chips; viability assay day3 and two-week culture outgrowth. Culture day is not postoperative healing day. rpm groups also carry protocol covariates; no isolated causal rpm effect.', time='HANDOVER')
    add('D-E-CULTURE-NOT-ISQ', ['cell_viability', 'cell_yield', 'isq'], ['K31'], r28, rv28, ['/decisive_observed_answers/predicted_ISQ_weeks', '/decisive_observed_answers/predicted_healing_delay_weeks'], 'POPULATION', 'ISQ;week', 'No quantified cell-assay to human ISQ(t) map', src('results/' + r28 + '/results.json')['external_referent'], 'Explicit missing handover from biological assay to human longitudinal stability.', status='UNKNOWN', time='HANDOVER')
    r44 = 'LANE_X44_MISSING_CHAINS'
    rv44 = 'LANE_XREVIEW_BATCH8/REVIEW_LANE_X44_MISSING_CHAINS.json'
    add('D-E-REGION-FIELD', ['region_id', 'anatomy_surface', 'sdf'], ['K09', 'K10', 'K13'], r44, rv44, ['/chains/K09/regions', '/chains/K09/ownership_mismatch_count', '/chains/K09/distance_max_errors_mm'], 'PER_POINT', 'mixed: mm;mm2;mm3;count', 'Named region signed distance/owner reproduces source geometry; volumes and interfaces checked separately against published VTK', src('results/' + r44 + '/results.json')['chains']['K09']['external_referent'], 'One published OFJ mandible, source piecewise-linear geometry only; point query errors plus regional integrals. Physical anatomy error remains unknown.', role='SOURCE_GEOMETRY_COMPARISON')
    add('D-E-PDL-CALIBRATION', ['pdl_E', 'tooth_force', 'tooth_displacement'], ['K13', 'K17'], r44, rv44, ['/chains/K13/conditional_E_interval_MPa', '/chains/K13/calibrated_F_at_0p15_N', '/chains/K13/conditioned_on', '/chains/K13/uncertainty'], 'PER_TOOTH', 'mixed: MPa;N;mm', 'Match incisor force at 0.15 mm using archived FE force-versus-modulus curves', src('results/' + r44 + '/results.json')['chains']['K13']['external_referent'], 'Upper-incisor effective population anchor; nu=0.45 and source geometry/rate. Interpolation remainder UNKNOWN; no patient or mandibular material modulus.', status='OPEN', unc='Conditional interval is not a rigorous physical enclosure', enclosure='MISSING_INTERPOLATION_REMAINDER')
    add('D-E-PDL-OBSERVATION', ['tooth_displacement', 'tooth_force', 'material_history'], ['K17', 'K50'], r44, rv44, ['/chains/K50_R2/held_pass_count', '/chains/K50_R2/held_count', '/chains/K50_R2/max_error_N', '/chains/K50_R2/protocol_guard/material_promotion'], 'POPULATION', 'N;count', 'Phase-conditioned interpolation tested at ten held digitized force-deflection points', src('results/' + r44 + '/results.json')['chains']['K50_R2']['external_referent'], 'Observed loading/unloading group curves only. ±0.3 N digitization. All points seen in legacy source before split; material interpretation refused.', unc='Observed held errors and digitization interval; continuum interpolation enclosure missing', enclosure='MISSING_CONTINUUM_REMAINDER')
    r26 = 'LANE_X26_DECIDABILITY'
    rv26 = 'LANE_XREVIEW_BATCH3/REVIEW_JOB_LANE_X26_DECIDABILITY.json'
    add('D-E-CANAL-NUMERIC-BUDGET', ['wall_sigma', 'mesh_pitch', 'numerical_error', 'nerve_distance'], ['K04', 'K06', 'K37', 'K42'], r26, rv26, ['/canal_conditional/floor_mm', '/canal_conditional/sigma_w_mm', '/canal_conditional/tau_mm', '/canal_conditional/numerical_operator', '/canal_conditional/class1', '/canal_conditional/class2', '/canal_conditional/physical_certificates'], 'PER_SURFACE_REGION', 'mixed: mm;count', 'Strict finite resampling budget: numerical distance error <= sqrt(3)*h; physical floor remains conditional', src('results/' + r26 + '/results.json')['external_referents'][3], '4132 dependent policy questions on 338 scans; scalar Gaussian parameters PHENOMENOLOGICAL, zero physical certificates.', status='OPEN', unc='No simultaneous anatomical surface coverage', enclosure='RIGOROUS_NUMERICAL_ONLY_PHYSICAL_MISSING')
    r31 = 'LANE_X31_CLINICAL_ANSWERS'
    rv31 = 'LANE_XREVIEW_BATCH3/REVIEW_JOB_LANE_X31_CLINICAL_ANSWERS.json'
    add('D-E-GUIDE-CLEARANCE', ['guide_play', 'apex_error', 'wall_sigma', 'nerve_distance'], ['K25', 'K37', 'K42'], r31, rv31, ['/guide_parameters/fully_guided'], 'PER_SURFACE_REGION', 'mm', 'Finite pose/support bound for whole-body clearance under declared profile versus apex-only model', src('results/' + r31 + '/results.json')['external_referent'], 'Rigid-pose engineering scenarios. Published means/SD are not signed clearance quantiles or patient margins.', status='OPEN', unc='Physical calibration UNKNOWN; boundary scenarios do not constitute a confidence interval', enclosure='CONDITIONAL_SUPPORT_BOUND_ONLY')
    r33 = 'LANE_X33_LASER_PULP'
    rv33 = 'LANE_XREVIEW_BATCH4/REVIEW_LANE_X33_LASER_PULP.json'
    add('D-E-PULP-HEAT', ['heat_flux', 'dentin_thickness', 'pulp_temperature'], ['K23'], r33, rv33, ['/external_validation', '/R1/refinement_gate', '/R1/physical_validity'], 'PER_POINT', 'degC', 'Transient heat solver compared to published sensor changes; numerical refinement and unmatched physical transfer remain failed', src('results/' + r33 + '/results.json')['external_referent'], 'Source sensor population mean is not local transient maximum; no physical heat partition calibration.', status='OPEN', enclosure='MISSING_NUMERICAL_AND_PHYSICAL_ENCLOSURE')
    r29 = 'LANE_X29_BRIDGE_CONNECTOR'
    rv29 = 'LANE_XREVIEW_BATCH3/REVIEW_JOB_LANE_X29_BRIDGE_CONNECTOR.json'
    add('D-E-BRIDGE-SYSTEM', ['connector_area', 'fracture_mean'], ['K36', 'K49'], r29, rv29, ['/R2/groups/2/area_mm2', '/R2/groups/2/mean_N', '/R2/groups/2/mean_CI95_N', '/R2/groups/2/connector_capacity'], 'PER_SURFACE_REGION', 'mixed: mm2;N', 'Reported group mean system break load and mean confidence interval at a fixed connector area', src('results/' + r29 + '/results.json')['R1']['external_referent'], 'KATANA STML 4Y 9 mm2, fixed geometry; origin not confirmed in connector. Mean load is not connector capacity or individual lower-tail quantile.')
    package_ports = [('D-E-SINTER-FIT', ['sinter_factor', 'marginal_gap'], ['K43'], 'demos/X14', None, ['/geometric_reconstruction', '/unknowns'], 'PER_SURFACE_REGION', 'um', 'Angular observations reconstruct published calculated gaps, not measured seated fit'), ('D-E-DEFORMABLE-FORCES', ['contact_gap', 'tooth_force', 'patch_force'], ['K16', 'K21'], 'batch14/demos/X54', 'LANE_XREVIEW_BATCH14/REVIEW_LANE_X54_DEFORMABLE_CONTACT.json', ['/decision_change_valid_loaded_cases', '/rigorous_enclosure', '/physical_measurement_status'], 'PER_SURFACE_REGION', '1', 'Contact calculations compared to jaw-matched relative sensor signals; absolute same-specimen patch loads remain unknown'), ('D-E-MICROMOTION-OBSERVATION', ['tooth_force', 'micromotion'], ['K19', 'K27'], 'batch14/demos/X56', 'LANE_XREVIEW_BATCH14/REVIEW_LANE_X56_MICROMOTION.json', ['/measurement_counts', '/per_implant_type_bone_class_interface_threshold'], 'PER_POINT', 'um', 'Loaded system displacements do not identify local interface-slip thresholds'), ('D-E-CANAL-RELEASE', ['anatomy_surface', 'wall_sigma', 'nerve_distance'], ['K04', 'K06', 'K37'], 'batch15/demos/X58', 'LANE_XREVIEW_BATCH15/REVIEW_LANE_X58_CANAL_WALL_SPREAD.json', ['/outcome', '/regional_shared_wall', '/radial_sufficiency_test'], 'PER_POINT', 'mm', 'Same-image release revisions change local canal distances; not independently labelled anatomical truth'), ('D-E-REPLICA-CT', ['cement_gap', 'marginal_gap'], ['K43', 'K48'], 'batch15/demos/X59', 'LANE_XREVIEW_BATCH15/REVIEW_LANE_X59_REPLICA_MICROCT.json', ['/source_bias_replica_minus_ct_um', '/individual_bland_altman', '/facit_raw_specimen_pairs'], 'PER_SURFACE_REGION', 'um', 'Regional replica-versus-CT means mix apparatus and seating state; individual paired agreement unknown'), ('D-E-PRELOAD-HISTORY', ['insertion_torque', 'preload', 'assembly_gap'], ['K20', 'K38'], 'batch16/demos/X63', 'LANE_XREVIEW_BATCH16/REVIEW_LANE_X63_IMPLANT.json', ['/summary/source_force_first_mean_N', '/summary/source_force_tenth_mean_N', '/summary/source_mean_drop_N'], 'PER_TOOTH', 'N', 'Same repeated screw-tightening protocol changes measured group clamp force; not an implant insertion-torque law')]
    variables['tightening_torque'] = copy.deepcopy(variables['insertion_torque'])
    variables['tightening_torque'].update(desc='Abutment screw tightening torque', working_bindings=[working_binding('DENT-IF-IA-PRELOAD')], chains=['K20', 'K38'])
    variables['tightening_torque']['knowledge_debt']['replacement_measurement'] = 'Matched screw torque and axial preload, same friction and tightening history'
    for (id_, vs, ks, rel, rv, keys, level, unit, formula) in package_ports:
        if id_ == 'D-E-PRELOAD-HISTORY':
            vs = ['tightening_torque', 'preload', 'assembly_gap']
        name = 'results/DEMO48_PACKAGE/' + rel + '/results.json'
        d = src(name)
        e = {'id': id_, 'between': vs, 'inputs': vs[:-1], 'outputs': vs[-1:], 'constraint': formula, 'status': 'OPEN', 'status_reason': 'Scoped reviewed demo; missing matched physical transfer or ACCEPT* gate for this exact copied result', 'relation_role': 'PACKAGE_TRANSFER_PORT', 'claim_type': 'capability', 'consumer_chains': ks, 'resolution_level': level, 'evidence_resolution_level': level, 'timescale': 'HANDOVER' if id_ in ['D-E-SINTER-FIT', 'D-E-PRELOAD-HISTORY'] else 'SIMULTANEOUS', 'time_contract': {'input': 'source protocol', 'output': 'source protocol; repeated tightening retained'}, 'evidence': [evidence(name, key, level, unit) for key in keys], 'review': review_contract('results/' + rv, src('results/' + rv)['result_file']) if rv else None, 'package_record': str(path('results/DEMO48_PACKAGE/demos.json')), 'external_referent': d['external_referent'], 'scope': formula, 'uncertainty': 'Local physical predictive uncertainty is not identified by published group summaries', 'working_bindings': list({b['id']: b for v in vs for b in variables[v]['working_bindings']}.values()), 'rigorous_enclosure': {'status': 'MISSING_FOR_TRANSFER', 'reason': 'Source comparisons do not enclose unmeasured target specimens'}, 'knowledge_debt': {'status': 'UNKNOWN', 'replacement_measurement': 'Matched specimen, region and protocol observations for ' + ', '.join(vs)}, 'no_clinical_recommendation': True}
        edges.append(e)
        for c in chains:
            if c['id'] in ks:
                c['edge_ids'].append(id_)
                for v in vs:
                    if v not in c['variables']:
                        c['variables'].append(v)
                    if c['id'] not in variables[v]['chains']:
                        variables[v]['chains'].append(c['id'])
    for e in edges:
        e['parameter_sampling_resolution'] = 'POPULATION' if e['relation_role'] == 'SCOPED_OBSERVATION' else e['evidence_resolution_level']
    used_reviews = {e['review']['record'] for e in edges if e['review']}
    census = read(HERE / 'REVIEW_CENSUS.json')
    rows = []
    for r in census:
        if r['review_file'] in used_reviews:
            reason = 'KEPT_REVIEW_SCOPE'
        elif not str(r['decision']).startswith('ACCEPT'):
            reason = 'NO_ACCEPT_SCIENTIFIC_SCOPE'
        elif r['declared_sha256'] != r['actual_sha256']:
            reason = 'RESULT_HASH_MISMATCH_OR_MISSING'
        else:
            reason = 'OUTSIDE_CURATED_RELATION_SET_OR_DUPLICATE_SCOPE'
        rows.append({**r, 'selection': reason})
    counts = Counter((r['selection'] for r in rows))
    net = {'seed': 'PROOF_LANE-constraint-net / PREREG_R1.json', 'objective': 'K01–K52 quantity-level constraint queries with explicit evidence and unresolved transfer', 'schema_version': 'dental-constraint-net/1.0', 'claim_type': 'capability', 'schema_note': 'Vehicle schema core retained: variables object; edges list with id, between, constraint, status and evidence[{source_file,key,value}]. Dental status is stricter: TIGHT only exact reviewed external comparison, never whole-chain validation by association.', 'source_files': sorted({x['source_file'] for e in edges for x in e['evidence']}), 'variables': variables, 'edges': edges, 'chains': chains, 'flip_list': [], 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'no_scientific_admission': True, 'granularity_policy': 'Edge operates at finest named common support; evidence group means retain population sampling unit. Never manufacture PER_POINT values from a regional/population mean.', 'dropout': {'denominator_kind': 'all top-level LANE_XREVIEW*/REVIEW_*.json records, not all dental claims', 'candidates': len(rows), 'kept': counts['KEPT_REVIEW_SCOPE'], 'rejected_or_deferred': len(rows) - counts['KEPT_REVIEW_SCOPE'], 'fraction': (len(rows) - counts['KEPT_REVIEW_SCOPE']) / len(rows), 'reasons': dict(counts), 'ledger': 'CANDIDATE_LEDGER.json'}, 'stress_map_v1_note': 'CONSTRAINT_STRESS_MAP_DENTAL.json contains unresolved consumer sets. Scalar counts cannot choose sufficient measurement bundles.'}
    for e in edges:
        if e['status'] != 'TIGHT':
            net['flip_list'].append({'id': 'MEASURE-' + e['id'], 'question': e['knowledge_debt']['replacement_measurement'], 'flips': e['id'], 'cost_class': 'UNKNOWN', 'cost_note': 'See X71 prospective assumptions; no empirical measurement-cost claim.'})
    dump('CONSTRAINT_NET_DENTAL.json', net)
    dump('CANDIDATE_LEDGER.json', rows)
    dump('CURATION_CONTRACT.json', {'edges': [{k: e[k] for k in ['id', 'between', 'consumer_chains', 'resolution_level', 'timescale', 'status', 'external_referent', 'scope']} for e in edges]})
    print(json.dumps({'variables': len(variables), 'edges': len(edges), 'chains': len(chains), 'status_counts': dict(Counter((e['status'] for e in edges)))}))
if __name__ == '__main__':
    main()
