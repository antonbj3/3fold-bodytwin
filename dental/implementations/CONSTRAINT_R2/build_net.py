"""Extend R1 read-only. A missing residual is not a zero residual."""
from dental_release.paths import expand as _release_expand
import copy
from collections import Counter
from decimal import Decimal
from common import *
CORRECTED = [f'D-E-CONTRAST-{i:02}' for i in range(1, 9)] + ['D-E-CEMENT-TRANSFER-FAIL', 'D-E-BRIDGE-SYSTEM']

def decimal(x):
    return Decimal(str(x))

def make_gap(quantity, unit, level, operands, operation='reported_residual', scope='', scale=None, kind='residual_magnitude'):
    vals = [decimal(r['value']) for r in operands]
    if operation == 'reported_residual':
        signed = vals[0]
    elif operation == 'difference':
        signed = vals[0] - vals[1]
    elif operation == 'complement100':
        signed = Decimal(100) - sum(vals)
    elif operation == 'maximum':
        signed = max((abs(x) for x in vals))
    else:
        raise ValueError(operation)
    return {'value': float(abs(signed)), 'signed_value': float(signed), 'exact_decimal': str(abs(signed)), 'unit': unit, 'quantity': quantity, 'resolution_level': level, 'kind': kind, 'operation': operation, 'operands': operands, 'scope': scope, 'normalization_scale': scale, 'physical_transfer_residual': 'UNKNOWN', 'interpretation': 'Numerical residual for this stated observable only; no unseen-specimen error bound.', 'rigorous_enclosure': 'Exact decimal arithmetic on reported operands; measurement/model enclosure inherited, not newly proved.'}

def missing_gap(e, variables):
    return {'value': None, 'unit': variables[e['outputs'][0]]['unit'] if e.get('outputs') else 'UNSPECIFIED', 'quantity': e['constraint'], 'resolution_level': e['resolution_level'], 'kind': 'UNIDENTIFIED_RESIDUAL', 'operation': None, 'operands': [], 'reason': 'No commensurate observed response and known contribution for the full stated relation in the frozen sources.', 'required_measurement': e.get('knowledge_debt', {}).get('replacement_measurement', e['constraint']), 'source_of_missingness': e.get('evidence', []), 'not_complete': True}

def intake():
    lock = read(HERE / 'INPUT_LOCK.json')
    rows = []
    for job in lock['jobs']:
        reviews = []
        for p in lock['files']:
            if Path(p).name == 'REVIEW_' + job + '.json':
                r = source(p)
                actual = lock['files'].get(r['result_file'], {}).get('sha256')
                reviews.append({'record': p, 'record_sha256': lock['files'][p]['sha256'], 'result_file': r['result_file'], 'expected_sha256': r.get('result_sha256'), 'actual_sha256': actual, 'hash_match': actual == r.get('result_sha256'), 'decision': r.get('graph_decision', r.get('decision')), 'corrections': r.get('correction', [])})
        eligible = any((r['hash_match'] and r['decision'] in ['ACCEPT', 'ACCEPT_WITH_CORRECTION'] for r in reviews))
        has_result = source_path('results/' + job + '/results.json') in lock['files'] or any((r['actual_sha256'] for r in reviews))
        rows.append({'job': job, 'eligible_reviewed': eligible, 'reviews': reviews, 'disposition': 'REVIEWED' if eligible else 'DEFERRED_REVIEW_PENDING' if has_result else 'DEFERRED_NO_FINAL_RESULT', 'reason': None if eligible else 'Completed matching independent review absent at freeze; no admission from SAMPLE or producer status.'})
    return rows

def main():
    net = source('results/PROOF_LANE_CONSTRAINT_NET_DENTAL/CONSTRAINT_NET_DENTAL.json')
    net.update(schema_version='dental-constraint-net/r2', seed='PROOF_LANE-constraint-net-r2', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', objective='Sourced residuals, corrected observation aggregation, unit-safe partial stress and X71 decision incidence')
    original_ids = {e['id'] for e in net['edges']}
    vars = net['variables']
    byid = {e['id']: e for e in net['edges']}
    corrections = []
    for id_ in CORRECTED:
        e = byid[id_]
        old = e['resolution_level']
        e['anatomical_support_level'] = old
        e['resolution_level'] = e['evidence_resolution_level'] = 'POPULATION'
        e['observation_aggregation'] = 'Published group means / group contrasts; not individual teeth'
        for r in e['evidence']:
            r['anatomical_support_level'] = r['resolution_level']
            r['resolution_level'] = 'POPULATION'
        corrections.append({'edge': id_, 'before': old, 'after': 'POPULATION', 'source': 'XREV18-C04'})
    for e in net['edges']:
        e['r2_origin'] = 'INHERITED_R1'
        e['decision_ids'] = sorted(set(e['consumer_chains']))
        e['decision_count'] = len(e['decision_ids'])
        e['decision_count_scope'] = 'Distinct direct K questions, not independent measurements or all transitive descendants'
        if e['status'] in ['OPEN', 'UNKNOWN']:
            e['gap'] = missing_gap(e, vars)
        e['gap_components'] = []

    def ev(job, key, unit, level):
        return evidence('results/' + job + '/results.json', key, unit, level)

    def attach(ids, g, primary=True):
        for id_ in ids:
            e = byid[id_]
            e['gap_components'].append(copy.deepcopy(g))
            if primary:
                e['gap'] = copy.deepcopy(g)
    g = make_gap('maximum held Catphan calibration residual', 'HU_ref', 'PER_SURFACE_REGION', [ev('LANE_X24_CBCT_HU_CALIBRATION', '/key_results/R1_max_holdout_HU/value', 'HU_ref', 'PER_SURFACE_REGION')], scope='Four source CBCT scanners, held insert responses; no bone modulus or anatomical truth calibration')
    attach(['D-E-K02'], g)
    g = make_gap('unmatched published pulp-temperature proxy MAE', 'degC', 'POPULATION', [ev('LANE_X33_LASER_PULP', '/external_validation/proxy_transfer_MAE_C', 'degC', 'POPULATION')], scope='Explicitly unmatched retrospective group comparison; no pointwise transient bound')
    attach(['D-E-PULP-HEAT', 'D-E-K23'], g)
    job = 'LANE_X44_MISSING_CHAINS'
    rows = source('results/' + job + '/results.json')['chains']['K50_R2']['held_points']
    held = []
    for (i, row) in enumerate(rows):
        base = f'/chains/K50_R2/held_points/{i}'
        gg = make_gap('PDL held group force residual', 'N', 'POPULATION', [ev(job, base + '/F_mean_N', 'N', 'POPULATION'), ev(job, base + '/prediction/F_mean_N', 'N', 'POPULATION')], 'difference', scope=f"{row['phase']}, displacement {row['deflection_mm']} mm; digitized group curve, interpolation enclosure missing", scale=ev(job, base + '/tolerance_N', 'N', 'POPULATION'))
        held.append(gg)
    worst = max(held, key=lambda x: x['value'])
    for id_ in ['D-E-K17', 'D-E-K50']:
        attach([id_], worst)
        byid[id_]['residual_field'] = held
    p = 'results/DEMO48_PACKAGE/demos/X14/results.json'
    gg = []
    for i in range(2):
        gg.append(make_gap('published calculated sinter gap reconstruction residual', 'um', 'PER_SURFACE_REGION', [evidence(p, f'/geometric_reconstruction/values_um/{i}', 'um', 'PER_SURFACE_REGION'), evidence(p, f'/geometric_reconstruction/reference_um/{i}', 'um', 'PER_SURFACE_REGION')], 'difference', scope='Published calculated geometry; measured seated-fit residual is UNKNOWN'))
    attach(['D-E-SINTER-FIT'], max(gg, key=lambda x: x['value']), primary=False)
    byid['D-E-SINTER-FIT']['residual_field'] = gg
    p = 'results/DEMO48_PACKAGE/batch15/demos/X59/results.json'
    biases = source(p)['source_bias_replica_minus_ct_um']
    gs = []
    for (manufacturer, regions) in biases.items():
        for region in regions:
            gs.append(make_gap('replica minus dry CT regional group contrast', 'um', 'POPULATION', [evidence(p, f'/source_bias_replica_minus_ct_um/{manufacturer}/{region}', 'um', 'POPULATION')], scope=manufacturer + '; ' + region + '; mixed observation state, no individual agreement or correction factor'))
    attach(['D-E-REPLICA-CT'], max(gs, key=lambda x: x['value']))
    byid['D-E-REPLICA-CT']['residual_field'] = gs
    byid['D-E-REPLICA-CT']['anatomical_support_level'] = 'PER_SURFACE_REGION'
    byid['D-E-REPLICA-CT']['resolution_level'] = byid['D-E-REPLICA-CT']['evidence_resolution_level'] = 'POPULATION'
    for r in byid['D-E-REPLICA-CT']['evidence']:
        r['anatomical_support_level'] = r['resolution_level']
        r['resolution_level'] = 'POPULATION'
    p = 'results/DEMO48_PACKAGE/batch16/demos/X63/results.json'
    g = make_gap('preload loss after repeated tightening', 'N', 'POPULATION', [evidence(p, '/summary/source_force_first_mean_N', 'N', 'POPULATION'), evidence(p, '/summary/source_force_tenth_mean_N', 'N', 'POPULATION')], 'difference', scope='Group means, first minus tenth screw tightening; no implant insertion torque or new-system prediction')
    attach(['D-E-PRELOAD-HISTORY'], g)
    byid['D-E-PRELOAD-HISTORY']['anatomical_support_level'] = 'PER_TOOTH'
    byid['D-E-PRELOAD-HISTORY']['resolution_level'] = byid['D-E-PRELOAD-HISTORY']['evidence_resolution_level'] = 'POPULATION'
    for r in byid['D-E-PRELOAD-HISTORY']['evidence']:
        r['anatomical_support_level'] = r['resolution_level']
        r['resolution_level'] = 'POPULATION'
    ledger = intake()
    eligible = {r['job']: r for r in ledger}
    new = []
    pending = []

    def var(id_, desc, unit, level):
        if id_ not in vars:
            vars[id_] = {'desc': desc, 'unit': unit, 'resolution_level': level, 'value': None, 'chains': [], 'value_status': 'SCOPED_OBSERVATIONS_ONLY', 'working_bindings': [], 'knowledge_debt': {'status': 'UNKNOWN', 'replacement_measurement': 'Matched specimen ' + desc}}
    var('relative_tooth_force', 'Tooth relative force share', 'percentage_points', 'PER_TOOTH')
    var('fracture_closure_mean', 'Weibull-derived expected group fracture load', 'N', 'POPULATION')
    var('measurement_plan', 'Named measurement and specimen obligations', '1', 'PHENOMENOLOGICAL')
    var('contact_count', 'Named tooth observed contact count', 'count', 'PER_TOOTH')
    var('crown_anatomy_residual', 'Full crown surface reconstruction residual', 'mm', 'PER_TOOTH')
    var('bond_strength_group', 'Coupon nominal shear bond strength group mean', 'MPa', 'POPULATION')

    def add(id_, job, between, ks, constraint, gap=None, level='PER_POINT', time='SIMULTANEOUS', fields=None, keys=()):
        source_name = 'results/' + job + '/results.json'
        d = source(source_name) if source_path(source_name) in read(HERE / 'INPUT_LOCK.json')['files'] else {}
        item = eligible[job]
        e = {'id': id_, 'between': between, 'inputs': between[:-1], 'outputs': between[-1:], 'consumer_chains': ks, 'decision_ids': ks, 'decision_count': len(set(ks)), 'decision_count_scope': 'Distinct direct K consumers; no clinical utility interpretation', 'constraint': constraint, 'status': 'OPEN' if gap else 'UNKNOWN', 'claim_type': 'information_link', 'resolution_level': level, 'evidence_resolution_level': level, 'timescale': time, 'time_contract': {'input': 'source-defined same specimen/protocol', 'output': 'same time' if time == 'SIMULTANEOUS' else 'terminal state passed to next process'}, 'evidence': [ev(job, key, '1', level) for key in keys], 'source_review_state': item['disposition'], 'review': item['reviews'], 'external_referent': d.get('external_referent', {'kind': 'external_review', 'locator': item['reviews'][0]['record'] if item['reviews'] else source_name, 'compared_quantity': constraint, 'refutes_us': True}), 'rigorous_enclosure': {'status': 'MISSING_FOR_PHYSICAL_TRANSFER', 'reason': 'No unmeasured physical error bound from source summaries'}, 'knowledge_debt': {'status': 'UNKNOWN', 'replacement_measurement': 'Same specimen/protocol observation of ' + vars[between[-1]]['desc']}, 'working_bindings': [], 'gap_components': [], 'r2_origin': job, 'scope': constraint, 'no_scientific_admission': True}
        e['gap'] = gap or missing_gap(e, vars)
        if fields:
            e['residual_field'] = fields
        if item['eligible_reviewed']:
            net['edges'].append(e)
            new.append(id_)
            byid[id_] = e
        else:
            pending.append(e)
        return e
    job = _release_expand('X69')
    d = source('results/' + job + '/results.json')
    for (j, jaw) in enumerate(d['r2']['jaws']):
        fs = []
        for (i, row) in enumerate(jaw['tooth_results']):
            fs.append(make_gap('paired IOS-CBCT tooth p95 surface residual', 'mm', 'PER_TOOTH', [ev(job, f'/r2/jaws/{j}/tooth_results/{i}/p95_mm', 'mm', 'PER_TOOTH')], scope=f"Hao2023 Demo1 {jaw['jaw']} FDI{row['fdi']}; fit-held={row['heldout_from_fit']}; STL mm conditional, physical TRE UNKNOWN"))
        e = add('D-E-X69-' + jaw['jaw'].upper(), job, ['scan_surface', 'anatomy_surface', 'pose', 'signed_error'], ['K05'], 'Paired surface residual after fitted rigid registration; retain every tooth and held/fit identity', max(fs, key=lambda x: x['value']), level='PER_TOOTH', fields=fs, keys=('/status', '/r2/heldout_gate_pass'))
    byid['D-E-K05']['numeric_component_edges'] = ['D-E-X69-LOWER', 'D-E-X69-UPPER']
    job = _release_expand('X70')
    d = source('results/' + job + '/results.json')
    for (i, row) in enumerate(d['external_force']['rows']):
        q = 'empirical group mean fracture load residual' if i < 3 else 'Weibull-closure expected fracture load residual'
        g = make_gap(q, 'N', 'POPULATION', [ev(job, f'/external_force/rows/{i}/observed_N', 'N', 'POPULATION'), ev(job, f'/external_force/rows/{i}/predicted_N', 'N', 'POPULATION')], 'difference', scope=row['locator'] + '; retrospective held-source comparison; XREV18-C01 preserves empirical vs closure distinction')
        e = add(f'D-E-X70-FORCE-{i + 1}', job, ['wall_thickness', 'fracture_mean' if i < 3 else 'fracture_closure_mean'], ['K36', 'K49'], q, g, level='POPULATION')
        e['quantity_origin'] = 'EMPIRICAL_PUBLISHED_GROUP_MEAN' if i < 3 else 'CONSTITUTIVE_CLOSURE_WEIBULL_EXPECTATION'
        e['comparison_error'] = ev(job, f'/external_force/rows/{i}/error_log', '1', 'POPULATION')
        e['comparison_tolerance'] = ev(job, '/external_force/tolerance_log', '1', 'POPULATION')
        e['external_referent'] = {'kind': 'independent_measurement' if i < 3 else 'closed_form', 'locator': row['locator'], 'compared_quantity': q, 'refutes_us': True}
    add('D-E-X71-PLAN', _release_expand('X71'), ['asbuilt_error', 'cement_gap', 'fracture_mean', 'measurement_plan'], ['K43', 'K48', 'K49'], 'M01+M02+M04 makes five specified observation questions evaluable on the same 12-crown pilot; no observed physical outcome', level='PER_TOOTH', time='HANDOVER', keys=('/first_week_crown72_only', '/physically_new_decisions', '/full_chains_closed'))
    job = _release_expand('X72')
    p = 'results/' + job + '/raw/X18.json'
    fs = []
    for (i, row) in enumerate(source(p)['rows']):
        fs.append(make_gap('conditional local friction sandwich width', 'MPa', 'PER_TOOTH', [evidence(p, f'/rows/{i}/coulomb_peak_interval_MPa/1', 'MPa', 'PER_TOOTH'), evidence(p, f'/rows/{i}/coulomb_peak_interval_MPa/0', 'MPa', 'PER_TOOTH')], 'difference', scope=f"case {row['case']} FDI{row['fdi']}, 20 N scenario; geometric cone sandwich only; FE rounding/physical bound missing", kind='conditional_enclosure_width'))
    add('D-E-X72-FRICTION', job, ['contact_gap', 'patch_force', 'friction', 'stress'], ['K18', 'K21'], 'Width between inner and outer local-friction stress optimizations; retain site allocations and scenario', max(fs, key=lambda x: x['value']), level='PER_TOOTH', fields=fs, keys=('/rigorous_enclosure',))
    job = _release_expand('FALT_TANDLAST')
    p = 'results/LANE_X2_OCCLUSION_B2B/raw/EXTERNAL_REFERENTS.json'
    ref = source(p)['studies'][0]
    g = make_gap('unreported relative force mass on 16-FDI support', 'percentage_points', 'POPULATION', [evidence(p, f'/studies/0/mean_percent/{i}', 'percentage_points', 'POPULATION') for i in range(len(ref['mean_percent']))], 'complement100', scope='Ferrato Table3 reports 14 upper FDI; remaining FDI18/28 unobserved. Arithmetic remainder only; Table5 totals 99.9%, no joint mechanical completion claimed.')
    e = add('D-E-TANDLAST-PARTIAL-REFERENCE', job, ['contact_gap', 'relative_tooth_force'], ['K16', 'K21'], 'XREV18-C03: 14-coordinate reference cannot refute normalization on 16 coordinates', g, level='POPULATION')
    e['missing_fdi'] = [18, 28]
    e['anatomical_support_level'] = 'PER_TOOTH'
    add('D-E-X73-WALL', _release_expand('X73'), ['anatomy_surface', 'wall_sigma', 'nerve_distance'], ['K04', 'K06', 'K37'], 'Dense same-image release revisions affect clearance; independent annotator physical floor remains unidentified', keys=('/independent_spread', '/unique_tooth_consumer_result'))
    job = _release_expand('X74')
    d = source('results/' + job + '/results.json')
    fs = []
    for (i, row) in enumerate(d['chains']['R05']['covered_system_screens']):
        fs.append(make_gap('coupon challenge minus nominal SBS group mean', 'MPa', 'POPULATION', [ev(job, f'/chains/R05/covered_system_screens/{i}/lab_challenge_MPa', 'MPa', 'PHENOMENOLOGICAL'), ev(job, f'/chains/R05/covered_system_screens/{i}/mean_MPa', 'MPa', 'POPULATION')], 'difference', scope=row['system'] + '; nominal coupon group; challenge is an engineering assumption, no cohesive patient law'))
    add('D-E-X74-BOND', job, ['material_history', 'bond_strength_group'], ['K13', 'K34'], 'Cure/cover protocol to nominal coupon bond group observations; full local law unknown', level='POPULATION', time='HANDOVER', fields=fs, keys=('/chains/R05/patient_result',))
    add('D-E-X75-BONE', _release_expand('X75'), ['cbct_gray', 'cortex_thickness', 'bone_E', 'insertion_torque'], ['K12', 'K26'], 'Apparent same-site contrast profiles do not calibrate density or torque', keys=('/variance_explained_torque', '/variance_explained_temperature', '/physical_cortex_density_accuracy'))
    job = _release_expand('FULL_CROWN')
    d = source('results/' + job + '/results.json')
    fs = []
    for (i, row) in enumerate(d['rounds']['R7']['rows']):
        if 'reconstruction_p95_mm' in row:
            fs.append(make_gap('full crown reconstruction p95 residual', 'mm', 'PER_TOOTH', [ev(job, f'/rounds/R7/rows/{i}/reconstruction_p95_mm', 'mm', 'PER_TOOTH')], scope=row['key'] + '; sampled bidirectional p95 includes virtual reference closure; continuous p95 enclosure missing'))
    add('D-E-FULL-CROWN-ANATOMY', job, ['anatomy_surface', 'design_choice', 'crown_anatomy_residual'], ['K36'], 'R7 full crown anatomy error on fixed v5 metric; every scored tooth retained', max(fs, key=lambda x: x['value']), level='PER_TOOTH', fields=fs)
    job = _release_expand('GENCAD_V6')
    b = '/rounds/R4/largest_reported_difference'
    g = make_gap('force-share change unexplained by identical contact count', 'percentage_points', 'PER_TOOTH', [ev(job, b + '/A/force_pp', 'percentage_points', 'PER_TOOTH'), ev(job, b + '/B/force_pp', 'percentage_points', 'PER_TOOTH')], 'difference', scope='Same reported count 2 at FDI11; time-specific source observations; not patient force prediction')
    add('D-E-GENCAD-V6-FORCE', job, ['contact_count', 'relative_tooth_force'], ['K16'], 'Equal observed contact count does not identify tooth relative force', g, level='PER_TOOTH', keys=(b + '/A', b + '/B'))
    for e in net['edges']:
        for v in e['between']:
            for k in e['consumer_chains']:
                if k not in vars[v]['chains']:
                    vars[v]['chains'].append(k)
        for c in net['chains']:
            if c['id'] in e['consumer_chains'] and e['id'] not in c['edge_ids']:
                c['edge_ids'].append(e['id'])
    net['corrections_applied'] = corrections
    net['additional_aggregation_repairs'] = ['D-E-REPLICA-CT', 'D-E-PRELOAD-HISTORY']
    net['pending_unreviewed_edges'] = pending
    net['post_r1_intake'] = ledger
    net['source_files'] = sorted(set(net['source_files']) | set(read(HERE / 'INPUT_LOCK.json')['files']))
    net['r2_source_lock'] = 'INPUT_LOCK.json'
    net['granularity_policy'] = 'Relations retain finest common supported observation. Group aggregation is POPULATION, anatomical support is a separate field. Residual fields retain each supplied site/phase.'
    net['gap_policy'] = 'Finite sourced scoped residuals only; unidentified quantities stay null and MUST fail strict validation. No invented scalar closes a physical chain.'
    net['dropout'] = {'new_deliveries_considered': len(ledger), 'reviewed_retained': sum((r['eligible_reviewed'] for r in ledger)), 'deferred': sum((not r['eligible_reviewed'] for r in ledger)), 'fraction_deferred': sum((not r['eligible_reviewed'] for r in ledger)) / len(ledger), 'details': ledger, 'edge_candidates': len(new) + len(pending), 'reviewed_edges_added': len(new), 'pending_edges': len(pending)}
    dump('CONSTRAINT_NET_DENTAL_R2.json', net)
    dump('CORRECTIONS_APPLIED.json', corrections)
    dump('INTAKE_LEDGER.json', ledger)
    dump('GAP_ACQUISITION_LEDGER.json', [{'edge': e['id'], 'decision_ids': e['decision_ids'], 'gap': e['gap']} for e in net['edges'] if e['status'] in ['OPEN', 'UNKNOWN'] and e['gap']['value'] is None])
    print(json.dumps({'variables': len(vars), 'edges': len(net['edges']), 'new_reviewed': len(new), 'pending_edges': len(pending), 'status': dict(Counter((e['status'] for e in net['edges'])))}))
if __name__ == '__main__':
    main()
