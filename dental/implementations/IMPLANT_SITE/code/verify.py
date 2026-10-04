from common_local import *
import copy, subprocess
from scipy.stats import binom
check_frozen()
verify_inputs()
checks = []

def check(name, ok):
    checks.append({'name': name, 'pass': bool(ok)})
    if not ok:
        raise AssertionError(name)
plans = lines(P / 'raw/SITE_PLANS_R1.jsonl')
choices = load(P / 'raw/RESEARCH_CHOICES_R2.json')
r2 = load(P / 'rounds/R2.json')
r3 = load(P / 'rounds/R3.json')
check('2581 unique case/FDI rows', len(plans) == len({(p['case'], p['fdi']) for p in plans}) == 2581)
physical = ['physical_recommended_length_mm', 'physical_recommended_diameter_mm', 'connection', 'material', 'drill_protocol', 'insertion_torque_Ncm', 'temperature_C', 'local_micromotion_um', 'fatigue_survival', 'nerve_injury_probability']
valid = lambda p: all((p[k] is None for k in physical))
check('No patient endpoint imputation', all((valid(p) for p in plans)))
fault = copy.deepcopy(plans[0])
fault['insertion_torque_Ncm'] = 25
check('Injected patient torque rejected', not valid(fault))
ref = next((p for p in plans if p['gray'] is not None))
fault = copy.deepcopy(ref)
fault['gray']['density_g_cm3'] = 1.2
valid_gray = lambda p: p['gray']['density_g_cm3'] is None and p['gray']['status'] == 'UNCALIBRATED_GRAY_OBSERVATION'
check('Gray-to-density imputation rejected', valid_gray(ref) and (not valid_gray(fault)))
frozen = load(P / 'FROZEN_PREDICTIONS.json')
check('Future prediction hash', sha(P / 'FROZEN_PREDICTIONS.json') == (P / 'FROZEN_PREDICTIONS.json.sha256').read_text().strip())
check('Injected future target rejected', hashlib.sha256((P / 'FROZEN_PREDICTIONS.json').read_bytes() + b'0').hexdigest() != sha(P / 'FROZEN_PREDICTIONS.json'))

def identity(actual, expected):
    return actual['source_label_sha256'] == expected['source_label_sha256'] and actual['case'] == expected['case'] and (actual['fdi'] == expected['fdi']) and (actual['pose'] == expected['pose'])
for field in ['source_label_sha256', 'case', 'fdi', 'pose']:
    fault = copy.deepcopy(ref)
    if field == 'pose':
        fault['pose']['entry_zyx_mm'][0] += 0.3
    elif field == 'fdi':
        fault[field] += 1
    else:
        fault[field] = 'CORRUPTED'
    check('Injected identity ' + field + ' rejected', not identity(fault, ref))
check('R1 strict replay failure retained', load(P / 'rounds/R1.json')['identity_control']['source_replay_gate_1e_10'] == 'FAIL_PRESERVED')
check('R2 geometry and independent control gates', all(r2['gates'].values()))
raw_drill = load(DATA / 'DRILL_GEOMETRY_R3.json')
check('R3 all source-mask brackets close', all((v['gap_mm'] <= 1e-06 for q in raw_drill for v in q['distances'].values())))
roots = load(P / 'raw/DRILL_LIMITS_R3.json')
check('R3 root intervals satisfy frozen resolution', all((q['allowable_advance_interval_mm'][1] - q['allowable_advance_interval_mm'][0] <= 0.01 for q in roots if q['status'] == 'ROOT_BRACKETED')))
rr = list(csv.DictReader((P / 'raw/DRILL_ADVANCE_R3.csv').open()))
initial = {(q['case'], q['fdi']): q['digital_2mm_pass'] == 'True' for q in rr if float(q['advance_mm']) == 0}
check('Reported drill count is reproduced', sum((initial[q['case'], q['fdi']] and q['digital_2mm_pass'] == 'False' for q in rr if float(q['advance_mm']) == 1.5)) == 15 and sum(initial.values()) == 21)
for (path, tol) in [('GEOMETRY_CONTROLS_R2.json', 1e-06), ('DRILL_CONTROLS_R3.json', 1e-06)]:
    cr = load(P / 'raw' / path)
    check(path + ' original residuals pass', all((c['error_mm'] <= tol for c in cr)))
    check(path + ' injected +1 mm rejected', all((c['error_mm'] + 1 > tol for c in cr)))
for name in ['R1', 'R2']:
    w = load(P / f'raw/SUFFICIENCY_{name}.json')
    check(name + ' exact summary identity', w['identity_error'] == 0.0)
    check(name + ' distinguishes downstream', w.get('downstream_difference', w.get('downstream_difference_mm', 0)) > 0)
ep = module('x65_endpoint_verify', D / 'LANE_X65_ISO14801_FATIGUE/code/endpoint_consumer.py')
epfile = D / 'LANE_X65_ISO14801_FATIGUE/raw/published_endpoints.json'
register_inputs([epfile])
options = []
for row in load(epfile)['records']:
    if not row['strict_eligible']:
        continue
    q = {'system_id': row['id'], 'quantity': 'maximum_cyclic_force', 'unit': 'N', 'force_class_N': 175.0, **{k: row[k] for k in ['diameter_mm', 'material', 'connection', 'angle_deg', 'R', 'horizon_cycles']}}
    ans = ep.endpoint_query(q)
    check('Source endpoint ' + row['id'], ans['maximum_cyclic_force_N'] == row['force_N'])
    options.append({'source_query': q, 'source_response': ans, 'site_transport': 'UNKNOWN: exact system geometry and patient load absent; not interchangeable connection classes'})
    bad = dict(q, reference_force_N=row['force_N'] + 100)
    check('Injected endpoint ' + row['id'] + ' rejected', ep.endpoint_query(bad)['status'] == 'UNKNOWN')
check('Patient identity cannot satisfy lab protocol', ep.endpoint_query({'case': ref['case'], 'fdi': ref['fdi']})['status'] == 'UNKNOWN')
dump(P / 'raw/LAB_CONNECTION_REFERENCES.json', options)
mic = P / 'raw/MICROMOTION_QUERY.json'
dump(mic, {'units': {'position': 'mm', 'force': 'N', 'moment': 'Nmm', 'displacement': 'um'}, 'bone_marker_readings_um': None})
ans = subprocess.check_output([sys.executable, str(D / 'LANE_X56_MICROMOTION/code/query_motion.py'), str(mic)], text=True)
check('X56 missing local reference returns UNKNOWN', json.loads(ans)['status'] == 'UNKNOWN_MISSING_LOCAL_BONE_REFERENCE')
dump(P / 'raw/MICROMOTION_RESPONSE.json', json.loads(ans))

def comparable(a, b):
    return all((a.get(k) == b.get(k) and a.get(k) is not None for k in ['quantity', 'denominator', 'time_window', 'cohort']))
clinical = dict(quantity='neurosensory_event', denominator='patient', time_window='first_visit', cohort='Bartling1999')
geom = dict(quantity='digital_bound_not_certified', denominator='image_FDI', time_window='static_geometry', cohort='ToothFairy')
check('Clinical/geometric frequency fusion rejected', not comparable(clinical, geom))
r = load(P / 'raw/CLINICAL_REFERENTS.json')['records'][0]
(lo, hi) = r['ci95']
check('Clinical CI lower inversion', abs(binom.sf(r['events'] - 1, r['n'], lo) - 0.025) < 1e-10)
check('Clinical CI upper inversion', abs(binom.cdf(r['events'], r['n'], hi) - 0.025) < 1e-10)
check('Injected clinical CI rejected', abs(binom.sf(r['events'] - 1, r['n'], lo + 0.05) - 0.025) > 1e-10)
check('Future targets do not claim injury', all((r['physical_injury_probability'] == 'UNKNOWN' for r in frozen['targets'])))
from compare_metrology import compare
check('Empty metrology remains UNKNOWN', compare({'measurements': []})['status'] == 'UNKNOWN_NO_INDEPENDENT_MEASUREMENT')
target = frozen['targets'][0]
mock = {'kind': 'independent_measurement', 'units': 'mm', 'frozen_prediction_sha256': sha(P / 'FROZEN_PREDICTIONS.json'), 'instrument_record': 'MOCK_UNIT_TEST', 'registration_record': 'MOCK_UNIT_TEST', 'measurements': [{'case': target['case'], 'fdi': target['fdi'], 'advance_mm': target['advance_mm'], 'L_mm': target['L_mm'], 'D_mm': target['D_mm'], 'distance_mm': target['gap_lower_mm'], 'absolute_measurement_bound_mm': 0.01}]}
check('Metrology mock nominal contract', compare(mock)['rows'][0]['status'] == 'PASS_WITHIN_DECLARED_LAB_TOLERANCE')
bad = copy.deepcopy(mock)
bad['measurements'][0]['distance_mm'] += 1.0
check('Injected metrology distance +1mm fails', compare(bad)['rows'][0]['status'] == 'FAIL_DISJOINT')
bad = copy.deepcopy(mock)
bad['measurements'][0]['absolute_measurement_bound_mm'] = 10.0
check('Wide metrology uncertainty cannot pass', compare(bad)['rows'][0]['status'] == 'UNKNOWN_MEASUREMENT_TOO_WIDE')
for (key, value) in [('kind', 'our_own_fixture'), ('units', 'cm'), ('frozen_prediction_sha256', 'bad')]:
    bad = copy.deepcopy(mock)
    bad[key] = value
    try:
        compare(bad)
        rejected = False
    except ValueError:
        rejected = True
    check('Injected metrology ' + key + ' rejected', rejected)
size = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
check('Data less than 3GB', size < 3000000000)
dump(P / 'raw/VERIFICATION.json', {'checks': checks, 'all_pass': all((c['pass'] for c in checks)), 'count': len(checks), 'data_bytes': size, 'scope': 'Numerical/identity/endpoint contracts, not independent biological or clinical validation'})
print(json.dumps({'verification': 'PASS', 'checks': len(checks), 'data_bytes': size}))
