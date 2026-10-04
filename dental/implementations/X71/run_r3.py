import copy, datetime, tempfile, time
from planning import *
from measurement_port import validate
verify_lock()
start = time.perf_counter()
reg = read('PREREG_R3.json')
if reg['input_lock_sha256'] != sha(HERE / 'INPUT_LOCK.json'):
    raise ValueError('PREREG_INPUT_DRIFT')
ids = []
ledger = []
for design in ['D1', 'M1', 'M2']:
    for angle in [0, 30]:
        for n in range(1, 13):
            sid = f'{design}_{angle}_{n:02d}'
            ids.append(sid)
            ledger.append(dict(specimen_id=sid, design=design, angle_deg=angle, final_endpoint='fracture', state_history=['AS_BUILT', 'DRY', 'CEMENTED', 'LOADED', 'FRACTURED'], pilot=n <= 2, noninterference_validation='UNKNOWN'))
fpath = HERE / 'FROZEN_PREDICTIONS.json'
if not fpath.exists():
    contract = next((x for x in read('SOURCE_MANIFEST.json') if x['path'].endswith('X1B/LAB_PROTOCOL.md')))
    write('FROZEN_PREDICTIONS.json', dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), freeze_scope='Prospective acquisition/allocation contract; no absolute physical predictions', claim_type='information_link', input_lock_sha256=sha(HERE / 'INPUT_LOCK.json'), specimen_ids=ids, geometry_source_sha256=contract['sha256'], geometry_binding_scope='Protocol source hash; raw CAD/STL hashes must be added in a new freeze before new lab pieces', Q_primary=0.4743756, Q_primary_window=[0.351426, 0.64034], separable_Q=1.0, separable_window=[0.740818, 1.349859], absolute_force_N='UNKNOWN', retention_aging_ratio='UNKNOWN', costs='PHENOMENOLOGICAL', physical_measurements='NOT_PERFORMED', prediction_source_locator=contract['path'], prediction_scope='Inherited rounded X1B protocol contrasts; reuse requires exact geometry/batch/support/contact checks'))
freeze = read('FROZEN_PREDICTIONS.json')
fs = sha(fpath)
write('FROZEN_PREDICTIONS.sha256', dict(sha256=fs))
write('SPECIMEN_LEDGER.json', ledger)
alternative = copy.deepcopy(ledger)
for x in alternative:
    if int(x['specimen_id'].split('_')[-1]) > 8:
        x['final_endpoint'] = 'retention'
        x['state_history'] = ['AS_BUILT', 'DRY', 'CEMENTED', 'AGED', 'PULLED_OFF']

def summary(xs):
    return [len(xs)] + [sum((x['design'] == d and x['angle_deg'] == a for x in xs)) for d in ['D1', 'M1', 'M2'] for a in [0, 30]]

def fracture_slots(xs):
    return sum((x['final_endpoint'] == 'fracture' for x in xs))

def endpoint_gate(xs):
    return all((len([x for x in xs if x['design'] == d and x['angle_deg'] == a and (x['final_endpoint'] == 'fracture')]) == 12 for d in ['D1', 'M1', 'M2'] for a in [0, 30]))

def ownership(xs):
    return len({x['specimen_id'] for x in xs}) == len(xs) and all((not ('PULLED_OFF' in x['state_history'] and 'FRACTURED' in x['state_history']) for x in xs))
double = copy.deepcopy(ledger)
double[0]['state_history'].insert(-1, 'PULLED_OFF')
dup = copy.deepcopy(ledger)
dup[1]['specimen_id'] = dup[0]['specimen_id']
template = dict(quantity='fracture', unit='N', measurement_state='FRACTURED', resolution_level='PER_TOOTH', specimen_id=ids[0], freeze_sha256=fs, geometry_sha256=freeze['geometry_source_sha256'], value=None, absolute_error_bound=None, raw_path='', raw_sha256='', acquired_at_utc='', registration_locator='', calibration_locator='', operator_id='', fracture_origin='unknown', provenance_kind='independent_measurement')
write('MEASUREMENT_TEMPLATE.json', template)
mutations = []
with tempfile.TemporaryDirectory(prefix='X71_fixture_', dir=HERE / 'raw') as tmp:
    rawp = Path(tmp) / 'fixture_force.json'
    base = copy.deepcopy(template)
    base.update(value=100.0, absolute_error_bound=2.0, raw_path=str(rawp), acquired_at_utc=(datetime.datetime.fromisoformat(freeze['timestamp_utc']) + datetime.timedelta(days=1)).isoformat(), registration_locator='fixture_registration', calibration_locator='fixture_calibration', operator_id='fixture', provenance_kind='our_own_fixture', fracture_origin='intaglio_tensile_zone')
    rawp.write_text(json.dumps({k: base[k] for k in ['value', 'unit', 'specimen_id', 'measurement_state', 'quantity', 'resolution_level']}))
    base['raw_sha256'] = sha(rawp)
    nominal = validate(base, freeze, fs, test_only=True)
    for (key, bad) in [('unit', 'MPa'), ('measurement_state', 'DRY'), ('value', 110.0), ('raw_sha256', '0' * 64), ('freeze_sha256', '0' * 64), ('specimen_id', 'foreign_specimen'), ('geometry_sha256', '0' * 64), ('calibration_locator', ''), ('registration_locator', ''), ('acquired_at_utc', '2000-01-01T00:00:00+00:00'), ('absolute_error_bound', -1), ('fracture_origin', 'assumed_tensile'), ('resolution_level', 'POPULATION')]:
        mutant = copy.deepcopy(base)
        mutant[key] = bad
        mutations.append(dict(field=key, nominal_status=nominal['status'], mutant_result=validate(mutant, freeze, fs, test_only=True)))
    disguised = copy.deepcopy(base)
    disguised['provenance_kind'] = 'independent_measurement'
    ready_freeze = dict(freeze, physical_geometry_ready=True)
    fake_result = validate(disguised, ready_freeze, fs, test_only=False)
    unknown_origin = copy.deepcopy(base)
    unknown_origin['fracture_origin'] = 'unknown'
    unknown_result = validate(unknown_origin, freeze, fs, test_only=True)
out = dict(round='R3', claim_type='information_link', outcome='STATEFUL_ACQUISITION_CONTRACT_DELIVERED_LAB_NOT_RUN', sufficiency=dict(summary_A=summary(ledger), summary_B=summary(alternative), summary_identity_error=max((abs(a - b) for (a, b) in zip(summary(ledger), summary(alternative)))), fracture_endpoint_slots_A=fracture_slots(ledger), fracture_endpoint_slots_B=fracture_slots(alternative), downstream_difference=fracture_slots(ledger) - fracture_slots(alternative), full_Q_protocol_A=endpoint_gate(ledger), full_Q_protocol_B=endpoint_gate(alternative), resolution_level='POPULATION of planned independent specimens', external_referent=dict(kind='our_own_fixture', locator='SPECIMEN_LEDGER.json; X47 and X1B independent protocols', compared_quantity='state-aware planned fracture endpoint ownership, not observed fracture/retention force', refutes_us=True), smallest_extension='One destructive endpoint/state tag per specimen ID; full observation-state sequence needed for measurement reuse'), tests=dict(nominal_test_only_accepts=nominal['status'] == 'TEST_ONLY_ACCEPTED', live_record_mutations=mutations, all_live_mutations_rejected=all((x['mutant_result']['status'] == 'REJECT' for x in mutations)), double_destructive_reuse_rejected=not ownership(double), duplicate_id_rejected=not ownership(dup), baseline_ownership_accepts=ownership(ledger), retention_alternative_is_valid_separate_protocol=ownership(alternative), retention_alternative_rejected_as_unchanged_Q=not endpoint_gate(alternative), disguised_fixture_rejected=fake_result['status'] == 'REJECT', unknown_origin_stays_unknown=unknown_result.get('mechanism_status') == 'UNKNOWN', empty_real_port=validate(template, freeze, fs)), same_information_control=dict(summary_A_control=[72, 12, 12, 12, 12, 12, 12], fracture_slots_control_A=6 * 12, fracture_slots_control_B=6 * 8), physical_measurements_performed=0, freeze_sha256=fs, runtime_s=time.perf_counter() - start)
write('raw/RESULTS_R3.json', out)
write('raw/RETENTION_ALTERNATIVE_LEDGER.json', alternative)
if not all([out['tests']['all_live_mutations_rejected'], out['tests']['double_destructive_reuse_rejected'], out['tests']['disguised_fixture_rejected'], out['tests']['unknown_origin_stays_unknown'], out['sufficiency']['summary_identity_error'] == 0]):
    raise SystemExit('R3_CONTROL_FAILED')
print('summary identity0; fracture endpoint slots72 vs48; difference24;13 live record mutations rejected')
checkpoint('R3_COMPLETE', 'EXACT_SUMMARY_COUNTEREXAMPLE_AND_LIVE_RECORD_MUTATIONS_PASS', 'R4: actual readiness versus crown72-only and extra-fixture profiles')
