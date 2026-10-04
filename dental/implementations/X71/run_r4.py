import copy, time
from planning import *
verify_lock()
start = time.perf_counter()
reg = read('PREREG_R4.json')
if reg['input_lock_sha256'] != sha(HERE / 'INPUT_LOCK.json'):
    raise ValueError('PREREG_INPUT_DRIFT')
ds = read('inputs/DECISIONS.json')
ms = read('inputs/MEASUREMENTS.json')
ass = read('inputs/LAB_ASSUMPTIONS.json')
crown_ids = {'M01', 'M02', 'M03', 'M04', 'M05'}
allp = list(portfolios(ds, ms))

def select(available, scratch=False):
    candidates = []
    for p in allp:
        if not set(p['ids']) <= set(available):
            continue
        x = copy.deepcopy(p)
        if scratch and p['ids']:
            x['operator_hours'] = [v + ass['shared_manufacture_operator_hours'][i] for (i, v) in enumerate(p['operator_hours'])]
            x['instrument_hours'] = [v + ass['shared_manufacture_machine_hours'][i] for (i, v) in enumerate(p['instrument_hours'])]
            x['shared_sinter_elapsed_hours'] = ass['shared_sinter_elapsed_hours']
        if x['operator_hours'][1] <= 40 and x['instrument_hours'][1] <= 40:
            candidates.append(x)
    candidates.sort(key=lambda p: (-len(p['evaluation_ready']), p['new_specimens'], p['operator_hours'][1], p['ids']))
    return candidates[0]
ready = select(crown_ids)
extra = select(crown_ids | {'M08'})
scratch = select(crown_ids, True)
unknown = select([])
template = dict(status='NOT_MEASURED', currency_cost='UNKNOWN', operator_hours_available='UNKNOWN', ready_specimen_ids=[], instrument_access=[], verified_cost_intervals={}, matched_foursupport_specimens='UNKNOWN', calibration_and_registration_evidence='UNKNOWN', specimen_noninterference='UNKNOWN', actual_CAD_and_batch_hashes='UNKNOWN', note='Fill with independently checked lab inventory/time quotes; unknown fields do not authorize physical readiness. No human input was requested in this autonomous run.')
write('LAB_READINESS_TEMPLATE.json', template)
out = dict(round='R4', claim_type='information_link', outcome='READY_SPECIMEN_CONDITIONAL_WEEK_PLAN_ACTUAL_READINESS_UNKNOWN', profiles=dict(ready_crown72_only=ready, ready_crown72_and_matched_foursupport=extra, crown72_manufacture_from_scratch=scratch, actual_readiness_unknown=unknown), test=dict(unknown_readiness_returns_no_plan=unknown['ids'] == [], original72_retention_conflict_excluded='M06' not in ready['ids'], scratch_high_cost_rejects_full_plan=scratch['ids'] == [], eligible_measurements_have_available_profile_resources=all((x in crown_ids for x in ready['ids'])), readiness_live_mutation_rejects_m08='M08' not in select(crown_ids)['ids'], same_information_control=[p for p in allp if p['ids'] == ready['ids']][0]['evaluation_ready'] == ready['evaluation_ready']), actual_lab_inventory='UNKNOWN', actual_costs='UNKNOWN', physically_new_decisions=0, limits='Ready profiles assume exact matched fixtures, independent calibration and available instruments; upper cost intervals are estimates and need quotes. Every actual response can remain undecidable at a threshold.', runtime_s=time.perf_counter() - start)
write('raw/RESULTS_R4.json', out)
if not all(out['test'].values()):
    raise SystemExit('R4_READINESS_CONTROL_FAILED')
print('crown72-only', ready['ids'], 'potential queries', len(ready['evaluation_ready']), 'operator_h', ready['operator_hours'])
print('actual readiness UNKNOWN; from-scratch conservative40h plan unavailable')
checkpoint('R4_COMPLETE', 'READINESS_CONTROLS_PASS_PHYSICAL_READINESS_UNKNOWN', 'Finalize one-command artifacts, scoped graph feedback and laboratory handoff')
