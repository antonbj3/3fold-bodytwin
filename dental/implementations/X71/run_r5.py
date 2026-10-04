from planning import *
import time
start = time.perf_counter()
verify_lock()
reg = read('PREREG_R5.json')
if reg['input_lock_sha256'] != sha(HERE / 'INPUT_LOCK.json'):
    raise ValueError('PREREG_INPUT_DRIFT')
ledger = read('SPECIMEN_LEDGER.json')
allids = {x['specimen_id'] for x in ledger}
pilot = {x['specimen_id'] for x in ledger if x['pilot']}
other = {x['specimen_id'] for x in ledger if int(x['specimen_id'].split('_')[-1]) in [3, 4]}

def records(scan, film, force):
    return [{'port': port, 'specimen_id': sid, 'state': state} for (port, ids, state) in [('scan', scan, 'AS_BUILT'), ('film', film, 'CEMENTED'), ('force', force, 'FRACTURED')] for sid in sorted(ids)]
a = records(pilot, pilot, pilot)
b = records(pilot, other, pilot)

def summary(rows):
    return [(port, sum((x['port'] == port for x in rows))) for port in ['scan', 'film', 'force']]

def usable(rows):
    expected = {'scan': 'AS_BUILT', 'film': 'CEMENTED', 'force': 'FRACTURED'}
    by = {port: {x['specimen_id'] for x in rows if x['port'] == port and x['state'] == state} for (port, state) in expected.items()}
    return set.intersection(*by.values())

def full72(rows):
    return usable(rows) == allids
wrong = [dict(x) for x in a]
wrong[12]['state'] = 'DRY'
foreign = [dict(x) for x in a]
foreign[0]['specimen_id'] = 'different_specimen'
full = records(allids, allids, allids)
out = dict(round='R5', claim_type='information_link', outcome='SPECIMEN_COVERAGE_INDEX_REQUIRED_FULL72_COST_NOT_IDENTIFIED', sufficiency=dict(summary_A=summary(a), summary_B=summary(b), summary_identity_error=0, downstream_joint_coverage_A=len(usable(a)), downstream_joint_coverage_B=len(usable(b)), downstream_difference=len(usable(a)) - len(usable(b)), resolution_level='PER_TOOTH identities across PER_POINT field acquisitions', smallest_extension='Specimen identity plus observation state/support on each port; retain intersection rather than counts', external_referent=dict(kind='our_own_fixture', locator='raw/COVERAGE_COUNTEREXAMPLE.json; source X1B LAB_PROTOCOL.md', compared_quantity='joint coverage of fixed planned specimens, not measured geometry or strength', refutes_us=True)), tests=dict(nominal_pilot_joint12=len(usable(a)) == 12, disjoint_film_joint0=len(usable(b)) == 0, wrong_state_rejected=len(usable(wrong)) == 11, foreign_id_rejected=len(usable(foreign)) == 11, pilot12_rejected_as_complete72=not full72(a), full_nominal72_accepts=full72(full), direct_set_control_matches=usable(a) == pilot and (not usable(b))), full72_cost_status='UNKNOWN: R2 M05 row36–60h prices only the12-piece metrology prefix. Additional60 shape/film/setup records need a new cost quote, so it is not a full-protocol cost estimate.', impact='Crown72-only first-week pilot M01/M02/M04 still uses identical12 identities and keeps5 potential local queries; full72 completion and its cost were never established.', external_protocol_locator=next((x['path'] for x in read('SOURCE_MANIFEST.json') if x['path'].endswith('X1B/LAB_PROTOCOL.md'))), physical_measurements_performed=0, runtime_s=time.perf_counter() - start)
write('raw/RESULTS_R5.json', out)
write('raw/COVERAGE_COUNTEREXAMPLE.json', dict(state_A=a, state_B=b, wrong_state=wrong, foreign_id=foreign))
if not all(out['tests'].values()):
    raise SystemExit('R5_COVERAGE_FAILURE')
checkpoint('R5_COMPLETE', 'EXACT_PORT_COUNT_IDENTITY_WITH_JOINT_COVERAGE12_VS0', 'Finalize with full72 cost UNKNOWN; require same-specimen actual acquisition and cost logs')
print('identical port counts; joint coverage12 vs0; full72 setup missing60 metrology/film records')
