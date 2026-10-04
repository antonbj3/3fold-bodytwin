import csv, re, time
from planning import *
start = time.perf_counter()
verify_lock()
reg = read('PREREG_R1.json')
if reg['input_lock_sha256'] != sha(HERE / 'INPUT_LOCK.json'):
    raise ValueError('PREREG_INPUT_DRIFT')
decisions = read('inputs/DECISIONS.json')
measures = read('inputs/MEASUREMENTS.json')
seen = {}
witness = None
n = 0
for p in portfolios(decisions, measures):
    n += 1
    key = (len(p['touched']), *p['operator_hours'], *p['instrument_hours'], p['new_specimens'])
    if key in seen and len(seen[key]['evaluation_ready']) != len(p['evaluation_ready']):
        witness = dict(summary=key, summary_identity_error=0, state_A=seen[key], state_B=p, downstream_difference=abs(len(seen[key]['evaluation_ready']) - len(p['evaluation_ready'])))
        break
    seen.setdefault(key, p)
man = read('SOURCE_MANIFEST.json')
src = next((x for x in man if x['path'].endswith('X59/raw/cunali2017.txt')))
primary_text = Path(src['path']).read_text()
lines = [x for x in primary_text.splitlines() if re.search('(Replica|Micro-CT)\\s+\\d+\\.\\d+±', x)]
primary = []
for line in lines:
    row = re.findall('(\\d+\\.\\d+)±(\\d+\\.\\d+)', line)
    primary.append(row)
if len(primary) != 4 or any((len(r) != 4 for r in primary)):
    raise ValueError('PRIMARY_TABLE_PARSE')
rows = read('inputs/EXTERNAL_FACIT.json')

def table_gate(imported):
    errors = []
    for (i, x) in enumerate(imported):
        a = 2 * (i // 4)
        j = i % 4
        expected = float(primary[a][j][0]) - float(primary[a + 1][j][0])
        observed = float(x['mean_1_um']) - float(x['mean_2_um'])
        errors.extend([abs(observed - expected), abs(float(x['sd_1_um']) - float(primary[a][j][1])), abs(float(x['sd_2_um']) - float(primary[a + 1][j][1]))])
    return (max(errors) <= 0.011, max(errors))
(good, err) = table_gate(rows)
bad = [dict(x) for x in rows]
bad[0]['mean_1_um'] = str(float(bad[0]['mean_1_um']) + 10)
(bad_pass, baderr) = table_gate(bad)
out = dict(round='R1', claim_type='information_link', candidate='affected count and summed interval cost sufficient', outcome='COUNT_SUMMARY_REFUTED' if witness else 'NO_COUNTEREXAMPLE_IN_FROZEN_DOMAIN', witness=witness, portfolios_examined_until_witness=n, smallest_extension='No minimum established by R1. Cost tuple may encode the finite portfolio identity accidentally; test specimen-state summary in R3. Named required-port incidence remains needed for scope.', primary_table=dict(locator='doi:10.1590/0103-6440201601531 Table1, local primary text', passed=good, max_error_um=err, live_plus10um_rejected=not bad_pass, bad_error_um=baderr, compared_quantity='replica minus dry microCT regional group means and reported SD; PER_SURFACE_REGION', mean_differences_um=[float(x['mean_1_um']) - float(x['mean_2_um']) for x in rows]), physical_measurements_performed=0, runtime_s=time.perf_counter() - start)
write('raw/RESULTS_R1.json', out)
if not good or bad_pass:
    raise SystemExit('TABLE_CONTROL_FAILED')
write('HANDOFF_R1.json', dict(outcome=out['outcome'], next_operation='R2 retain AND prerequisites, interval cost, and specimen compatibility', frozen_prereg_sha256=sha(HERE / 'PREREG_R1.json')))
print(out['outcome'], 'identity', witness['summary_identity_error'] if witness else None, 'difference', witness['downstream_difference'] if witness else None)
checkpoint('R1_COMPLETE', out['outcome'], 'R2: dependency-aware full portfolio with declared cost intervals')
