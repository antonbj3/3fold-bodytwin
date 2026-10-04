import csv, json, time, math, re
from pathlib import Path
import numpy as np
from physics import ROOT, legacy, draws, forward, PROFILE
from endpoint_consumer import endpoint_query

def save(p, x):
    (ROOT / p).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    data = json.loads((ROOT / 'raw/published_endpoints.json').read_text())
    byid = {r['id']: r for r in data['records']}
    queries = []
    mutations = []
    for aid in ['K1-A3a', 'K1-A3b', 'K1-A3c']:
        row = byid[aid]
        q = {'system_id': aid, **{k: row[k] for k in ['diameter_mm', 'material', 'connection', 'angle_deg', 'R', 'horizon_cycles']}, 'quantity': 'maximum_cyclic_force', 'unit': 'N', 'force_class_N': 175}
        queries.append({'query': q, 'answer': endpoint_query(q)})
        for (k, v) in [('horizon_cycles', 2000000.0), ('angle_deg', 35), ('quantity', 'failure_probability'), ('unit', 'MPa'), ('system_id', 'unmeasured_new_connection')]:
            altered = dict(q)
            altered[k] = v
            ans = endpoint_query(altered)
            mutations.append({'gate': 'EXACT_ENDPOINT_IDENTITY', 'nominal_pass': True, 'injected_error_rejected': ans['status'] == 'UNKNOWN', 'injection': {k: v}})
        bad = dict(q)
        del bad['connection']
        mutations.append({'gate': 'EXACT_ENDPOINT_IDENTITY', 'nominal_pass': True, 'injected_error_rejected': endpoint_query(bad)['status'] == 'UNKNOWN', 'injection': 'remove connection'})
        bad = dict(q)
        bad['reference_force_N'] = 2 * row['force_N']
        mutations.append({'gate': 'PRIMARY_FORCE_VALUE', 'nominal_pass': True, 'injected_error_rejected': endpoint_query(bad)['status'] == 'UNKNOWN', 'injection': 'double asserted primary force'})
    tables = json.loads((ROOT / 'raw/primary_tables.json').read_text())
    source = tables['PMC10820107']
    tt = next((t for t in source['tables'] if t['id'] == 'materials-17-00434-t002'))
    cell = next((row[1] for row in tt['cells'] if row and row[0] == 'V'))
    (mean, sd) = [float(re.search('-?\\d+(?:\\.\\d+)?', s).group()) for s in cell.replace('−', '-').split(' ± ')]
    stress = {'system': 'Klockner VEGA; 4mm, grade3 reported in article abstract', 'mean_MPa': mean, 'sd_MPa': sd, 'region': 'measured titanium surface; not proven identical to K33 critical root', 'locator': source['locator'], 'table_locator': 'materials-17-00434-t002 row V (pre-fatigue control)', 'resolution': 'PER_SURFACE_REGION', 'transport_to_K33': 'UNVALIDATED; different configuration/material/region'}
    allrows = {r['id']: r for r in csv.DictReader((ROOT / 'inputs/anchors_legacy.csv').open())}
    rows = [allrows[i] for i in legacy.S_SET]
    dd = draws(rows)
    probes = []
    for r in rows:
        d = dd[r['id']]
        (base, _) = forward(r, d)
        entries = []
        for rs in [mean - sd, mean, mean + sd]:
            (F, ports) = forward(r, d, residual_MPa=rs)
            C = np.minimum(1.43 * (d['HV'] + 120) / d['sqa'] ** (1 / 6), 1.6 * d['HV'])
            invalid = []
            for (g, s0) in ports:
                invalid.append((0.45 * g * F > C + 1e-06) & (g > 0))
            entries.append({'residual_MPa': rs, 'median_force_N': float(np.median(F)), 'ratio_to_default': float(np.median(F) / np.median(base)), 'fraction_exceeding_tensile_resistance_ceiling': float(np.mean(np.any(invalid, axis=0)))})
        probes.append({'id': r['id'], 'published_force_N': float(r['F_lim_N']), 'default_N': float(np.median(base)), 'mean_plusminus_sd_scenarios_not_confidence_interval': entries})
    inverse = json.loads((ROOT / 'raw/inverse_mechanisms.json').read_text())
    preload = []
    for x in inverse:
        if x['changed_input'] != 'preload_scale' or x['value'] is None:
            continue
        r = allrows[x['id']]
        d = dd[x['id']]
        (F, ports) = forward(r, d, preload_scale=x['value'])
        C = np.minimum(1.43 * (d['HV'] + 120) / d['sqa'] ** (1 / 6), 1.6 * d['HV'])
        invalid = np.any([(0.45 * g * F > C + 1e-06) & (g > 0) for (g, s0) in ports], axis=0)
        preload.append({'id': x['id'], 'preload_multiplier': x['value'], 'fraction_outside_tensile_fatigue_ceiling': float(np.mean(invalid)), 'causal_calibration': 'REJECTED_WITHOUT_MEASUREMENT_AND_BRANCH_VALIDITY'})
    outcome = {'claim_type': 'information_link', 'endpoint_queries': queries, 'published_XRD': stress, 'residual_stress_scenarios': probes, 'inverse_preload_branch_check': preload, 'measured_cause_of_0_608': 'UNKNOWN', 'matched_endpoint_delivered': all((x['answer']['status'] == 'PUBLISHED_ENDPOINT' for x in queries)), 'mutations': mutations, 'wall_s': time.perf_counter() - start, 'external_referent': {'kind': 'independent_measurement', 'locator': [source['locator'], 'https://doi.org/10.3390/ijerph17238988'], 'compared_quantity': 'pre-fatigue surface XRD stress and configuration-specific5e6 maximum cyclic force', 'refutes_us': True}}
    save('results_R3.json', outcome)
    print(json.dumps({k: v for (k, v) in outcome.items() if k not in ['residual_stress_scenarios', 'mutations', 'inverse_preload_branch_check', 'endpoint_queries']}))
if __name__ == '__main__':
    main()
