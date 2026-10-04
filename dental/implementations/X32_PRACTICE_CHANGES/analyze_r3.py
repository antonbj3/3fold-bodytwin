import json, csv, hashlib, copy
from pathlib import Path
R = Path(__file__).resolve().parent

def gate(original, observed):
    keys = [(x['case'], str(x['fdi'])) for x in original]
    assert len(keys) == len(set(keys)), 'duplicate source key'
    assert set(observed) <= set(keys), 'observed source key absent'
    return set(keys)

def main():
    binding = json.loads((R / 'raw/R3_INPUT_BINDING.json').read_text())
    p = Path(binding['path'])
    assert hashlib.sha256(p.read_bytes()).hexdigest() == binding['sha256'], 'R3 source drift'
    source = json.loads(p.read_text())
    rows = source['sites']
    obs = {(x['patient'], x['object']): x['changed'] == 'True' for x in csv.DictReader((R / 'raw/decision_ledger.csv').open()) if x['demo'] == 'X8'}
    keys = gate(rows, obs)
    missing = keys - set(obs)
    changes = {k for (k, v) in obs.items() if v}
    patients = {k[0] for k in keys}
    already = {k[0] for k in changes}
    uncertain = {k[0] for k in missing}
    faults = {}
    for (name, data) in [('duplicate_original_key', rows + [rows[0]]), ('omitted_observed_key', [x for x in rows if (x['case'], str(x['fdi'])) != next(iter(obs))])]:
        try:
            gate(data, obs)
        except AssertionError:
            faults[name + '_rejected'] = True
    lower = sum(obs.values())
    upper = lower + len(missing)
    assert lower == sum({**obs, **{k: False for k in missing}}.values())
    assert upper == sum({**obs, **{k: True for k in missing}}.values())
    try:
        assert upper <= lower - 1
    except AssertionError:
        faults['inverted_bound_rejected'] = True
    incomplete = [x for x in rows if (x['case'], str(x['fdi'])) in missing]
    out = {'round': 'R3', 'claim_type': 'information_link', 'source_binding': binding, 'source_declared_n_cases': source['n_cases'], 'represented_original_patient_keys': len(patients), 'represented_original_sites': len(keys), 'complete_sites': len(obs), 'unmatched_original_sites': len(missing), 'unmatched_rows_with_no_L_plan': sum((x.get('L_plan') is None for x in incomplete)), 'unmatched_rows_with_no_defined_plan_note': 'A binary decision on these rows is a hypothetical completion, not an observed plan comparison', 'already_changed_patients': len(already), 'patients_with_unmatched_sites': len(uncertain), 'additional_potentially_changed_patients': len(uncertain - already), 'site_completion_bound_per100': [100 * lower / len(keys), 100 * upper / len(keys)], 'patient_any_change_completion_bound_per100': [100 * len(already) / len(patients), 100 * len(already | uncertain) / len(patients)], 'bound_type': 'Logical binary-completion identification bound, not a sampling confidence interval', 'original_cohort_rate': 'UNKNOWN when an excluded plan is undefined; no complete-case generalization', 'same_information_control': 'Explicit extreme binary completions match set bounds', 'fault_injections': faults, 'resolution_levels': {'sites': 'PER_TOOTH', 'patient_aggregation': 'PER_ARCH', 'reported_cohort_rates': 'POPULATION'}, 'timescale': 'SIMULTANEOUS'}
    (R / 'rounds/RESULTS_R3.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    main()
