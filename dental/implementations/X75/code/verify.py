import hashlib, json, pathlib
import numpy as np
P = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def main():
    checks = {}
    counts = {}
    for k in ['R1', 'R2', 'R3']:
        p = P / f'PREREG_{k}.json'
        checks['frozen_' + k] = sha(p) == (P / f'PREREG_{k}.sha256').read_text().strip()
    for n in ['FROZEN_PREDICTIONS', 'FROZEN_PREDICTIONS_R2', 'FROZEN_PREDICTIONS_R3']:
        checks[n] = sha(P / (n + '.json')) == (P / (n + '.sha256')).read_text().strip()
    for k in ['R1', 'R2']:
        r = json.loads((P / f'results_{k}.json').read_text())
        rows = [json.loads(l) for l in (P / f'raw/PER_SITE_{k}.jsonl').open()]
        checks['site_count_' + k] = len(rows) == 95 == r['sites']
        records = json.loads((P / f'raw/CONTROLS_{k}.json').read_text())
        checks['control_tolerances_' + k] = len(records) == 665 and all((v['pass'] and v['error'] <= (0 if v['kind'] == 'gray' else 1e-08) for v in records))
        checks['raw_hashes_' + k] = all((sha(v['raw_file']['path']) == v['raw_file']['sha256'] for v in rows))
        corevalid = 0
        quant = []
        for v in rows:
            with np.load(v['raw_file']['path']) as raw:
                g = raw['core_gray']
                corevalid += len(g) >= 30
                q = v['trabecular_candidate_roi']['gray_quantiles_05_25_50_75_95']
                quant.append(q is None if len(g) < 30 else np.array_equal(np.quantile(g, [0.05, 0.25, 0.5, 0.75, 0.95]), q))
        checks['core_quantiles_' + k] = all(quant)
        counts['valid_core_' + k] = corevalid
        pp = [v for r in rows for q in r['rays'] for v in q['profiles'].values()]
        counts['retained_profiles_' + k] = sum((v['rejection_reason'] is None for v in pp))
        checks['profile_counts_' + k] = len(pp) == 570 and counts['retained_profiles_' + k] == r['apparent_layer']['retained_contrast_profiles']
    R3 = json.loads((P / 'results_R3.json').read_text())
    checks['original_payloads_exact'] = R3['exact_shared_payload_pairs'] == 8 and all((r['max_abs_difference_gray'] == 0 and r['different_voxels'] == 0 and (r['original_payload_sha256'] == r['TF2_payload_sha256']) for r in R3['rows']))
    checks['original_fault_rejected'] = all((all(r['fault_injection'].values()) for r in R3['rows']))
    s = json.loads((P / 'SUMMARY_SUFFICIENCY.json').read_text())
    checks['identity_exact'] = s['geometry']['identity_max_error_mm'] == 0 and s['heat']['identity_error_mean'] == 0 and (s['heat']['identity_error_histogram'] == 0)
    checks['separation_and_model_enclosure'] = all(s['heat']['gates'].values())
    checks['injected_source_and_kernel_faults'] = all(s['fault_injection'].values())
    checks['physical_unknown_preserved'] = all((v['density_g_cm3'] is None and all((x is None for x in v['cortex_mm'].values())) and (v['X68']['prediction_Ncm'] is None) and (v['X67']['temperature_C'] is None) for v in map(json.loads, (P / 'SITE_INPUT_PORTS.jsonl').read_text().splitlines())))
    checks['consumer_admission_rejects'] = all((v['X68']['actual_operator_result']['status'] == 'UNKNOWN_REGIME_COMPATIBILITY' and v['X67']['actual_operator_admission'] == 'REJECTED' for v in map(json.loads, (P / 'SITE_INPUT_PORTS.jsonl').read_text().splitlines())))
    checks['own_data_budget'] = sum((pathlib.Path(v['path']).stat().st_size for k in ['R1', 'R2'] for v in json.loads((P / f'results_{k}.json').read_text())['data_files'])) < 3000000000
    out = {'pass': all(checks.values()), 'checks': checks, 'counts': counts, 'independent_scientific_review': False, 'checks_do_not_validate_anatomical_cortex_density_or_patient_response': True}
    (P / 'VERIFICATION.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))
    assert out['pass']
if __name__ == '__main__':
    main()
