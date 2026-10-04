from crownbench import *

def verify(write_receipt=True):
    split = json.load(open(P / 'SPLIT.json'))
    sets = [set(split['patients'][k]) for k in ('train', 'test_R1', 'reserve_R2')]
    assert not (sets[0] & sets[1] or sets[0] & sets[2] or sets[1] & sets[2])
    checks = []
    for rn in ['R1', 'R2']:
        pr = P / f'PREREG_{rn}.json'
        frp = P / f'FROZEN_PREDICTIONS_{rn}.json'
        fr = json.load(open(frp))
        res = json.load(open(P / f'RESULTS_{rn}.json'))
        assert sha(pr) == fr['prereg_sha256']
        assert sha(P / 'SPLIT.json') == fr['split_sha256']
        assert sha(frp) == res['frozen_predictions_sha256']
        assert sha(P / f'RAW_METRICS_{rn}.jsonl') == res['raw_metrics_sha256']
        assert sha(fr['predictions_index']) == fr['index_sha256']
        for d in fr['prediction_files']:
            assert sha(d['path']) == d['sha256']
        assert fr['frozen_utc'] < res['evaluated_utc']
        group = 'test_R1' if rn == 'R1' else 'reserve_R2'
        contexts = json.load(open(DATA / 'public' / f'{group}_INDEX.json'))
        for c in contexts:
            if c['eligible']:
                with np.load(c['context']) as d:
                    assert int(d['target_fdi']) not in np.unique(d['labels'])
                assert sha(c['context']) == c['context_sha256']
        controls = json.load(open(P / f'METRIC_CONTROLS_{rn}.json'))
        assert all((not a['gate_pass'] for a in controls['injected_threshold_values'].values()))
        assert controls['all_pass'] and all(controls['checks'].values()), 'GEOMETRY_SENSITIVITY_CHECK_FAILED'
        assert res['metric_controls_all_pass'], 'RESULT_METRIC_CHECK_FAILED'
        checks.append({'round': rn, 'prediction_hashes': len(fr['prediction_files']), 'scheduled_contexts': len(contexts), 'generation_failures': fr['failure_count'], 'geometry_injection_checks_all_pass': controls['all_pass'], 'frozen_before_evaluation': True})
    if (P / 'RESULTS_R3.json').exists():
        rr = json.load(open(P / 'RESULTS_R3.json'))
        fr = json.load(open(P / 'FROZEN_PREDICTIONS_R3.json'))
        assert sha(P / 'PREREG_R3.json') == fr['prereg_sha256']
        assert sha(P / 'SPLIT_R3.json') == fr['split_sha256']
        assert sha(fr['predictions_index']) == fr['index_sha256']
        assert sha(P / 'FROZEN_PREDICTIONS_R3.json') == rr['frozen_predictions_sha256']
        assert sha(P / 'RAW_METRICS_R3.jsonl') == rr['raw_metrics_sha256']
        assert fr['frozen_utc'] < rr['evaluated_utc']
        for r in fr['prediction_files']:
            assert sha(r['path']) == r['sha256']
        new = json.load(open(P / 'SPLIT_R3.json'))
        assert not set(new['cases']) & set(sum(split['patients'].values(), []))
        controls = json.load(open(P / 'METRIC_CONTROLS_R3.json'))
        assert controls['all_pass'] and all(controls['checks'].values()), 'R3_GEOMETRY_SENSITIVITY_CHECK_FAILED'
        assert all((not a['gate_pass'] for a in controls['injected_threshold_values'].values())), 'R3_THRESHOLD_CHECK_FAILED'
    assert (P / 'benchmark_figure.png').stat().st_size > 1000
    assert (P / 'README_DEMO.md').exists()
    assert disk_guard() < 3000000000
    if write_receipt:
        put(P / 'VERIFICATION.json', {'verified_utc': now(), 'patient_disjoint': True, 'target_vertex_exclusion': True, 'hashes_valid': True, 'rounds': checks, 'data_bytes': disk_guard(), 'scientific_review': 'PENDING_INDEPENDENT_REVIEW', 'geometry_gate_failures_retained': True})
    print('VERIFIED hashes, split, target exclusion, freeze ordering, injected-value rejection and disk budget', flush=True)
if __name__ == '__main__':
    verify()
