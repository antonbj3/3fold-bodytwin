"""Same frozen anatomical metrics, mathematically exact centroid-enclosure kernel.
Old Rtree trial/failure preserved; no sampling, threshold or prediction changed.
"""
import crownbench as C
from exact_distance import exact_distances
import pathlib, json, time

def run_eval(round_name):
    P = C.P
    validation = json.load(open(P / 'NUMERICAL_KERNEL_VERIFICATION_K2.json'))
    assert validation['all_exactness_pass']
    impl = P / f'EVALUATION_IMPLEMENTATION_{round_name}.json'
    if not impl.exists():
        C.put(impl, {'frozen_utc': C.now(), 'prediction_freeze_sha256': C.sha(P / f'FROZEN_PREDICTIONS_{round_name}.json'), 'scientific_metric_changed': False, 'sampling_changed': False, 'thresholds_changed': False, 'kernel': 'exact centroid radius enclosure; exhaustive triangle comparator1e-10mm tolerance', 'kernel_validation_sha256': C.sha(P / 'NUMERICAL_KERNEL_VERIFICATION_K2.json'), 'code_hashes': {str(p): C.sha(p) for p in [pathlib.Path(__file__).resolve(), pathlib.Path(__file__).resolve().parent / 'exact_distance.py', pathlib.Path(__file__).resolve().parent / 'evaluate.py']}, 'old_Rtree_failed_kernel_gate_retained': True})
    C.distances = exact_distances
    from evaluate import evaluate
    result = evaluate(round_name)
    result['evaluation_implementation'] = {'path': str(impl), 'sha256': C.sha(impl), 'verification': 'exact minimum agrees with exhaustive projection<=1e-10mm', 'initial_Rtree_attempt': json.load(open(P / 'ABORTED_RTREE_EVALUATION.json')) if (P / 'ABORTED_RTREE_EVALUATION.json').exists() else None, 'kernel_K1_failed_gate': json.load(open(P / 'NUMERICAL_KERNEL_FAILED_V1.json'))}
    C.put(P / f'RESULTS_{round_name}.json', result)
    return result
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--round', choices=['R1', 'R2'], default='R1')
    a = ap.parse_args()
    run_eval(a.round)
