import json, hashlib, time, itertools
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.optimize import OptimizeWarning
import warnings
warnings.filterwarnings('ignore', category=OptimizeWarning, message='Unrecognized options detected')
from calibration_port import matrix
R = Path(__file__).resolve().parents[1]

def put(n, x):
    (R / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def envelope(rows, H, compression_error, wrench_error):
    lo = []
    hi = []
    y = []
    scales = np.array([1, 1, 1, 1 / H, 1 / H, 1 / H])
    for row in rows:
        (_, compression) = matrix(row)
        sig = -np.asarray(row['region_signatures']).T
        small = sig * np.maximum(compression - compression_error, 0)
        large = sig * np.maximum(compression + compression_error, 0)
        lo.append(np.minimum(small, large) * scales[:, None])
        hi.append(np.maximum(small, large) * scales[:, None])
        y.append(np.r_[row['force_N'], row['moment_Nmm']] * scales)
    low = np.vstack(lo)
    high = np.vstack(hi)
    target = np.concatenate(y)
    upper = np.vstack([low, -high])
    rhs = np.r_[target + wrench_error, -target + wrench_error]
    return (upper, rhs)

def bounds(A, b):
    low = []
    high = []
    statuses = []
    for i in range(5):
        c = np.eye(5)[i]
        a = linprog(c, A_ub=A, b_ub=b, bounds=[(0, None)] * 5, method='highs', options={'threads': 4, 'parallel': False})
        z = linprog(-c, A_ub=A, b_ub=b, bounds=[(0, None)] * 5, method='highs', options={'threads': 4, 'parallel': False})
        statuses.append([a.status, z.status])
        if a.status == 2 or z.status == 2:
            return dict(status='INFEASIBLE', LP_statuses=statuses, low=None, high=None)
        low.append(float(a.fun) if a.success else None)
        high.append(float(-z.fun) if z.success else None)
    return dict(status='BOUNDED' if all((v is not None for v in low + high)) else 'UNBOUNDED', low=low, high=high, LP_statuses=statuses)

def main():
    start = time.perf_counter()
    pr = json.loads((R / 'PREREG_R4.json').read_text())
    ge = pr['gates']
    assert hashlib.sha256((R / 'PREREG_R4.json').read_bytes()).hexdigest() == json.loads((R / 'PREREG_R4.sha256.json').read_text())['sha256']
    observations = json.loads((R / 'raw/R3_FIXTURE_OBSERVATIONS.json').read_text())
    r3 = json.loads((R / 'round3/results.json').read_text())
    predictions = json.loads((R / 'FROZEN_PREDICTIONS_R3.json').read_text())['predictions']
    rows = []
    injections = []
    for case in predictions:
        obs = [o for o in observations if o['specimen_id'] == 'OUR_OWN_FIXTURE_' + case['case']]
        design = next((r for r in r3['rows'] if r['case'] == case['case']))
        true = np.array(case['fixture_stiffness_N_per_mm'])
        selected = [o for o in obs if o['probe_id'] in design['selected_probe_ids']]
        for (label, group) in [('R3 selected3', selected), ('all16 frozen probes', obs)]:
            for (compression, wrench) in itertools.product(pr['error_scenarios']['compression_bound_mm'], pr['error_scenarios']['wrench_bound_after_moment_scaling_N']):
                (A, b) = envelope(group, case['H_mm'], compression, wrench)
                bound = bounds(A, b)
                truth_ok = bool(np.all(A @ true - b <= ge['LP_feasibility_tolerance']))
                relative = None
                if bound['status'] == 'BOUNDED':
                    relative = float(np.max(np.maximum((true - np.array(bound['low'])) / true, (np.array(bound['high']) - true) / true)))
                passed = relative is not None and relative <= ge['max_regional_relative_bound']
                rows.append(dict(case=case['case'], queries=label, compression_error_mm=compression, scaled_wrench_error_N=wrench, bounds=bound, max_regional_relative_bound=relative, gate=passed, frozen_true_k_feasible=truth_ok, source_kind='our_own_fixture', resolution_level='PER_SURFACE_REGION'))
            corrupt = [dict(o) for o in group]
            corrupt[0]['force_N'] = list(corrupt[0]['force_N'])
            corrupt[0]['force_N'][0] += 10
            (A, b) = envelope(corrupt, case['H_mm'], 0, 0)
            injected = bounds(A, b)
            injections.append(dict(case=case['case'], queries=label, injected_status=injected['status'], rejected=injected['status'] == 'INFEASIBLE'))
    requirements = []
    for label in ['R3 selected3', 'all16 frozen probes']:
        for wrench in pr['error_scenarios']['wrench_bound_after_moment_scaling_N']:
            passing = [e for e in pr['error_scenarios']['compression_bound_mm'] if all((r['gate'] for r in rows if r['queries'] == label and r['scaled_wrench_error_N'] == wrench and (r['compression_error_mm'] == e)))]
            requirements.append(dict(queries=label, scaled_wrench_error_N=wrench, largest_passing_compression_bound_mm=max(passing) if passing else None, meaning='Across all five effective stiffnesses and all three own fixtures on the frozen error grid; not a proven hardware tolerance'))
    zero = [r for r in rows if r['compression_error_mm'] == 0 and r['scaled_wrench_error_N'] == 0]
    result = dict(round='R4', claim_type='capability', outcome='REGIONAL_CALIBRATION_METROLOGY_REQUIREMENT_COMPUTED_CONDITIONALLY', external_referent=pr['external_referent'], measurement_resolution_referent=pr['measurement_resolution_referent'], fixture_referent=pr['fixture_referent'], requirements=requirements, rows=rows, injections=injections, gates=dict(frozen_k_within_all_envelopes=all((r['frozen_true_k_feasible'] for r in rows)), zero_error_recovery=all((r['max_regional_relative_bound'] is not None and r['max_regional_relative_bound'] <= ge['exact_zero_error_recovery_relative'] for r in zero)), injected_wrong_wrench_rejected=all((r['rejected'] for r in injections))), scenarios=len(rows), wall_s=time.perf_counter() - start, physical_validation='UNKNOWN; no matched metrology/force measurement', next_operation='Matched rigid-typodont metrology/force/DIC experiment; use source-conditioned probe set and independently test the regional model rather than assume voxel spacing is its error bound')
    put('round4/results.json', result)
    put('CURRENT_WORK_STATE.json', dict(lane='X27-aligner-regional', state='R4_COMPLETE', latest_gate=result['outcome'], next_operation='Fit single central-incisor gain with two unfit tooth-type means'))
    (R / 'round4/HANDOFF.md').write_text('# R4 metrology handoff\n\nRegional rank alone does not close calibration. Conditional LP outer bounds computed for all five stiffnesses and three fixtures. Voxel size is not an accuracy guarantee.\n\n' + json.dumps(requirements, indent=2) + '\n\n' + result['next_operation'] + '\n')
    print(json.dumps(dict(outcome=result['outcome'], requirements=requirements, gates=result['gates']), indent=2))
if __name__ == '__main__':
    main()
