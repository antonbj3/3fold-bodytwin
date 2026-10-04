"""An observed postcut-dentin event, without inventing enamel thickness."""
import csv, hashlib, json, resource, time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
ROOT = Path(__file__).resolve().parent

def event_interval(total_lower, total_upper, reduction_lower, reduction_upper, observed=False):
    if not observed:
        return {'status': 'UNKNOWN_EVENT_NOT_OBSERVED', 'dentin_interval_mm': None}
    (tl, tu, rl, ru) = (total_lower, total_upper, reduction_lower, reduction_upper)
    if not np.all(np.isfinite([tl, tu, rl, ru])) or min(tl, tu, rl, ru) < 0 or tl > tu or (rl > ru):
        raise ValueError('Invalid physical intervals')
    return {'status': 'CONDITIONAL_LOCAL_POSTCUT_DENTIN_EVENT', 'dentin_interval_mm': [max(0, tl - ru), max(0, tu - rl)], 'event_definition': 'DEJ was crossed along the actual ray, E<=actual r. A dentin endpoint is a sufficient witness; an endpoint already in pulp is a different material state', 'resolution': 'PER_POINT', 'whole_horn_coverage': 'UNKNOWN unless event and geometry cover all relevant rays', 'conditions': 'Independent same-ray event E<=actual r; total tissue/preparation error bounds must be physically established'}

def lp_extrema(tl, tu, r):
    if tu <= r:
        return (0.0, 0.0)
    bounds = [(max(tl, r), tu), (0, r), (0, None)]
    A = [[-1, 1, 0]]
    b = [0]
    low = linprog([0, 0, 1], A_ub=A, b_ub=b, A_eq=[[1, 0, -1]], b_eq=[r], bounds=bounds, method='highs')
    high = linprog([0, 0, -1], A_ub=A, b_ub=b, A_eq=[[1, 0, -1]], b_eq=[r], bounds=bounds, method='highs')
    if not low.success or not high.success:
        raise ValueError('LP event infeasible')
    return (0 if tl <= r else low.x[2], high.x[2])

def main():
    start = time.perf_counter()
    p = ROOT / 'PREREG_R4_POSTCUT_TISSUE_EVENT.json'
    if hashlib.sha256(p.read_bytes()).hexdigest() != p.with_suffix('.sha256').read_text().split()[0]:
        raise ValueError('R4 prereg drift')
    pr = json.loads(p.read_text())
    eps = pr['metrics']['digital_band_mm']
    teeth = [x for x in csv.DictReader((ROOT / 'inputs/teeth.csv').open()) if x['domain_truncated'] == 'False' and x['tooth_type'].startswith(('premolar', 'molar')) and x['occlusal_vertical_min_mm']]
    combos = sorted(set(((float(x['occlusal_vertical_min_mm']), r) for x in teeth for r in [1, 1.5, 2])))
    control = []
    for (h, r) in combos:
        tl = max(0, h - eps)
        tu = h + eps
        (lo, hi) = lp_extrema(tl, tu, r)
        expected = event_interval(tl, tu, r, r, True)['dentin_interval_mm']
        err = max(abs(lo - expected[0]), abs(hi - expected[1]))
        wrong = max(abs(lo - (expected[0] + 1)), abs(hi - expected[1]))
        control.append({'H_mm': h, 'r_mm': r, 'LP_min_D_mm': lo, 'LP_max_D_mm': hi, 'error_mm': err, 'valid_pass': err <= 1e-07, 'injected_1mm_lower_bound_rejected': wrong > 1e-07})
    if not all((x['valid_pass'] and x['injected_1mm_lower_bound_rejected'] for x in control)):
        raise AssertionError('Event reduction failed')
    absent = event_interval(3, 3.2, 2, 2, False)
    if absent['status'] != 'UNKNOWN_EVENT_NOT_OBSERVED':
        raise AssertionError('Unobserved event promoted')
    wrong_world = {'T_mm': 4, 'E_mm': 2.5, 'r_mm': 1, 'correct_D_mm': 1.5, 'naive_T_minus_r_mm': 3, 'event_E_le_r': False, 'naive_equality_refuted': True, 'kind': 'our_own_fixture'}
    summary = []
    per = []
    for r in [1, 1.5, 2]:
        for threshold in [0.5, 1]:
            groups = []
            for t in teeth:
                h = float(t['occlusal_vertical_min_mm'])
                (lo, hi) = event_interval(max(0, h - eps), h + eps, r, r, True)['dentin_interval_mm']
                cls = 'ABOVE_ALL' if lo >= threshold else 'BELOW_ALL' if hi < threshold else 'UNKNOWN_INTERVAL'
                groups.append(cls)
                per.append({'case': t['case'], 'fdi': t['fdi'], 'reduction_mm': r, 'threshold_mm': threshold, 'counterfactual_event_dentin_lower_mm': lo, 'counterfactual_event_dentin_upper_mm': hi, 'counterfactual_class': cls, 'event_actually_observed': False, 'actual_dentin_class': 'UNKNOWN', 'resolution': 'PER_TOOTH proxy minimum; requires PER_POINT coverage of all horn rays'})
            summary.append({'reduction_mm': r, 'threshold_mm': threshold, 'teeth_n': len(teeth), 'counterfactual_event_above_per100': 100 * groups.count('ABOVE_ALL') / len(groups), 'counterfactual_event_below_per100': 100 * groups.count('BELOW_ALL') / len(groups), 'counterfactual_event_unresolved_per100': 100 * groups.count('UNKNOWN_INTERVAL') / len(groups), 'observed_events_n': 0, 'physical_dentin_answer': 'UNKNOWN', 'resolution': 'POPULATION'})
    for (name, rows) in [('raw/R4_EVENT_PER_TOOTH.csv', per), ('raw/R4_EVENT_LP.csv', control), ('tables/R4_MINIMUM_OBSERVATION.csv', summary)]:
        with (ROOT / name).open('w') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    result = {'round': 'R4', 'claim_type': 'capability', 'verdict': 'EVENT_SUFFICIENT_UNDER_COLUMN_ASSUMPTIONS', 'external_referent': pr['external_referent'], 'LP_combinations_n': len(combos), 'all_LP_valid_pass': True, 'all_LP_faults_rejected': True, 'absent_event_status': absent['status'], 'counterexample_without_event': wrong_world, 'counterfactual_summaries': summary, 'observed_events_n': 0, 'physical_dentin': 'UNKNOWN', 'wall_seconds': time.perf_counter() - start, 'maxrss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024}
    result['event_semantics'] = 'Formal condition E<=r means the cut crossed DEJ. Observing final surface dentin is one sufficient witness. A cut reaching pulp must be recorded as pulp exposure, not a dentin surface. Counterfactual rows apply to the weaker DEJ-crossed event; no endpoint material was observed.'
    (ROOT / 'raw/R4_EVENT_RESULTS.json').write_text(json.dumps(result, indent=2) + '\n')
    (ROOT / 'HANDOFF_R4.md').write_text("# R4 : minimum local material observation\n\nAn independent observation that the current cut surface of the same horn beam is dentin (E ≤ r) is sufficient for the residual hard tissue to be the residue of the layered column model. Actual HiGHS -extremization tests the limits; a 1mm error falls in each query. Without observation, T − r can give dentin error; the counter example can be preserved. No event was observed in the dataset, and the whole horn geometry/measurement error still had to be covered. The next design is a blind same-specimen-specimen with frozen material event and measured total gap/preparation, rather than full DEJ segmentation.\n")
    print('R4 postcut event sufficient under model; ' + str(len(combos)) + ' LP combinations PASS; 0 events measured.')
if __name__ == '__main__':
    main()
