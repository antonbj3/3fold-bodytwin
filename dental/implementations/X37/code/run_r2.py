import json, csv, time, itertools, resource, datetime
import numpy as np
from scipy.spatial.transform import Rotation
from geometry import ROOT, apply, mat, angle
start = time.monotonic()
rows = json.loads((ROOT / 'raw/R1_POSES.json').read_text())
ix = {(r['patient'], r['jaw'], r['source'], r['week'], r['fdi']): r for r in rows}
pairs = []
checks = []
rng = np.random.default_rng(3702)
Q = mat(Rotation.random(random_state=rng).as_matrix(), np.array([100, -50, 25]))

def pose(p, j, s, t, f, ref='bilateral_molars'):
    return np.array(ix[p, j, s, t, f]['H'][ref])

def center(p, j, s, t, f, ref='bilateral_molars'):
    return apply(pose(p, j, s, t, f, ref), np.array(ix[p, j, s, t, f]['centroid_base_mm']))

def distance(p, j, s, t, i, k, ref='bilateral_molars'):
    return float(np.linalg.norm(center(p, j, s, t, i, ref) - center(p, j, s, t, k, ref)))
max_dist_err = 0.0
max_mat_err = 0.0
rejected = 0
for p in ['3485', '6457']:
    for j in ['OK', 'UK']:
        ids = sorted({r['fdi'] for r in rows if r['patient'] == p and r['jaw'] == j})
        for (i, k) in itertools.combinations(ids, 2):
            for t in range(1, 10):
                ds = distance(p, j, 'Sirona', t, i, k)
                ds0 = distance(p, j, 'Sirona', 0, i, k)
                dsp = distance(p, j, 'Sirona', t - 1, i, k)
                dp = distance(p, j, 'Bottmedical', t, i, k)
                dp0 = distance(p, j, 'Bottmedical', 0, i, k)
                dpp = distance(p, j, 'Bottmedical', t - 1, i, k)
                Hsi = pose(p, j, 'Sirona', t, i)
                Hsk = pose(p, j, 'Sirona', t, k)
                Hi0 = pose(p, j, 'Sirona', t - 1, i)
                Hk0 = pose(p, j, 'Sirona', t - 1, k)
                relative = np.linalg.inv(Hsk) @ Hsi
                relative_previous = np.linalg.inv(Hk0) @ Hi0
                error = 0.0
                for src in ['Sirona', 'Bottmedical']:
                    for tt in [t, 0]:
                        d = distance(p, j, src, tt, i, k)
                        refR = np.linalg.inv(pose(p, j, src, tt, k)) @ pose(p, j, src, tt, i)
                        ci = np.array(ix[p, j, src, tt, i]['centroid_base_mm'])
                        ck = np.array(ix[p, j, src, tt, k]['centroid_base_mm'])
                        direct = float(np.linalg.norm(apply(refR, ci) - ck))
                        error = max(error, abs(d - direct))
                        for ref in ['all_crowns', 'left_molars', 'right_molars']:
                            dr = distance(p, j, src, tt, i, k, ref)
                            rR = np.linalg.inv(pose(p, j, src, tt, k, ref)) @ pose(p, j, src, tt, i, ref)
                            max_dist_err = max(max_dist_err, abs(d - dr))
                            max_mat_err = max(max_mat_err, float(np.max(np.abs(refR - rR))))
                        ai = Q @ pose(p, j, src, tt, i)
                        ak = Q @ pose(p, j, src, tt, k)
                        error = max(error, abs(d - np.linalg.norm(apply(ai, ci) - apply(ak, ck))))
                        max_mat_err = max(max_mat_err, float(np.max(np.abs(refR - np.linalg.inv(ak) @ ai))))
                max_dist_err = max(max_dist_err, error)
                req = [ix[p, j, src, tt, f] for src in ['Sirona', 'Bottmedical'] for tt in [t, 0] for f in [i, k]]
                rmove = [ix[p, j, 'Sirona', tt, f] for tt in [t, t - 1] for f in [i, k]]
                u = 0.8 + sum((x['sampling_spread_mm'] for x in req))
                uv = 0.4 + sum((x['sampling_spread_mm'] for x in rmove))
                good = all((x['surface_fit_pass'] for x in req + rmove))
                if not good:
                    rejected += 1
                wrong_reject = abs(ds + 1 - np.linalg.norm(apply(relative, np.array(ix[p, j, 'Sirona', t, i]['centroid_base_mm'])) - np.array(ix[p, j, 'Sirona', t, k]['centroid_base_mm']))) > 1e-07
                pairs.append({'patient': p, 'jaw': j, 'fdi_i': i, 'fdi_j': k, 'week': t, 'resolution': 'PER_TOOTH', 'observed_separation_mm': ds, 'observed_weekly_separation_change_mm': ds - dsp, 'planned_weekly_separation_change_mm': dp - dpp, 'observed_cumulative_separation_change_mm': ds - ds0, 'planned_cumulative_separation_change_mm': dp - dp0, 'relative_tracking_error_mm': abs(ds - ds0 - (dp - dp0)), 'movement_sensitivity_budget_mm': uv, 'tracking_sensitivity_budget_mm': u, 'budget_kind': 'PHENOMENOLOGICAL worst-sum floor + operational sampling, not clinical CI', 'observed_relative_weekly_rotation_deg': angle(relative @ np.linalg.inv(relative_previous)), 'relative_H': relative.tolist(), 'surface_fit_pass': good, 'relative_weekly_motion_resolved': bool(good and abs(ds - dsp) > uv), 'relative_tracking_gap_resolved': bool(good and abs(ds - ds0 - (dp - dp0)) > u), 'individual_injected_1mm_rejected': bool(wrong_reject), 'absolute_motion': 'UNKNOWN', 'independent_population_units': 2})
result = {'round': 'R2', 'claim_type': 'information_link', 'external_referent': json.loads((ROOT / 'PREREG_R2.json').read_text())['external_referent'], 'pair_week_rows': len(pairs), 'patient_count': 2, 'resolution': 'PER_TOOTH', 'surface_fit_rejected': rejected, 'max_reference_distance_error_mm': max_dist_err, 'max_reference_relative_matrix_error': max_mat_err, 'gauge_gate_pass': max_dist_err <= 1e-07 and max_mat_err <= 1e-09, 'injected_wrong_scalar_rejected_all': all((r['individual_injected_1mm_rejected'] for r in pairs)), 'resolved_weekly_motion': sum((r['relative_weekly_motion_resolved'] for r in pairs)), 'resolved_tracking_gap': sum((r['relative_tracking_gap_resolved'] for r in pairs)), 'elapsed_s': time.monotonic() - start, 'peak_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'by_patient': {p: {'rows': sum((r['patient'] == p for r in pairs)), 'resolved_tracking_rows': sum((r['patient'] == p and r['relative_tracking_gap_resolved'] for r in pairs)), 'max_tracking_mm': max((r['relative_tracking_error_mm'] for r in pairs if r['patient'] == p))} for p in ['3485', '6457']}}
for (n, d) in [('R2_PAIRS.json', pairs), ('R2_RESULTS.json', result)]:
    (ROOT / 'raw' / n).write_text(json.dumps(d, indent=2) + '\n')
flat = [{k: v for (k, v) in r.items() if k != 'relative_H'} for r in pairs]
with (ROOT / 'raw' / 'relative_movement.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(flat[0]))
    w.writeheader()
    w.writerows(flat)
print(json.dumps(result, indent=2))
state = json.loads((ROOT / 'CURRENT_WORK_STATE.json').read_text())
state.update(status='R2_MEASURED', latest_gate='Gauge-invariant relative observable computed; absolute motion remains unknown', next_operation='R3 rolling stage geometry -> X27 conditional wrench -> observed weekly outcome, plus power assumptions', updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2) + '\n')
