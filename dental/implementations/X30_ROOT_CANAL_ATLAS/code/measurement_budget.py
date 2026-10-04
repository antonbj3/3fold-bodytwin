import json, time, resource, numpy as np
from operators import circumcircle
from uncertainty import radius_interval
from atlas import ROOT, DATA, dump, csvwrite

def precision_budget(tr, r, tol=0.2):

    def passes(eps):
        (lo, up) = radius_interval(tr, eps)
        return up is not None and lo >= (1 - tol) * r and (up <= (1 + tol) * r)
    (a, b) = (0.0, 1.0)
    for _ in range(40):
        m = (a + b) / 2
        if passes(m):
            a = m
        else:
            b = m
    return (a, passes(a), not passes(max(a * 1000, 0.1)))

def main():
    start = time.monotonic()
    cpu = time.process_time()
    rows = []
    for row in map(json.loads, (DATA / 'R1_paths.jsonl').open()):
        p = np.array(row['smooth_points_mm'])
        for window in [3.0, 6.0, 9.0]:
            n = round(window / 0.3 / 2)
            candidates = []
            for i in range(n, len(p) - n):
                tr = p[[i - n, i, i + n]]
                (r, angle) = circumcircle(*tr)
                if np.isfinite(r):
                    candidates.append((r, angle, tr))
            if not candidates:
                continue
            (r, angle, tr) = min(candidates, key=lambda q: q[0])
            (eps, ok, poison) = precision_budget(tr, r)
            rows.append({'case': row['case'], 'fdi': row['fdi'], 'branch_id': row['branch_id'], 'resolution_level': 'PER_SURFACE_REGION', 'window_axial_mm': window, 'nominal_radius_mm': r, 'window_turn_deg': angle, 'required_point_precision_mm': eps, 'compatible_with_voxel_half_diagonal': eps >= np.sqrt(3) * 0.3 / 2, 'budget_gate_pass': ok, 'poison_precision_x1000_rejected': poison, 'bending_strain_per_diameter_per_mm': 1 / (2 * r), 'fatigue_cycles_to_failure': None, 'clinical_fracture_probability': None})
    csvwrite(ROOT / 'raw/R3_measurement_budget.csv', rows)
    summary = []
    for window in [3.0, 6.0, 9.0]:
        rs = [r for r in rows if r['window_axial_mm'] == window]
        eps = np.array([r['required_point_precision_mm'] for r in rs])
        frac = float(np.mean([r['compatible_with_voxel_half_diagonal'] for r in rs])) if rs else 0
        summary.append({'window_axial_mm': window, 'n_branches': len(rs), 'required_point_precision_P5_P50_P95_mm': np.quantile(eps, [0.05, 0.5, 0.95]).tolist() if rs else None, 'compatible_fraction': frac, 'frozen_80percent_gate_pass': frac >= 0.8, 'resolution_level': 'POPULATION'})
    res = {'claim_type': 'capability', 'summary': summary, 'all_budget_gate_pass': all((r['budget_gate_pass'] for r in rows)), 'all_precision_poison_rejected': all((r['poison_precision_x1000_rejected'] for r in rows)), 'wall_s': time.monotonic() - start, 'cpu_s': time.process_time() - cpu, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fatigue_link': {'claim_type': 'information_link', 'status': 'QUALITATIVE_SOURCE_LINK_ONLY', 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://doi.org/10.1016/S0099-2399(97)80250-6', 'compared_quantity': 'Cycles-to-failure decreases for radius5->2mm at fixed instrument;angle alone insufficient', 'refutes_us': False}, 'absolute_cycles': 'UNKNOWN', 'absolute_probability': 'UNKNOWN', 'patient_risk': 'UNKNOWN'}}
    dump(ROOT / 'raw/R3_result.json', res)
    print(json.dumps(res), flush=True)
if __name__ == '__main__':
    main()
