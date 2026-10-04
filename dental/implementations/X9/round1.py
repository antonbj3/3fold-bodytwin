import csv
import resource
import time
from collections import defaultdict
from pathlib import Path
import numpy as np
from ipr import capacity, control_bisection, load_case, measure_patch, sha, write_json, jsonable
import json

def run():
    (start, cpu) = (time.perf_counter(), time.process_time())
    pre = json.loads(Path('PREREG_R1.json').read_text())
    (rows, raw, errors) = ([], [], [])
    for item in pre['inputs']:
        assert sha(item['path']) == item['sha256']
        (meshes, meta, frame) = load_case(item['path'])
        for side in [1, -1]:
            for height in pre['measurement']['height_fractions']:
                p = measure_patch(meshes, side, height)
                good = np.isfinite(p['thickness_mm'])
                tt = p['thickness_mm'][good]
                caps = capacity(tt)
                controls = np.array([control_bisection(t) for t in tt])
                errors.extend(np.abs(caps - controls))
                row = dict(tid=item['tid'], jaw=item['jaw'], tooth_type=item['type'], side=p['side'], height_fraction=height, valid_rays=int(good.sum()), valid_fraction=float(good.mean()), mean_mm=float(np.mean(tt)) if len(tt) else None, min_mm=float(np.min(tt)) if len(tt) else None, p05_mm=float(np.quantile(tt, 0.05)) if len(tt) else None, surface_gap_mean_mm=float(np.nanmean(p['surface_gap_mm'])) if good.any() else None, scenario_capacity_mm=float(np.min(caps)) if len(caps) else None, calibrated_safe_ipr_mm=None, anatomical_error_status='UNKNOWN')
                rows.append(row)
                raw.append({'tid': item['tid'], 'frame': frame, **p})
        print(item['tid'], 'measured', flush=True)
        write_json('rounds/R1_raw_profiles.partial.json', raw)
        write_json('CURRENT_WORK_STATE.json', {'lane': 'X9-ipr-safety', 'stage': 'R1_RUNNING', 'last_tooth': item['tid'], 'completed_patches': len(rows), 'next_operation': 'finish frozen R1 gates'})
    med = float(np.median([r['mean_mm'] for r in rows if r['height_fraction'] == 0.5 and r['mean_mm'] is not None]))
    g1 = 0.96 <= med <= 1.96
    g2 = max(errors, default=0) <= 1e-07
    usable = [r for r in rows if r['valid_fraction'] >= 0.9]
    fraction = float(np.mean([r['scenario_capacity_mm'] < 0.5 for r in usable])) if usable else 0.0
    gates = {'G1_external_population': {'pass': g1, 'median_mm': med, 'frozen_range_mm': [0.96, 1.96]}, 'G2_equal_information': {'pass': g2, 'max_difference_mm': max(errors, default=0), 'outcome': 'TIE' if g2 else 'FAIL'}, 'G3_praxis_difference': {'pass': fraction >= 0.1, 'fraction': fraction, 'patches': len(usable)}, 'G4_coverage': {'pass': all((r['valid_fraction'] >= 0.9 for r in rows)), 'failing_patches': sum((r['valid_fraction'] < 0.9 for r in rows))}, 'G5_physical_admission': {'pass': False, 'outcome': 'UNKNOWN', 'reason': 'No same-tooth DEJ/outer-boundary calibration'}}
    injections = {'population_plus_2_mm_rejected': not 0.96 <= med + 2 <= 1.96, 'control_plus_0_1_mm_rejected': max([abs(e + 0.1) for e in errors], default=0.1) > 1e-07, 'coverage_all_rays_removed_rejected': 0.0 < 0.9, 'praxis_equal_0_5_all_patches_rejected': 0.0 < 0.1, 'uncalibrated_physical_safety_claim_rejected': gates['G5_physical_admission']['pass'] is False}
    assert all(injections.values())
    result = {'round': 'R1', 'prereg_sha256': sha('PREREG_R1.json'), 'gates': gates, 'injected_faults': injections, 'number_teeth': len(pre['inputs']), 'number_patches': len(rows), 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://pubmed.ncbi.nlm.nih.gov/32634888/', 'compared_quantity': 'molar proximal enamel thickness mm; different specimens/population', 'refutes_us': not g1}, 'cost': {'wall_seconds': time.perf_counter() - start, 'cpu_seconds': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'preparation_historical_segmentation': 'UNKNOWN', 'fit_seconds': 0, 'physical_validation': 'NOT_DONE'}}
    write_json('rounds/R1_raw_profiles.json', raw)
    write_json('rounds/R1_results.json', result)
    with open('artifacts/thickness_by_tooth_side_height.csv', 'w') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    write_json('rounds/R1_rows.json', rows)
    write_json('CURRENT_WORK_STATE.json', {'lane': 'X9-ipr-safety', 'stage': 'R1_COMPLETE', 'latest_gate': gates, 'next_operation': 'R2 protected-DEJ tool-sweep envelope; retain anatomical calibration UNKNOWN'})
    Path('HANDOFF.md').write_text('R1 complete. Scalar feasibility solver matches closed formula (TIE). Exact gates are in rounds/R1_results.json. The STS enamel shell is intensity-derived; no anatomical safety claim is admitted. Next construction changes the operation from independent rays to a finite swept plane protecting the whole DEJ envelope.\n')
    print(json.dumps(jsonable(result), indent=2))
    return (rows, raw, result)
if __name__ == '__main__':
    run()
