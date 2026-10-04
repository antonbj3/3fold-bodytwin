"""R5: preserve the normal-ray cut component at the target, not a global root.

The old first-root reduction can cross retained material to a cutting ball
on the opposite side of a wall. This local union remains a sparse witness
set and is not a complete CAM remaining-stock model.
"""
from dental_release.paths import expand as _release_expand
import time, math
from fractions import Fraction as Q
from common import *
from collision import first_union_ball_entry
from run_r1 import lib

def local_entry(points, normals, centres, radii, window=(-0.4, 0.8)):
    if not len(centres):
        return np.full(len(points), np.nan)
    v = np.asarray(centres)[None, :, :] - np.asarray(points)[:, None, :]
    a = np.sum(v * np.asarray(normals)[:, None, :], axis=2)
    rho = np.asarray(radii)[None, :] ** 2 - (np.sum(v * v, axis=2) - a * a)
    root = np.sqrt(np.maximum(rho, 0))
    left = a - root
    right = a + root
    good = (rho >= 0) & (right >= window[0]) & (left <= window[1])
    out = np.full(len(points), np.nan)
    for i in range(len(points)):
        intervals = sorted(((max(float(left[i, j]), window[0]), min(float(right[i, j]), window[1])) for j in np.flatnonzero(good[i])))
        merged = []
        for (lo, hi) in intervals:
            if merged and lo <= merged[-1][1]:
                merged[-1][1] = max(hi, merged[-1][1])
            else:
                merged.append([lo, hi])
        relevant = next((x for x in merged if x[1] >= 0), None)
        if relevant is not None:
            out[i] = relevant[0]
    return out

def exact_coaxial_oracle(centres, radii):
    intervals = sorted(((Q(str(c)) - Q(str(r)), Q(str(c)) + Q(str(r))) for (c, r) in zip(centres, radii)))
    merged = []
    for (lo, hi) in intervals:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(hi, merged[-1][1])
        else:
            merged.append([lo, hi])
    relevant = next((x for x in merged if x[1] >= 0), None)
    return None if relevant is None else float(relevant[0])

def controls():
    p = np.zeros((1, 3))
    n = np.array([[0.0, 0.0, 1.0]])
    configs = [([-0.3125, 0.1875], [0.0625, 0.0625]), ([-0.0625], [0.3125])]
    rows = []
    for (cs, rs) in configs:
        centres = np.c_[np.zeros(len(cs)), np.zeros(len(cs)), cs]
        old = first_union_ball_entry(p, n, centres, rs)[0]
        new = local_entry(p, n, centres, rs)[0]
        exact = exact_coaxial_oracle(cs, rs)
        assert new == exact and abs(new + 0.05 - exact) > 1e-06
        rows.append(dict(global_first_entry_mm=old, local_component_entry_mm=new, exact_rational_entry_mm=exact, injected50um_rejected=True))
    assert rows[0]['global_first_entry_mm'] == rows[1]['global_first_entry_mm']
    difference = abs(rows[0]['local_component_entry_mm'] - rows[1]['local_component_entry_mm'])
    assert difference >= 0.1
    empty = local_entry(p, n, [[0, 0, -0.3125]], [0.0625])[0]
    assert np.isnan(empty)
    return dict(rows=rows, summary_identity_error_mm=0.0, downstream_difference_mm=difference, minimum_extension='Connected interval component registered to the target t=0; global first ray root is insufficient', all_exact_oracle_errors_mm=0.0, all_injections_rejected=True, isolated_negative_branch_stays_UNKNOWN=True, resolution='PER_POINT', external_referent=dict(kind='closed_form', locator='https://mathworld.wolfram.com/Sphere.html equation7; exact_coaxial_oracle in code/local_ray_branch.py', compared_quantity='line/sphere roots and their exact rational local-component union; mathematical facit, not physical stock', refutes_us=True))

def run():
    path = ROOT / 'PREREG_R5.json'
    if not path.exists():
        pre = dict(round='R5', claim_type='capability', capability='Distinguish local intaglio/exterior cut boundary from a disconnected cutting-ball component on another wall', obstacle='R2 global first ray root can silently cross retained material; one safe-envelope exterior sample was labelled >25um overcut', changed_operation='Retain all sphere-ray intervals, form connected components, select the first component intersecting t>=0; disconnected negative components never define local stock', consumer=['DENT-MFG-ASBUILT-DEVIATION', _release_expand('X13')], resolution='PER_POINT', timescale='HANDOVER', metrics=dict(source_points=7680, summary_identity_error_mm_max=0.0, downstream_difference_mm_min=0.1, root_oracle_abs_mm_max=1e-06, corruption_mm=0.05, overcut_reporting_um=25.0, normal_ray_window_mm=[-0.4, 0.8]), decision_criteria='All source points get a local component or UNKNOWN; exact identical global-first-root pair must differ locally by >=.1mm; rational oracle and injected50um gates; original R2/coverage files unchanged', strongest_equal_information_control='Exact rational interval union for dyadic coaxial spheres with identical information, no method-superiority claim', falsifiers=['Wrong selected component compared with exact rational union', 'A disconnected entirely-negative component receives local stock instead of UNKNOWN', 'Any original prediction file changes'], external_referent=dict(kind='closed_form', locator='https://mathworld.wolfram.com/Sphere.html equation7', compared_quantity='oriented line/sphere interval, mathematical facit only', refutes_us=True), frozen_coverage_sha256=sha(ROOT / 'FROZEN_PREDICTIONS_SURFACE_COVERAGE.json'), full_cost=dict(preparation='Read frozen sparse centres/radii, no new geometric queries', fit=0, discovery='Per point sphere roots and connected interval union', validation='Exact rational dyadic oracle, identical-summary pair, injected50um', questions=0, fallback='UNKNOWN_NO_LOCAL_RAY_CUT_COMPONENT', physical='No CAM jobs or scans'), limitations=['Still sparse final-ball union, not CAM swept stock', 'No rigorous floating enclosure for general non-dyadic roots; exact control and R4 proof are narrower', 'Does not validate physical seating or manufactured stock'])
        dump(path, pre)
        path.with_suffix('.sha256').write_text(sha(path) + '\n')
        dump(ROOT / 'DECOMPOSITION_R5.json', dict(idea='Preserve which cut component belongs to a surface point', mechanism='An oriented ray can intersect several disconnected void/cut branches', equation='Sphere interval [a-sqrt(r²-d_perp²),a+sqrt(r²-d_perp²)]; connected union; discard components with high<0', operation='Sorted interval union with target-origin registration', representation='PER_POINT local ray component, original global root retained', leaves=[dict(leaf='sphere roots', status='DERIVED', stop_argument='Closed sphere equation; exact dyadic oracle provided, general rigorous float enclosure missing'), dict(leaf='locality', status='CONSTITUTIVE_CLOSURE', stop_argument='Choose zero-containing or nearest positive component; complete CAM material occupancy and manufacturing still UNKNOWN'), dict(leaf='physical stock', status='UNKNOWN', stop_argument='Requires actual CAM sweep or registered post-mill measurement')]))
    assert sha(path) == path.with_suffix('.sha256').read_text().strip()
    pre = read(path)
    assert sha(ROOT / 'FROZEN_PREDICTIONS_SURFACE_COVERAGE.json') == pre['frozen_coverage_sha256']
    start = time.perf_counter()
    ctl = controls()
    coverage = read(ROOT / 'RESULTS_SURFACE_COVERAGE.json')
    rows = []
    hashes = {}
    total_changed = 0
    state('R5_LOCAL_RAY_BRANCH_RUNNING', 'Frozen before component comparison', 'Keep original scalar fields and compute local-component fields')
    for (file, h) in coverage['raw_files'].items():
        assert sha(file) == h
        a = read(file)
        newrecords = []
        summaries = []
        for region in ['intaglio', 'exterior_and_rim']:
            rec = [p for p in a['records'] if p['region'] == region]
            P = np.array([p['point'] for p in rec])
            N = np.array([p['normal'] for p in rec])
            w = [p['safe_pose_witness'] for p in rec if p['safe_pose_witness'] is not None]
            centres = [x['centre'] for x in w]
            radii = [x['diameter_design_mm'] / 2 for x in w]
            safe = local_entry(P, N, centres, radii)
            tools = lib('vhf', 5)
            forced = local_entry(P, N, P + tools[0]['diameter_mm'] / 2 * N, np.full(len(P), tools[0]['diameter_mm'] / 2))
            old = np.array([np.nan if p['stock_safe_ball_union_um'] is None else p['stock_safe_ball_union_um'] / 1000 for p in rec])
            changed = (np.isfinite(old) != np.isfinite(safe)) | np.isfinite(old) & np.isfinite(safe) & (abs(old - safe) > 1e-06)
            total_changed += int(changed.sum())
            for (i, p) in enumerate(rec):
                h0 = np.nan if p['normal_gap_nominal_um'] is None else p['normal_gap_nominal_um']
                newrecords.append(dict(**p, legacy_global_first_stock_um=p['stock_safe_ball_union_um'], local_component_stock_um=safe[i] * 1000, local_component_forced_stock_um=forced[i] * 1000, local_component_normal_gap_um=h0 - safe[i] * 1000, local_component_forced_normal_gap_um=h0 - forced[i] * 1000, local_component_status='FINITE_BALL_LOCAL_COMPONENT' if np.isfinite(safe[i]) else 'UNKNOWN_NO_LOCAL_RAY_CUT_COMPONENT', component_changed=bool(changed[i])))
            summaries.append(dict(id=a['id'], region=region, n=len(rec), changed_points=int(changed.sum()), local_stock_gt25um_fraction=float(np.mean(safe > 0.025)), local_safe_overcut_gt25um_fraction=float(np.mean(safe < -0.025)), local_forced_overcut_gt25um_fraction=float(np.mean(forced < -0.025)), local_unknown_stock_fraction=float(np.mean(~np.isfinite(safe))), old_safe_overcut_gt25um_fraction=float(np.mean(old < -0.025)), resolution='PER_SURFACE_REGION', physical_gate='UNKNOWN'))
        dest = DATA / 'R5' / (a['id'] + '.json')
        dump(dest, dict(id=a['id'], prereg_sha256=sha(path), source_coverage_sha256=h, records=newrecords, regions=summaries))
        rows.extend(summaries)
        hashes[str(dest)] = sha(dest)
    for (file, h) in coverage['raw_files'].items():
        assert sha(file) == h
    result = dict(round='R5', claim_type='capability', controls=ctl, points=7680, changed_safe_points=total_changed, rows=rows, raw_files=hashes, seconds=time.perf_counter() - start, physical_gate='UNKNOWN', construction_gate='PASS', external_referent=pre['external_referent'], all_old_coverage_hashes_unchanged=True, rigorous_general_root_enclosure='MISSING')
    dump(ROOT / 'RESULTS_R5.json', result)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R5.json', dict(prereg_sha256=sha(path), result_sha256=sha(ROOT / 'RESULTS_R5.json'), point_files=hashes, code_sha256=sha(Path(__file__)), physical_measurement_performed=False))
    print('local ray branch', total_changed, 'safe points changed', result['seconds'])
if __name__ == '__main__':
    run()
