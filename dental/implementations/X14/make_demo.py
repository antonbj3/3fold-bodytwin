"""Offline, one-command reproduction of frozen X14 constructions."""
import argparse
import csv
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import resource
import time
import warnings
import xml.etree.ElementTree as ET
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from source_tables import ROOT, expanded_rows, hirano_angles, shrinkage_cells
from sinter_position import angular_gap_um, distortion_decision
from seating import profile, relative_angle_deg, seat_vertices, seat_lp, measured_seating_bound

def read(name):
    return json.loads((ROOT / name).read_text())

def write(name, value):
    (ROOT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def frozen(name):
    checksum = (ROOT / name.replace('.json', '.sha256')).read_text().split()[0]
    assert sha(ROOT / name) == checksum, f'Frozen prediction drift: {name}'
    value = read(name)
    round_id = value.get('round', name.removeprefix('FROZEN_PREDICTIONS_').removesuffix('.json'))
    if round_id == 'R1':
        prereg = 'PREREG_R1.json'
    else:
        prereg = f'PREREG_{round_id}.json'
    assert value['prereg_sha256'] == sha(ROOT / prereg), f'Prereg drift: {prereg}'
    return value

def rmse(p, observed):
    return float(np.sqrt(np.mean((np.asarray(p) - np.asarray(observed)) ** 2)))

def gate_r1(p, y, base):
    (error, baseline) = (rmse(p, y), rmse(base, y))
    return error <= 0.25 and error <= 0.7 * baseline

def gate_r2(p, y, iso, axis):
    error = rmse(p, y)
    return (error <= 0.2 and error <= 0.7 * rmse(iso, y), error <= 0.8 * rmse(axis, y))

def gate_r3(p, y, raw, scalar):
    error = rmse(p, y)
    return (error <= 0.2 and error <= 0.5 * rmse(raw, y), error <= 0.8 * rmse(scalar, y))

def run_r1():
    f = frozen('FROZEN_PREDICTIONS.json')
    assert f['source_sha256'] == sha(ROOT / 'sources/Hirano2025.xml')
    (train, test) = (hirano_angles(3), hirano_angles(4))
    for (row, p) in zip(train, f['predictions']):
        assert row['area'] == p['area'] and 'c-MCL-' + row['manufacturer'] == p['material']
        assert abs(row['delta_deg'] - p['candidate_deg']) < 1e-12
        mean = np.mean([r['delta_deg'] for r in train if r['manufacturer'] == row['manufacturer']])
        assert abs(mean - p['pooled_control_deg']) < 1e-12
    y = [x['delta_deg'] for x in test]
    p = [x['candidate_deg'] for x in f['predictions']]
    base = [x['pooled_control_deg'] for x in f['predictions']]
    r = dict(candidate_deg=rmse(p, y), pooled_control_deg=rmse(base, y), isotropic_control_deg=rmse(np.zeros(len(y)), y))
    r['candidate_to_pooled_ratio'] = r['candidate_deg'] / r['pooled_control_deg']
    r['accuracy_gate'] = 'PASS' if gate_r1(p, y, base) else 'FAIL'
    r['specimen_interval_coverage'] = 'UNKNOWN: no specimen SD/raw list, table statistic unnamed'
    r['future_total_gap'] = 'UNKNOWN: no seated-gap measurement'
    design = np.eye(6)
    lookup = design @ np.linalg.lstsq(design, [x['delta_deg'] for x in train], rcond=None)[0]
    r['equally_informed_categorical_LS_max_difference_deg'] = float(np.max(np.abs(lookup - p)))
    r['strongest_control_outcome'] = 'TIE'
    (x, z) = profile(0.932, 0.5)
    r['isotropic_angle_invariance_error_deg'] = max((abs(relative_angle_deg(x * ef, z * ef) - 0.932) for ef in [1.225, 1.231]))
    write('R1_RESULTS.json', r)
    return (r, train + test, dict(r1_10x_angle_rejected=not gate_r1(np.asarray(p) * 10, y, base), r1_zero_angle_rejected=not gate_r1(np.zeros(len(y)), y, base), r1_control_injected_10x_prediction_rejected=rmse(np.asarray(base) * 10, y) > 0.25))

def run_r2(cells):
    f = frozen('FROZEN_PREDICTIONS_R2.json')
    assert f['source_sha256'] == sha(ROOT / 'sources/Shrinkage2025.xml')
    rows = f['predictions']
    y = [cells[x['cell_id']]['mean_shrinkage_pct'] for x in rows]
    keys = ['prediction_pct', 'isotropic_pct', 'anisotropic_no_position_pct']
    predictions = {k: [x[k] for x in rows] for k in keys}
    for p in rows:
        target = cells[p['cell_id']]
        group = [cells[i] for i in p['training_cell_ids']]
        assert all((c['position'] != target['position'] for c in group))
        axis = [c for c in group if c['axis'] == target['axis']]
        assert abs(sum((w * c['mean_shrinkage_pct'] for (w, c) in zip(p['weights'], axis))) - p['prediction_pct']) < 1e-10
    r = {k: rmse(predictions[k], y) for k in keys}
    r['position_to_isotropic_ratio'] = r['prediction_pct'] / r['isotropic_pct']
    r['position_to_anisotropic_ratio'] = r['prediction_pct'] / r['anisotropic_no_position_pct']
    (accuracy, position) = gate_r2(predictions[keys[0]], y, predictions[keys[1]], predictions[keys[2]])
    r.update(accuracy_gate='PASS' if accuracy else 'FAIL', position_gain_gate='PASS' if position else 'FAIL', interval_coverage='UNKNOWN: n/covariance per cell not supplied', heldout_cells=len(rows))
    r['by_method'] = {}
    for method in ['micrometer', 'light_microscopy', 'surface_scan']:
        ids = [i for (i, x) in enumerate(rows) if x['method'] == method]
        r['by_method'][method] = {k: rmse([predictions[k][i] for i in ids], [y[i] for i in ids]) for k in keys}
    write('R2_RESULTS.json', r)
    return (r, dict(r2_10x_shrinkage_rejected=not all(gate_r2(np.asarray(predictions[keys[0]]) * 10, y, predictions[keys[1]], predictions[keys[2]])), r2_controls_injected_10x_rejected=all((rmse(np.asarray(predictions[k]) * 10, y) > 0.2 for k in keys[1:]))))

def run_r3(cells):
    rows = frozen('FROZEN_PREDICTIONS_R3.json')['predictions']
    key = lambda x: (x['material'], x['position'], x['horizontal'], x['axis'])
    gauge = {key(c): c for c in cells if c['method'] == 'micrometer'}
    y = [gauge[tuple(p['key'])]['mean_shrinkage_pct'] for p in rows]
    max_ls = 0.0
    for material in sorted({p['key'][0] for p in rows}):
        train = [c for c in cells if c['material'] == material and c['method'] == 'surface_scan' and (c['position'] == 'Middle')]
        design = np.array([[int(c['axis'] == axis) for axis in ['x', 'y', 'z']] for c in train])
        bias = np.array([c['mean_shrinkage_pct'] - gauge[key(c)]['mean_shrinkage_pct'] for c in train])
        fit = np.linalg.lstsq(design, bias, rcond=None)[0]
        for p in rows:
            if p['key'][0] == material:
                pred = p['raw_scan_pct'] - fit[['x', 'y', 'z'].index(p['key'][3])]
                max_ls = max(max_ls, abs(pred - p['axis_corrected_pct']))
    keys = ['raw_scan_pct', 'axis_corrected_pct', 'scalar_corrected_pct']
    pp = {k: [p[k] for p in rows] for k in keys}
    r = {k: rmse(pp[k], y) for k in keys}
    r.update(axis_to_raw_ratio=r['axis_corrected_pct'] / r['raw_scan_pct'], axis_to_scalar_ratio=r['axis_corrected_pct'] / r['scalar_corrected_pct'], heldout_cells=len(rows), strongest_LS_control_max_difference_pct=max_ls, future_batch_transfer='UNKNOWN', premanufacture_prediction='NO: query scan required')
    (a, b) = gate_r3(pp[keys[1]], y, pp[keys[0]], pp[keys[2]])
    r.update(accuracy_gate='PASS' if a else 'FAIL', axis_gain_gate='PASS' if b else 'FAIL')
    write('R3_RESULTS.json', r)
    return (r, dict(r3_10x_corrected_scan_rejected=not all(gate_r3(np.asarray(pp[keys[1]]) * 10, y, pp[keys[0]], pp[keys[2]])), r3_controls_injected_10x_rejected=all((rmse(np.asarray(pp[k]) * 10, y) > 0.2 for k in [keys[0], keys[2]]))))

def run_r4():
    f = frozen('FROZEN_PREDICTIONS_R4.json')
    rows = []
    for case in f['cases']:
        (x, z) = (np.array(case['x_mm']), np.array(case['z_mm']))
        start = time.perf_counter()
        v = seat_vertices(x, z)
        tv = time.perf_counter() - start
        start = time.perf_counter()
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='Unrecognized options detected')
            lp = seat_lp(x, z)
        tl = time.perf_counter() - start
        assert abs(v['max_gap_mm'] - case['candidate_predicted_gap_mm']) < 1e-12
        rows.append(dict(**case, vertex_result=v, published_LP_result=lp, vertex_seconds=tv, LP_seconds=tl))
    diff = max((abs(c['vertex_result']['max_gap_mm'] - c['published_LP_result']['max_gap_mm']) for c in rows))
    gaps = [c['vertex_result']['max_gap_mm'] * 1000 for c in rows if c['alpha_deg'] != 0]
    angle_error = max((abs(relative_angle_deg(c['x_mm'], c['z_mm']) - c['alpha_deg']) for c in rows))
    r = dict(cases=rows, max_candidate_LP_difference_mm=diff, solver_gate='PASS' if diff <= 1e-08 else 'FAIL', same_angle_gap_range_um=[min(gaps), max(gaps)], ambiguity_gate='PASS' if max(gaps) - min(gaps) >= 20 else 'FAIL', relative_angle_error_deg=angle_error, strongest_control_outcome='TIE', physical_gap_validation='UNKNOWN', source_alpha_deg=0.932, source_D_mm=[7, 11], proposed_span_mm=24, span_provenance='our own lab fixture, not source measurement')
    write('R4_RESULTS.json', r)
    return (r, dict(r4_injected_zero_LP_answer_rejected=max(gaps) > 1e-05, r4_injected_10x_vertex_answer_rejected=max(gaps) * 9 / 1000 > 1e-08, r4_injected_10x_angle_rejected=abs(relative_angle_deg(rows[-1]['x_mm'], np.asarray(rows[-1]['z_mm']) * 10) - 0.932) > 1e-08))

def run_r5():
    f = frozen('FROZEN_PREDICTIONS_R5.json')
    errors = []
    covered = []
    raw = []
    for case in f['predictions']:
        (x, z) = (np.array(case['x_mm']), np.array(case['z_mm']))
        translated = z + 1.0 + 0.03 * x
        canonical = measured_seating_bound(x, translated, 0.01)
        lp = seat_lp(x, translated)
        errors.append(abs(canonical['max_gap_mm'] - lp['max_gap_mm']))
        assert abs(canonical['max_gap_mm'] - case['predicted_nominal_mm']) < 1e-10
        for signs in product([-1.0, 1.0], repeat=4):
            perturbed = z + 0.01 * np.array(signs)
            value = seat_lp(x, perturbed)['max_gap_mm']
            (lo, hi) = case['predicted_bound_mm']
            covered.append(lo - 1e-10 <= value <= hi + 1e-10)
            raw.append(dict(alpha_deg=case['alpha_deg'], fraction=case['fraction'], signs=signs, gap_mm=value, frozen_interval_mm=[lo, hi]))
    r = dict(canonical_vs_direct_LP_max_error_mm=max(errors), canonical_gate='PASS' if max(errors) <= 1e-08 else 'FAIL', covered_corners=sum(covered), total_corners=len(covered), bounded_error_gate='PASS' if all(covered) else 'FAIL', height_error_mm=0.01, height_error_provenance='chosen deterministic stress bound, not measured scanner error', strongest_equal_information_control='TIE: full-point SciPy LP', minimal_planar_observations='two chord-relative margin slopes; relative angle alone lacks one invariant', physical_three_dimensional_gap='UNKNOWN')
    write('R5_RESULTS.json', r)
    write('raw/R5_error_corners.json', raw)
    case = f['predictions'][3]
    (x, z) = (np.array(case['x_mm']), np.array(case['z_mm']))
    bad = seat_lp(x, z + np.array([0.15, -0.15, 0.15, -0.15]))['max_gap_mm']
    (lo, hi) = case['predicted_bound_mm']
    return (r, dict(r5_injected_150um_height_error_rejected=not lo - 1e-10 <= bad <= hi + 1e-10, r5_injected_missing_measurement_error_blocks=measured_seating_bound(x, z)['status'] == 'UNKNOWN_MEASUREMENT_ERROR', r5_LP_injected_10x_objective_rejected=case['predicted_nominal_mm'] * 9 > 1e-08))

def source_metadata():
    r = ET.parse(ROOT / 'sources/Hirano2025.xml').getroot()
    products = expanded_rows(r.find(".//table-wrap[@id='materials-18-04234-t001']/table"))
    programs = expanded_rows(r.find(".//table-wrap[@id='materials-18-04234-t002']/table"))
    write('raw/materials_and_programs.json', dict(product_table=products, program_table=programs, locator='doi:10.3390/ma18184234 Tables1/2', warning='Product-specific programs confounded; no separate program-effect coefficient'))

def decision_table(angles):
    out = []
    for row in angles:
        a = row['delta_deg']
        for diameter in [7.0, 11.0]:
            g = angular_gap_um(a, diameter)
            scenarios = {str(b): distortion_decision([a, a], diameter, [b, b])['status'] for b in [0, 30, 50]}
            out.append(dict(material=row['material'], area=row['area'], diameter_mm=diameter, observed_aggregate_delta_deg=a, equal_split_distortion_um=g, available_other_gap_budget_um=120 - g, scenario_other_0_um=scenarios['0'], scenario_other_30_um=scenarios['30'], scenario_other_50_um=scenarios['50'], robust_future_total_gap='UNKNOWN', locator=row['locator']))
    with (ROOT / 'NESTING_DECISIONS.csv').open('w', newline='') as file:
        w = csv.DictWriter(file, fieldnames=out[0])
        w.writeheader()
        w.writerows(out)
    write('NESTING_DECISIONS.json', out)
    return out

def plot(angles, r1, r2, r3, r4):
    (fig, axs) = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for material in ['n-MCL-A', 'c-MCL-A', 'n-MCL-B', 'c-MCL-B']:
        rows = [r for r in angles if r['material'] == material]
        axs[0, 0].plot(['Top I', 'Middle II', 'Bottom III'], [r['delta_deg'] for r in rows], marker='o', label=material)
    axs[0, 0].axhline(0, color='black', lw=0.6)
    axs[0, 0].set_ylabel('Published aggregate angular change (degrees)')
    axs[0, 0].set_title('Measured: position changes bending sign')
    axs[0, 0].legend(fontsize=8)
    for material in ['c-MCL-A', 'c-MCL-B']:
        rows = [r for r in angles if r['material'] == material]
        axs[0, 1].plot(['Top I', 'Middle II', 'Bottom III'], [angular_gap_um(r['delta_deg'], 11) for r in rows], marker='o', label=material)
    axs[0, 1].axhline(120, label='120 um total budget', color='red', ls='--')
    axs[0, 1].axhline(70, label='70 um remaining if other processes use 50', color='orange', ls=':')
    axs[0, 1].set_ylabel('Equal-split geometric contribution (um)')
    axs[0, 1].set_title('Calculated: gap contribution, no measured fit')
    axs[0, 1].legend(fontsize=8)
    values = [r2[k] for k in ['prediction_pct', 'isotropic_pct', 'anisotropic_no_position_pct']] + [r3[k] for k in ['raw_scan_pct', 'axis_corrected_pct', 'scalar_corrected_pct']]
    axs[1, 0].bar(range(6), values, color=['#c54b42', '#6589a5', '#6589a5', '#a78d68', '#c54b42', '#6589a5'])
    axs[1, 0].set_xticks(range(6), ['R2 position', 'R2 isotropic', 'R2 axes', 'R3 raw scan', 'R3 axes', 'R3 scalar'], rotation=25, ha='right', fontsize=8)
    axs[1, 0].set_ylabel('Held-out shrinkage RMSE (percentage points)')
    axs[1, 0].set_title('Negative: position fit and one-height calibration')
    cases = [c for c in r4['cases'] if c['alpha_deg'] != 0]
    axs[1, 1].plot([c['fraction'] for c in cases], [c['vertex_result']['max_gap_mm'] * 1000 for c in cases], 'o-')
    axs[1, 1].set_xlabel('Unobserved allocation relative to center chord')
    axs[1, 1].set_ylabel('Planar seated maximum opening (um)')
    axs[1, 1].set_title('Simulation: same 0.932 deg, different seated gap')
    fig.suptitle('X14 sinter position — measurements, closures and failures kept separate', fontsize=13)
    fig.savefig(ROOT / 'FIGURE.png', dpi=160)
    fig.savefig(ROOT / 'FIGURE.pdf')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--margin-csv', type=Path, help='Lab x_mm,z_mm CSV; optional measured planar seating query')
    parser.add_argument('--height-error-um', type=float, help='Externally established absolute z-error bound; fixed x assumed')
    args = parser.parse_args()
    start = time.perf_counter()
    assert sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file())) < 3000000000
    (ROOT / 'raw').mkdir(exist_ok=True)
    source_metadata()
    cells = shrinkage_cells()
    write('raw/shrinkage_cells.json', cells)
    timings = {}
    faults = {}
    a = time.perf_counter()
    (r1, angles, f1) = run_r1()
    timings['R1_seconds'] = time.perf_counter() - a
    faults.update(f1)
    a = time.perf_counter()
    (r2, f2) = run_r2(cells)
    timings['R2_seconds'] = time.perf_counter() - a
    faults.update(f2)
    a = time.perf_counter()
    (r3, f3) = run_r3(cells)
    timings['R3_seconds'] = time.perf_counter() - a
    faults.update(f3)
    a = time.perf_counter()
    (r4, f4) = run_r4()
    timings['R4_seconds'] = time.perf_counter() - a
    faults.update(f4)
    a = time.perf_counter()
    (r5, f5) = run_r5()
    timings['R5_seconds'] = time.perf_counter() - a
    faults.update(f5)
    write('raw/angles_all.json', angles)
    calculated = [angular_gap_um(0.932, d) for d in [11, 7]]
    geometric_errors = [abs(a - b) for (a, b) in zip(calculated, [89.4, 56.9])]
    faults['geometric_10x_angle_rejected'] = any((abs(angular_gap_um(9.32, d) - g) > 1 for (d, g) in zip([11, 7], [89.4, 56.9])))
    faults['geometric_zero_angle_rejected'] = any((abs(angular_gap_um(0, d) - g) > 1 for (d, g) in zip([11, 7], [89.4, 56.9])))
    faults['gap_policy_injected_10x_angle_exceeds'] = distortion_decision([9.32, 9.32], 11, [0, 0])['status'] == 'EXCEEDS_DISTORTION_BUDGET'
    faults['gap_policy_missing_uncertainty_blocks'] = distortion_decision(None, 11, None)['status'] == 'UNKNOWN'
    assert all(faults.values()), faults
    write('FAULT_INJECTION.json', faults)
    decisions = decision_table(angles)
    plot(angles, r1, r2, r3, r4)
    if args.margin_csv:
        with args.margin_csv.open() as file:
            data = list(csv.DictReader(file))
        (x, z) = ([float(r['x_mm']) for r in data], [float(r['z_mm']) for r in data])
        query = dict(input_path=str(args.margin_csv.resolve()), input_sha256=sha(args.margin_csv), vertex=seat_vertices(x, z), LP=seat_lp(x, z), two_invariant_bound=measured_seating_bound(x, z, None if args.height_error_um is None else args.height_error_um / 1000), empirical_certification='UNKNOWN: input error/3D seating validity must be supplied')
        write('LAB_QUERY_RESULT.json', query)
    results = dict(lane='X14-sinter-position', review_state='PENDING_INDEPENDENT_REVIEW', outcome='AGGREGATE_ANGLE_GAIN_WITH_TWO_FAILED_TRANSFER_CONSTRUCTIONS_AND_PLANAR_OBSERVABILITY_RESULT', protocol_limitations=read('PROTOCOL_LIMITATIONS.json')['limitations'], observation_cost_warning='R3 needs a new post-sinter scan for each query; R4/R5 need common-frame points and externally calibrated errors. No full physical cost or prospective accuracy claim.', external_referent=dict(kind='independent_measurement', locator='doi:10.3390/ma18184234 Tables3/4', compared_quantity='pre/post angular change of standardized four-unit FDP', refutes_us=False), external_referents=[dict(kind='independent_measurement', locator='doi:10.3390/ma18143217 Table2', compared_quantity='XYZ coupon shrinkage means and matched gauge/scan values', refutes_us=True), dict(kind='published_code', locator='https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html', compared_quantity='planar constrained maxgap minimum, numerical verification only', refutes_us=False), dict(kind='our_own_fixture', locator='PREREG_R4.json', compared_quantity='same-alpha four-point ambiguity; synthetic span24mm', refutes_us=False)], rounds={'R1': r1, 'R2': r2, 'R3': r3, 'R4': {k: v for (k, v) in r4.items() if k != 'cases'}, 'R5': r5}, geometric_reconstruction=dict(values_um=calculated, reference_um=[89.4, 56.9], errors_um=geometric_errors, gate='PASS' if max(geometric_errors) <= 1 else 'FAIL', kind='published calculation, not independent seated-gap measurement'), supported_regime='Published product/lot/program and standardized four-unit FDP for angular observations; no radial position interpolation or new-program prediction', unknowns=['future specimen angle uncertainty', 'total measured marginal fit', 'continuous height/radial interpolation', 'crown and long-span transfer', 'joint same-batch tensor plus curvature', '3D interference and seating'], costs=dict(**timings, total_demo_seconds=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, threads=1, GPU_seconds=0, physical_measurements_performed=0, source_preparation_wall_seconds=None, historical_discovery_cost='UNKNOWN', language_model_cost='UNKNOWN', physical_acquisition_cost='UNKNOWN', fallback='four measured margin heights and same-batch before/after XYZ baselines; otherwise UNKNOWN'), software=dict(numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__), fault_injection_pass=all(faults.values()), artifacts={}, large_arrays=[])
    for name in ['PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json', 'PREREG_R4.json', 'PREREG_R5.json', 'FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R2.json', 'FROZEN_PREDICTIONS_R3.json', 'FROZEN_PREDICTIONS_R4.json', 'FROZEN_PREDICTIONS_R5.json', 'source_tables.py', 'sinter_position.py', 'seating.py', 'make_demo.py', 'run_all.sh', 'sources/Hirano2025.xml', 'sources/Shrinkage2025.xml', 'raw/angles_all.json', 'raw/shrinkage_cells.json', 'NESTING_DECISIONS.csv', 'R4_RESULTS.json', 'R5_RESULTS.json', 'raw/R5_error_corners.json', 'FIGURE.png', 'FAULT_INJECTION.json']:
        results['artifacts'][name] = dict(sha256=sha(ROOT / name), bytes=(ROOT / name).stat().st_size)
    for name in ['README_DEMO.md', 'RESULTS.md', 'DECOMPOSITION.json', 'DECOMPOSITION_R2.json', 'DECOMPOSITION_R3.json', 'DECOMPOSITION_R4.json', 'DECOMPOSITION_R5.json', 'PROTOCOL_LIMITATIONS.json', 'LAB_MEASUREMENT_SPEC.json', 'LICENSES.json', 'FAILURES.json']:
        results['artifacts'][name] = dict(sha256=sha(ROOT / name), bytes=(ROOT / name).stat().st_size)
    source = ET.parse(ROOT / 'sources/Hirano2025.xml').getroot()
    assert '0.932' in ' '.join(source.find(".//sec[@id='sec3-materials-18-04234']").itertext()), 'Published maximum locator drift'
    write('results.json', results)
    print(f"R1 angle RMSE {r1['candidate_deg']:.6f} vs {r1['pooled_control_deg']:.6f} degrees: {r1['accuracy_gate']}")
    print(f"R2 position field {r2['prediction_pct']:.6f} vs isotropic {r2['isotropic_pct']:.6f} pct points: {r2['accuracy_gate']}")
    print(f"R3 one-height scan calibration {r3['axis_corrected_pct']:.6f} vs raw {r3['raw_scan_pct']:.6f} pct points: {r3['accuracy_gate']}")
    print(f"R4 same angle planar gap {min(r4['same_angle_gap_range_um']):.3f}..{max(r4['same_angle_gap_range_um']):.3f} um; LP {r4['strongest_control_outcome']}")
    print(f"R5 two invariants: {r5['covered_corners']}/{r5['total_corners']} corners inside frozen gap bound; physical 3D gap UNKNOWN")
    print('Total marginal fit: UNKNOWN. Saved results.json, NESTING_DECISIONS.csv and FIGURE.png.')
if __name__ == '__main__':
    main()
