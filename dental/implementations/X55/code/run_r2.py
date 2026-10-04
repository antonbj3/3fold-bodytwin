"""Second construction: external datum gauge and signed/location sufficiency."""
import datetime
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from metrology import analyze, sha, write_json, kabsch, transform
from run_demo import check_frozen, known_deformation
R = Path(__file__).resolve().parents[1]

def same_regional_summary():
    a = np.tile([40.0, -40.0], 4)
    b = -a
    summaries_a = [dict(bias_um=float(a[i:i + 2].mean()), rms_um=float(np.sqrt(a[i:i + 2] @ a[i:i + 2] / 2))) for i in range(0, 8, 2)]
    summaries_b = [dict(bias_um=float(b[i:i + 2].mean()), rms_um=float(np.sqrt(b[i:i + 2] @ b[i:i + 2] / 2))) for i in range(0, 8, 2)]
    nominal = np.array([30.0, 100.0])
    ca = nominal - a[:2]
    cb = nominal - b[:2]
    identity = max((abs(x[k] - y[k]) for (x, y) in zip(summaries_a, summaries_b) for k in x))
    return dict(summary='Mean signed bias and RMS in every named region', summaries_a=summaries_a, summaries_b=summaries_b, identity_error_um=identity, downstream_quantity='Minimum local clearance across two marginal sites with fixed nominal clearance', nominal_clearance_um=nominal.tolist(), clearance_a_um=ca.tolist(), clearance_b_um=cb.tolist(), minimum_a_um=float(ca.min()), minimum_b_um=float(cb.min()), difference_um=abs(float(ca.min() - cb.min())), sufficient=False, minimal_extension='Signed deviation tied to each decision-relevant position and its nominal clearance', resolution='PER_POINT', timescale='HANDOVER', reference='Exact fixed-correspondence local clearance formula; synthetic scenario')

def identifiability():
    u1 = np.full((5, 8), 40.0)
    beta1 = np.zeros_like(u1)
    u2 = np.zeros_like(u1)
    beta2 = np.full_like(u1, 40.0)
    noise = np.arange(-2, 3)[:, None] * np.ones((1, 8))
    observed1 = u1 + beta1 + noise
    observed2 = u2 + beta2 + noise
    return dict(observation_identity_error_um=float(abs(observed1 - observed2).max()), downstream_true_manufacture_difference_um=float(abs(u1 - u2).max()), same_repeatability_sd_um=float(observed1.std(0, ddof=1)[0]), identifiable=False, missing_measurement='Independent calibrated surface reference of this SAME physical part, plus an external datum', systematic_bias_debt=dict(resolution='PER_POINT', status='UNKNOWN', replaces='Phenomenological scanner bias/noise transfer from X10 summaries'), reference='Observation equation d=u+scanner_bias+epsilon; exact algebra', synthetic_not_empirical=True)

def run():
    start = time.perf_counter()
    for f in ['PREREG_R2.json', 'DECOMPOSITION_R2.json', 'FROZEN_PREDICTIONS_R2.json']:
        check_frozen(f)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out = R / 'replays' / ('R2_' + stamp)
    out.mkdir(parents=True)
    reference = np.array([[-7.0, -7.0, -1.0], [7.0, -7.0, -1.0], [-7.0, 7.0, -1.0], [7.0, 7.0, 3.0]])
    datum_results = {}
    controls = {}
    wrong = []
    for name in ['D1', 'M2']:
        path = R / 'raw' / name
        meta = json.loads((path / 'metadata.json').read_text())
        for row in meta['scans'].values():
            T = np.array(row['scan_to_reference'])
            a = transform(reference, np.linalg.inv(T))
            row['scan_datum_points_mm'] = a.tolist()
            row['reference_datum_points_mm'] = reference.tolist()
            row.pop('scan_to_reference')
        write_json(out / (name + '_datum_metadata.json'), meta)
        rep = analyze(path / 'design.stl', [path / f'repeat_{j}.stl' for j in range(5)], path / 'regions.json', out / name / 'datum', mode='datum', metadata=out / (name + '_datum_metadata.json'))
        row = meta['scans']['repeat_0.stl']
        a = np.array(row['scan_datum_points_mm'])
        b = np.array(row['reference_datum_points_mm'])
        T = kabsch(a, b)
        control_start = time.perf_counter()

        def residual(x):
            return (a @ Rotation.from_rotvec(x[:3]).as_matrix().T + x[3:] - b).ravel()
        fit = least_squares(residual, np.zeros(6), xtol=1e-13, ftol=1e-13, gtol=1e-13)
        controls[name] = dict(kabsch_rms_um=float(np.sqrt(np.mean((transform(a, T) - b) ** 2)) * 1000), nonlinear_rms_um=float(np.sqrt(np.mean(fit.fun ** 2)) * 1000), point_difference_um=float(abs(transform(a, T) - (a @ Rotation.from_rotvec(fit.x[:3]).as_matrix().T + fit.x[3:])).max() * 1000), seconds=time.perf_counter() - control_start)
        datum_results[name] = rep['regions']
        bad = json.loads((path / 'metadata.json').read_text())
        badT = np.array(bad['scans']['rigid.stl']['scan_to_reference'])
        badT[2, 3] += 0.1
        bad['scans']['rigid.stl']['scan_to_reference'] = badT.tolist()
        write_json(out / (name + '_bad_datum.json'), bad)
        rejected = analyze(path / 'design.stl', [path / 'rigid.stl'], path / 'regions.json', out / name / 'bad_datum', mode='datum', metadata=out / (name + '_bad_datum.json'))
        error = max((r['pooled_normal_rms_um'] for r in rejected['regions'].values()))
        wrong.append(dict(name='known_wrong_datum_100_um_' + name, measured_rms_um=error, rejected_by_frozen_rigid_gate=bool(error > 2)))
        print(name, 'datum occlusal bias', rep['regions']['occlusal']['signed_bias_um'], flush=True)
    sf = same_regional_summary()
    ident = identifiability()
    import sys
    sys.path.insert(0, str(R / 'tests'))
    from test_sources import primary_gate
    cells = json.loads((R / 'raw/LITERATURE_CELLS.json').read_text())['cells']
    source_gate = primary_gate(cells)
    broken = [dict(c) for c in cells]
    broken[0]['mean_um'] += 1
    source_fault = not primary_gate(broken)
    known = json.loads((R / 'R1_RESULTS.json').read_text())['known_deformations']
    gates = dict(datum_landmark_recovery=all((c['kabsch_rms_um'] <= 1e-06 for c in controls.values())), known_field_rmse=all((v['normal_field_rmse_um'] <= 8 for n in known.values() for v in n.values())), regional_summary_insufficient=sf['identity_error_um'] == 0 and sf['difference_um'] >= 30, repeat_bias_nonidentifiable=ident['observation_identity_error_um'] == 0, source_cell_reproduction=source_gate, injected_wrong_value_rejects=source_fault and all((f['rejected_by_frozen_rigid_gate'] for f in wrong)))
    r1 = json.loads((R / 'R1_RESULTS.json').read_text())
    result = dict(claim_type='capability', status='DATUM_POINT_FIELD_CAPABILITY_WITH_PHYSICAL_IDENTIFIABILITY_LIMIT', gates=gates, datum_regions=datum_results, datum_control=controls, known_field_rmse_from_R1_same_data=known, sufficiency=sf, scanner_manufacture_identifiability=ident, faults=wrong, source_cell_fault_rejected=source_fault, external_referent=r1['external_referent'], rigorous_uncertainty='Pose/reference/correspondence enclosure MISSING for real scans; repeat-only CIs conditional', physical_validation='UNKNOWN_NO_REAL_ASBUILT_PAIR', frozen_prereg_sha256=sha(R / 'PREREG_R2.json'), replay_output=str(out), full_cost={'executed_seconds': time.perf_counter() - start, 'acquisition_and_discovery': 'UNKNOWN'}, resolution=['PER_POINT', 'PER_SURFACE_REGION', 'PER_TOOTH'], timescale='HANDOVER')
    write_json(out / 'R2_RESULTS.json', result)
    write_json(R / 'R2_RESULTS.json', result)
    combined = dict(claim_type='capability', status='K48_CLI_IMPLEMENTED_PHYSICAL_VALIDATION_UNKNOWN', r1=r1, r2=result, external_referent=r1['external_referent'], resolution=['PER_POINT', 'PER_SURFACE_REGION', 'PER_TOOTH', 'POPULATION'], review_state='PENDING_INDEPENDENT_REVIEW')
    write_json(R / 'results.json', combined)
    write_json(R / 'CURRENT_WORK_STATE.json', dict(status='R2_COMPLETE', latest_gate=gates, next_operation='Independent actual-part reference and external datum with known uncertainty; then held-out K43 calibration', physical_measurements_acquired=0))
    print('R2 gates', gates, flush=True)
    return result
if __name__ == '__main__':
    run()
