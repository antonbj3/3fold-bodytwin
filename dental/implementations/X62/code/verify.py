import csv, json, copy, math
from fractions import Fraction
import numpy as np
from common import ROOT, DATA, save, sha, article, text
from patient import load
from beam import exact_binary_box
from measurement import assess_pose, assess_force
from run_round1 import source_cells, moment_certificate
from run_round3 import reciprocal_gate
from fixture_port import query, FIXTURE

def valid_cells(cells):
    for c in cells:
        if c['unit'] not in ['mm', 'degree']:
            return False
        if not all((math.isfinite(c[k]) for k in ['mean', 'sd', 'threshold'])):
            return False
        if c['sd'] < 0 or c['threshold'] <= 0:
            return False
    return True

def check_source(cells):
    gold = {('tab4', 1): [-0.26, 0.12, -0.29, 0.15, -0.01, 0.24], ('tab4', 2): [-0.15, 0.13, -0.24, 0.11, 0.01, 0.22], ('tab4', 3): [-0.1, 0.15, -0.12, 0.22, -0.01, 0.26], ('tab5', 1): [-0.06, 0.14, -0.06, 0.11, -0.18, 1.22], ('tab5', 2): [0.1, 0.39, -0.17, 0.4, 0.49, 1.39], ('tab5', 3): [0.16, 0.52, -0.05, 0.56, 0.25, 0.93]}
    if not valid_cells(cells):
        return False
    for ((tab, typ), a) in gold.items():
        rows = [c for c in cells if c['table'] == tab and c['tooth_type'] == typ]
        if len(rows) != 3:
            return False
        b = [x for c in rows for x in [c['mean'], c['sd']]]
        if max((abs(x - y) for (x, y) in zip(a, b))) > 5e-12:
            return False
    return True

def source_force_gate(pred, reference):
    return all((abs(p - o) / abs(o) <= 0.15 for (p, o) in zip(pred, reference)))

def prove_box(S, center, radius, lo, hi):
    for (i, row) in enumerate(S):
        for sign in [-1, 1]:
            point = [Fraction(float(u)) + sign * (1 if c >= 0 else -1) * Fraction(float(r)) for (c, u, r) in zip(row, center.ravel(), radius.ravel())]
            y = sum((Fraction(float(c)) * x for (c, x) in zip(row, point)))
            if not Fraction(float(lo[i])) <= y <= Fraction(float(hi[i])):
                return False
    return True

def run():
    (_, _, rows, meshes) = load()
    cells = source_cells()
    port = np.load(DATA / 'patient_port_R2.npz')
    S = port['S']
    r1 = json.loads((ROOT / 'raw/R1.json').read_text())
    r2 = json.loads((ROOT / 'raw/R2.json').read_text())
    tests = []

    def add(name, passed, kind):
        tests.append(dict(name=name, pass_gate=bool(passed), kind=kind))

    def rejected(fn):
        try:
            fn()
        except (ValueError, KeyError, AssertionError):
            return True
        return False
    add('independent36_source_mean_SD_cells', check_source(cells), 'normal')
    bad = copy.deepcopy(cells)
    bad[0]['mean'] += 10
    add('source_wrong_mean_rejected', not check_source(bad), 'injection')
    bad = copy.deepcopy(cells)
    bad[0]['sd'] = -1
    add('negative_SD_rejected', not valid_cells(bad), 'injection')
    bad = copy.deepcopy(cells)
    bad[0]['unit'] = 'N'
    add('incorrect_source_unit_rejected', not valid_cells(bad), 'injection')
    add('source_vertical_physical_closure_fails', not r1['O18']['tests'][1]['gate_pass'], 'normal_negative')
    add('injected_wrong_force_prediction_fails', not source_force_gate([12.7, 9.5], [1.27, 0.95]), 'injection')
    add('fixture_endpoint_replay', query('vertical', 0.02, FIXTURE)['force_linear_scenario_N'] == 0.95, 'normal_calibration')
    add('fixture_transfer_rejected', rejected(lambda : query('vertical', 0.02, 'titanium_bonded_patient')), 'injection')
    add('fixture_extrapolation_rejected', rejected(lambda : query('vertical', 0.2, FIXTURE)), 'injection')
    add('fixture_units_rejected', rejected(lambda : query('vertical', 0.02, FIXTURE, unit='m')), 'injection')
    with open(ROOT / 'raw/VIRTUAL_FORCE_IDENTITY.csv') as f:
        force = list(csv.DictReader(f))
    h = sha(ROOT / 'FROZEN_PREDICTIONS.json')
    add('virtual_force_identity', assess_force(force, h)['pass_gate'], 'normal_virtual')
    bad = copy.deepcopy(force)
    bad[0]['Fx_N'] = float(bad[0]['Fx_N']) + 1
    add('one_N_force_error_fails', not assess_force(bad, h)['pass_gate'], 'injection')
    add('wrong_prediction_hash_rejected', rejected(lambda : assess_force(force, '0' * 64)), 'injection')
    add('missing_mode_row_rejected', rejected(lambda : assess_force(force[:-1], h)), 'injection')
    bad = copy.deepcopy(force)
    bad[0]['Fx_N'] = 'nan'
    add('nonfinite_force_rejected', rejected(lambda : assess_force(bad, h)), 'injection')
    with open(ROOT / 'raw/VIRTUAL_POSE_IDENTITY.csv') as f:
        pose = list(csv.DictReader(f))
    add('nominal_full_pose_pass', assess_pose(pose, rows, meshes)['all_six_pass'], 'normal_virtual')
    bad = copy.deepcopy(pose)
    bad[2]['torque_deg'] = 4
    add('hidden_four_degree_rotation_fails', not assess_pose(bad, rows, meshes)['all_six_pass'], 'injection')
    bad = copy.deepcopy(pose)
    bad[0]['mesiodistal_mm'] = 0.6
    add('linear_limit_violation_fails', not assess_pose(bad, rows, meshes)['all_six_pass'], 'injection')
    add('missing_FDI_rejected', rejected(lambda : assess_pose(pose[:-1], rows, meshes)), 'injection')
    bad = copy.deepcopy(pose)
    bad[0]['FDI'] = bad[1]['FDI']
    add('duplicated_FDI_rejected', rejected(lambda : assess_pose(bad, rows, meshes)), 'injection')
    bad = copy.deepcopy(pose)
    bad[0]['rotation_deg'] = 'inf'
    add('nonfinite_pose_rejected', rejected(lambda : assess_pose(bad, rows, meshes)), 'injection')
    center = np.array(r2['force_box']['center_mm'])
    radius = np.array(r2['force_box']['radius_mm'])
    lo = np.array(r2['force_box']['lower_N'])
    hi = np.array(r2['force_box']['upper_N'])
    add('exact_binary_affine_enclosure', prove_box(S, center, radius, lo, hi), 'normal')
    broken = hi.copy()
    broken[0] = lo[0] - 1
    add('incorrect_affine_upper_bound_rejected', not prove_box(S, center, radius, lo, broken), 'injection')
    add('reciprocal_port_accepts_symmetric_PSD', reciprocal_gate(S), 'normal')
    bad = S.copy()
    bad[0, 5] += 10
    add('nonreciprocal_port_rejected', not reciprocal_gate(bad), 'injection')
    add('negative_stiffness_rejected', not reciprocal_gate(-np.eye(18)), 'injection')
    rng = np.random.default_rng(6202)
    max_error = -1e+100
    for trial in range(12):
        case = copy.deepcopy(pose)
        for row in case:
            for k in ['mesiodistal_mm', 'occlusogingival_mm', 'buccolingual_mm']:
                row[k] = float(rng.uniform(-0.5, 0.5))
            for k in ['torque_deg', 'rotation_deg', 'angulation_deg']:
                row[k] = float(rng.uniform(-6, 6))
        result = assess_pose(case, rows, meshes)
        max_error = max(max_error, max((x['whole_tooth_point_displacement_mm'] - x['rigorous_rigid_enclosure_mm'] for x in result['per_tooth'])))
    add('72_fullmesh_rigid_transforms_enclosed', max_error <= 0, 'normal')
    add('rigid_force_moment_balance', r1['O18']['numerical']['net_force_N'] < 1e-08 and r1['O18']['numerical']['net_moment_Nmm'] < 1e-08, 'normal')
    add('exact_summary_pairs', all([r1['O18']['sufficiency']['identity_error_mm'] == 0, r1['O20']['moment_sufficiency']['identity_error'] == 0, r1['O20']['pose_sufficiency']['center_summary_identity_error_mm'] == 0]), 'normal')
    inj = [t for t in tests if t['kind'] == 'injection']
    out = dict(tests=tests, passed=sum((t['pass_gate'] for t in tests)), total=len(tests), injected_faults_rejected=sum((t['pass_gate'] for t in inj)), injected_faults_total=len(inj), all_pass=all((t['pass_gate'] for t in tests)), worst_transform_minus_enclosure_mm=max_error, source_validation='manual local primary table reading, not independent clinician review', physical_measurement='NOT_RUN')
    save('raw/VERIFICATION.json', out)
    print(json.dumps(out, indent=2))
    if not out['all_pass']:
        raise SystemExit(1)
if __name__ == '__main__':
    run()
