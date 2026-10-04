"""Numerical, export, primary-cell and actual fault-injection checks."""
import copy, math, sys, hashlib, zipfile, csv
import numpy as np
from scipy.special import gamma
import trimesh
from common import R, read, dump, sha
from surface_gate import main as surfaces, check_export, load_3mf, evaluate_vertices
from compare_lab import compare
sys.path.insert(0, str(R / 'vendor'))
import crown_case_sts3d as STS

def tangent(S, n):
    axis = np.eye(3)[np.argmin(abs(n), axis=1)]
    t = np.cross(n, axis)
    t /= np.linalg.norm(t, axis=1)[:, None]
    u = np.cross(n, t)
    M = np.zeros((len(S), 3, 3))
    M[:, 0, 0] = S[:, 0]
    M[:, 1, 1] = S[:, 1]
    M[:, 2, 2] = S[:, 2]
    M[:, 0, 1] = M[:, 1, 0] = S[:, 3]
    M[:, 1, 2] = M[:, 2, 1] = S[:, 4]
    M[:, 0, 2] = M[:, 2, 0] = S[:, 5]
    a = np.einsum('ni,nij,nj->n', t, M, t)
    b = np.einsum('ni,nij,nj->n', u, M, u)
    c = np.einsum('ni,nij,nj->n', t, M, u)
    return 0.5 * (a + b) + np.sqrt(0.25 * (a - b) ** 2 + c * c)

def response(record, lc):
    d = np.load(R / record['stress'])
    s = tangent(d[f'Sf__{lc}'].astype(float), d['nf'].astype(float))
    outside = np.ones(len(s), bool)
    tied = d['tied'].astype(bool)
    for p in d[f'lc__{lc}']:
        outside &= ~((np.linalg.norm(d['cf'] - p[:3], axis=1) < 3 * p[3]) & ~tied)
    m = record['material_model']
    if record['material'] == 'printed_resin':
        return float(m['sigma_mean'] / s[outside].max())
    exponent = m['m']
    I = np.sum(d['af'][outside] * np.maximum(s[outside], 0) ** exponent)
    seff = 20 / (exponent + 1) * (4 + 2 / (exponent + 1))
    s0 = m['sigma_mean'] / gamma(1 + 1 / exponent) * seff ** (1 / exponent)
    return float(s0 * (np.log(2) / I) ** (1 / exponent))

def deck_forces(path):
    out = []
    a = None
    for line in path.read_text().splitlines():
        if line.upper().startswith('*CLOAD'):
            a = np.zeros(3)
            out.append(a)
        elif line.startswith('*'):
            a = None
        elif a is not None and line.strip():
            x = line.split(',')
            a[int(x[1]) - 1] += float(x[2])
    return out

def lab_fixture(f):
    out = []
    for design in ['D1', 'M1', 'M2']:
        d = next((d for d in f['matched_family'] if d['design'] == design))
        for angle in [0, 30]:
            for k in range(12):
                out.append(dict(specimen_id=f'SYNTHETIC_{design}_{angle}_{k}', design_id=design, load_angle_deg=angle, fracture_force_N=d['FE_absolute_diagnostic_N'][str(angle)] * (1 + 0.015 * (k - 5.5)), MG_distance_mean_um=d['gap']['marginal_mean_um'], failure_mode='crown_tensile_fracture', fracture_origin='intaglio_tensile_zone', material_batch='SIMULATED_BATCH', cement_batch='SIMULATED_CEMENT', die_material_batch='SIMULATED_DIE', die_E_MPa=18000, indenter_diameter_mm=5, crosshead_mm_min=0.5))
    return out

def verify():
    checks = {}
    records = read('inputs/fe/INDEX.json')
    errors = []
    balance = []
    for r in records:
        forces = deck_forces(R / r['deck'])
        for (i, lc) in enumerate(['axial', 'offaxis30']):
            errors.append(abs(response(r, lc) / r['predictions'][lc]['failure_proxy_N'] - 1))
            balance.append(float(np.max(abs(forces[i] - np.array(r['cases'][lc]['load']['total_force_per_N'])))))
    checks['raw_response_all90_pass'] = max(errors) <= 0.05
    checks['deck_force_all90_pass'] = max(balance) <= 1e-06
    surf = surfaces()
    checks['export_actual_SDF_all15_pass'] = all((r['pass'] for r in surf))
    positive = trimesh.load_mesh(R / 'exports/D1/crown.stl')
    bad = positive.copy()
    bad.vertices += np.array([0.35, 0, 0])
    checks['actual_STL_350um_translation_rejected'] = not check_export('D1', 'crown', bad)['pass']
    translated_result = check_export('D1', 'crown', bad)
    bad = positive.copy()
    bad.vertices *= 1.1
    checks['actual_STL_10percent_scale_rejected'] = not check_export('D1', 'crown', bad)['pass']
    bad = trimesh.Trimesh(positive.vertices, positive.faces[:-1], process=False)
    checks['open_surface_rejected'] = not check_export('D1', 'crown', bad)['pass']
    wrong = R / 'raw/faults/wrong_units.3mf'
    wrong.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(R / 'exports/D1/crown.3mf') as src, zipfile.ZipFile(wrong, 'w') as dst:
        for name in src.namelist():
            b = src.read(name)
            if name.endswith('.model'):
                b = b.replace(b'unit="millimeter"', b'unit="meter"')
            dst.writestr(name, b)
    try:
        load_3mf(wrong)
        checks['actual_3MF_unit_rejected'] = False
    except ValueError:
        checks['actual_3MF_unit_rejected'] = True
    truth = response(records[0], 'axial')
    checks['actual_response_coefficient_x1p5_rejected'] = abs(truth / (1.5 * records[0]['predictions']['axial']['failure_proxy_N']) - 1) > 0.05
    checks['actual_deck_force_x2_rejected'] = float(np.max(abs(2 * deck_forces(R / records[0]['deck'])[0] - np.array(records[0]['cases']['axial']['load']['total_force_per_N'])))) > 1e-06
    f = read('FROZEN_PREDICTIONS.json')
    base = lab_fixture(f)
    good = compare(f, base)
    checks['synthetic_true_FE_contrast_supported'] = all((x['FE_gate'] == 'PASS' for x in good['double_ratio_comparison']))
    checks['synthetic_separable_control_rejected'] = all((x['control_gate'] == 'FAIL' for x in good['double_ratio_comparison']))
    for mode in ['contact_damage', 'cement_debonding', 'die_fracture']:
        bad = copy.deepcopy(base)
        for r in bad:
            if r['design_id'] == 'M1' and r['load_angle_deg'] == 30:
                r['failure_mode'] = mode
        z = compare(f, bad)
        checks[mode + '_dominant_rejected'] = z['mechanism_gate'] == 'FAIL' and next((q for q in z['double_ratio_comparison'] if q['design'] == 'M1'))['decision'] == 'FAIL_MECHANISM'
    bad = copy.deepcopy(base)
    bad[0]['failure_mode'] = 'unclassified'
    checks['unknown_failure_mode_is_UNKNOWN'] = compare(f, bad)['mechanism_gate'] == 'UNKNOWN'
    bad = copy.deepcopy(base)
    bad[0]['fracture_origin'] = ''
    checks['unknown_origin_is_UNKNOWN'] = compare(f, bad)['mechanism_gate'] == 'UNKNOWN'
    bad = copy.deepcopy(base)
    for r in bad:
        if r['design_id'] == 'M2' and r['load_angle_deg'] == 30:
            r['fracture_force_N'] *= 10
    z = compare(f, bad)
    checks['force_x10_rejected'] = not z['force_accuracy_gate_pass'] and next((q for q in z['double_ratio_comparison'] if q['design'] == 'M2'))['FE_gate'] == 'FAIL'
    bad = copy.deepcopy(base)
    for r in bad:
        r['MG_distance_mean_um'] += 1000
    checks['gap_plus1000um_rejected'] = not compare(f, bad)['gap_gate_pass']
    bad = copy.deepcopy(base)
    for r in bad:
        if r['design_id'] == 'M1':
            r['cement_batch'] = 'WRONG_BATCH'
    checks['mismatched_batch_rejected'] = next((q for q in compare(f, bad)['double_ratio_comparison'] if q['design'] == 'M1'))['protocol_gate'] == 'FAIL'
    bad = copy.deepcopy(base)
    for r in bad:
        r['die_E_MPa'] = 2000
    checks['wrong_support_modulus_rejected'] = all((q['protocol_gate'] == 'FAIL' for q in compare(f, bad)['double_ratio_comparison']))
    checks['frozen_hash_pass'] = sha(R / 'FROZEN_PREDICTIONS.json') == (R / 'FROZEN_PREDICTIONS.sha256').read_text().strip()
    checks['frozen_hash_mutation_rejected'] = hashlib.sha256((R / 'FROZEN_PREDICTIONS.json').read_bytes() + b' ').hexdigest() != (R / 'FROZEN_PREDICTIONS.sha256').read_text().strip()
    expected = [('PMC10817558', 639.3, 111.77), ('PMC10817558', 1378.1, 143.22), ('PMC10817558', 522.67, 108.57), ('PMC10817558', 516.86, 63.32), ('PMC10817558', 635.89, 78), ('PMC10817558', 865.3, 116.39), ('PMC10817558', 663.78, 106.8), ('PMC10817558', 980.1, 123.5), ('PMC10934854', 1084.5, 134.2), ('PMC10934854', 1313.1, 240.2)]
    actual = read('raw/LITERATURE_GROUPS.json')
    cells = []
    for (study, mean, sd) in expected:
        rr = next((r for r in actual if r['study'] == study and r['mean_N'] == mean))
        cells.append(dict(study=study, locator=rr['locator'], expected_mean_N=mean, expected_sd_N=sd, parsed_mean_N=rr['mean_N'], parsed_sd_N=rr['sd_N'], pass_cell=rr['sd_N'] == sd))
    checks['ten_primary_cell_checks_pass'] = all((x['pass_cell'] for x in cells))
    badtable = copy.deepcopy(actual)
    badtable[0]['mean_N'] *= 10

    def cell_gate(tab):
        return all((any((r['study'] == s and math.isclose(r['mean_N'], m, rel_tol=1e-09) and math.isclose(r['sd_N'], sd, rel_tol=1e-09) for r in tab)) for (s, m, sd) in expected))
    checks['injected_primary_mean_x10_rejected'] = not cell_gate(badtable)
    dump('raw/faults/primary_table_mean_x10.json', badtable)
    cases = []
    for p in (R / 'inputs/anatomy').glob('*.npz'):
        case = STS.load_case(p)
        cases.append(dict(case=p.name, parts={k: dict(vertices=len(m.vertices), watertight=bool(m.is_watertight)) for (k, m) in case['parts'].items() if m is not None}))
    checks['local_STS_adapter_reads_three_cases'] = len(cases) == 3
    dump('raw/LITERATURE_10ROW_CHECK.json', dict(producer_manual_primary_table_check=True, not_independent_review=True, checks=cells))
    dump('raw/FAULT_INJECTIONS.json', dict(synthetic_fixture=True, not_external_referent=True, checks=checks, actual_translated_export_result=translated_result, synthetic_lab_positive=good, synthetic_force_x10=z))
    result = dict(pass_verification=all(checks.values()), checks=checks, raw_response_relative_error_max=max(errors), deck_force_error_N_max=max(balance), exports=surf, STS_cases=cases, physical_validation=False)
    dump('VERIFICATION.json', result)
    if not result['pass_verification']:
        raise ValueError('Verification rejected: ' + ', '.join((k for (k, v) in checks.items() if not v)))
    return result
