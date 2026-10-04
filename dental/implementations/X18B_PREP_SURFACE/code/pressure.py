"""Population force scenarios on contact-force simplex; no invented pressure shares."""
from common import *
import resource
from scipy.optimize import linprog
SOURCE = H.parents[1] / 'notes/occlusion_support_intervals/SUPPORT_INTERVALS.json'

def tables():
    h = dict(types=['I1', 'I2', 'C', 'P1', 'P2', 'M1', 'M2'], mean_N=[15.0, 12.6, 19.7, 23.1, 30.5, 130.4, 211.8], sd_N=[25.7, 22.1, 28.6, 25.0, 31.8, 67.0, 111.3], M3_contact='absent', n_people=42, n_sides=65, locator='10.7144/sgf.2.111 Table1 printed page113, wo/M3 group')
    f = dict(fdi=[17, 16, 15, 14, 13, 12, 11, 21, 22, 23, 24, 25, 26, 27], mean_percent=[14.5, 14.7, 7.6, 6.8, 2.2, 1.0, 1.5, 2.2, 0.9, 3.2, 5.4, 6.3, 14.0, 12.7], sd_percent=[4.7, 5.2, 3.1, 3.9, 1.5, 1.1, 1.5, 1.6, 0.8, 2.3, 2.3, 3.0, 3.2, 3.8], n_people=28, locator='10.11138/ads/2017.8.2.079 Table3 printed page83, control group', total_force_N=None)
    return (h, f)

def lookup(fdi, pr):
    (h, f) = tables()
    k = pr['interval_multiplier_SD']
    if fdi // 10 in [3, 4]:
        i = fdi % 10 - 1
        if not 0 <= i < 7:
            raise ValueError('No Hattori no-M3 reference for this tooth')
        mean = h['mean_N'][i]
        sd = h['sd_N'][i]
        return dict(mean=mean, lo=max(0.0, mean - k * sd), hi=mean + k * sd, unit='N', resolution='POPULATION', study='Hattori1996', scope='Maximum voluntary clench without M3 contact; individual applicability UNKNOWN')
    if fdi not in f['fdi']:
        raise ValueError('No Ferrato reference for this upper FDI')
    i = f['fdi'].index(fdi)
    mean = f['mean_percent'][i]
    sd = f['sd_percent'][i]
    return dict(mean=mean, lo=max(0.0, mean - k * sd), hi=min(100.0, mean + k * sd), unit='percent of unmeasured arch force', resolution='POPULATION', study='Ferrato2017', scope='T-Scan relative signal, no absolute-force calibration in this geometry; table sum93% retained')

def predict():
    pr = json.loads((H / 'PREREG_R2.json').read_text())
    src = json.loads((X18 / 'rounds/R3.json').read_text())
    r1 = json.loads((H / 'rounds/R1.json').read_text())
    rows = []
    for r in r1['rows']:
        if r['status'] != 'SCORED':
            continue
        patches = r['contacts']['0.1']['reference']['patches']
        force = lookup(r['fdi'], pr)
        upper = force['hi'] if force['unit'] == 'N' else force['hi'] / 100
        pp = []
        for (j, p) in enumerate(patches):
            a = p['area_mm2']
            pp.append(dict(patch=j, resolution='PER_SURFACE_REGION', area_mm2=a, centroid_xyz_mm=p['centroid_xyz_mm'], average_pressure_lo=0.0, average_pressure_hi=upper / a, pressure_unit='MPa' if force['unit'] == 'N' else 'MPa per N of unknown arch force', physical_pressure_status='UNKNOWN; scan proximity is not pressure observation'))
        rows.append(dict(case=r['case'], fdi=r['fdi'], jaw=r['jaw'], tooth_force_scenario=force, patch_count=len(patches), patches=pp, force_fraction_constraints='nonnegative patch loads summing to tooth resultant; within-patch uniform pressure closure', joint_box_warning='Marginal patch upper bounds cannot occur simultaneously; use simplex', physical_pose_admissible=r['physical_pose_admissible']))
    stress = []
    for r in src['rows']:
        p = Path(r['basis_file'])
        z = np.load(p)
        basis = z['stress_tensors_MPa']
        areas = z['area_shares']
        force = lookup(r['fdi'], pr)
        assert force['unit'] == 'N'
        mean_scale = force['mean'] / 100
        hi_scale = force['hi'] / 100
        stress.append(dict(case=r['case'], fdi=r['fdi'], resolution='PER_TOOTH', tooth_force_scenario=force, old_area_peak_MPa=r['area_weighted_peak_MPa'], old_sharp_upper_MPa=r['sharp_upper_peak_MPa'], old_shape_ratio=r['upper_over_area_peak'], mean_area_peak_MPa=r['area_weighted_peak_MPa'] * mean_scale, mean_sharp_upper_MPa=r['sharp_upper_peak_MPa'] * mean_scale, scenario_peak_interval_MPa=[0.0, r['sharp_upper_peak_MPa'] * hi_scale], new_shape_ratio=r['upper_over_area_peak'], remaining_independent_shape_channels=len(areas) - 1, basis=artifact(p), site_provenance='Inherited X18 original roof/proximity patches, not new exact exterior specimen FE', physical_stress_status='UNKNOWN: no registered pressure/strain validation, support/material closures retained'))
    (h, f) = tables()
    files = [artifact(P / 'raw_sources/Hattori1996.pdf'), artifact(P / 'raw_sources/Hattori1996_page3.png'), artifact(P / 'raw_sources/Ferrato2017.pdf'), artifact(P / 'raw_sources/Ferrato2017_page5.png'), artifact(SOURCE), artifact(X18 / 'rounds/R3.json')] + [r['basis'] for r in stress]
    output = dict(round='R2', tables=dict(Hattori=h, Ferrato=f), pressure_predictions=rows, stress_predictions=stress)
    dump(H / 'raw/PREDICTIONS_R2.json', output)
    files.append(artifact(H / 'raw/PREDICTIONS_R2.json'))
    dump(H / 'FROZEN_PREDICTIONS_R2.json', dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R2.json'), files=files, code_sha256={str(p): sha(p) for p in [H / 'code/common.py', H / 'code/pressure.py']}, physical_pressure_measurement_status='NOT_RUN', new_patient_specific_force_data=False))
    state('R2_PREDICTIONS_FROZEN', 'Jaw-specific force and patch pressure intervals frozen', 'Verify source contracts, load conservation and scaling invariance')

def peak(b):
    t = np.zeros(b.shape[:-1] + (3, 3))
    t[..., 0, 0] = b[..., 0]
    t[..., 1, 1] = b[..., 1]
    t[..., 2, 2] = b[..., 2]
    t[..., 0, 1] = t[..., 1, 0] = b[..., 3]
    t[..., 1, 2] = t[..., 2, 1] = b[..., 4]
    t[..., 0, 2] = t[..., 2, 0] = b[..., 5]
    return float(max(0, np.linalg.eigvalsh(t)[..., -1].max()))

def evaluate():
    verify_freeze('R2')
    pr = json.loads((H / 'PREREG_R2.json').read_text())
    pred = json.loads((H / 'raw/PREDICTIONS_R2.json').read_text())
    start = time.perf_counter()
    checks = []
    stress = []
    ext = json.loads((P / 'EXTERNAL_REFERENTS.json').read_text())['studies']
    f = next((s for s in ext if s['id'] == 'Ferrato2017'))
    h = pred['tables']['Hattori']
    ft = pred['tables']['Ferrato']
    ferrato_parity = max((np.max(np.abs(np.asarray(ft[k]) - np.asarray(f[k]))) for k in ['mean_percent', 'sd_percent']))
    checks.append(dict(control='Independent Ferrato primary transcription versus earlier saved table', passed=ferrato_parity <= pr['metrics']['source_table_tolerance'], injected_wrong_value_rejected=abs(ft['mean_percent'][0] * 0.1 - f['mean_percent'][0]) > pr['metrics']['source_table_tolerance']))
    hf = next((s for s in ext if s['id'] == 'Hattori1996'))
    force_share = 100 * np.asarray(h['mean_N']) / 443.1
    fig_share = 2 * np.asarray(hf['single_side_mean_percent'])
    difference = np.abs(force_share - fig_share)
    checks.append(dict(control='Hattori Table1 force shares versus published Figure1 percentages (ratio of means need not equal mean ratio)', passed=None, max_difference_pp=float(difference.max()), tolerance_note='No admission gate from this descriptive cross-observable check', M1_M2_sum_difference_pp=float(abs(force_share[-2:].sum() - fig_share[-2:].sum()))))
    checks.append(dict(control='Hattori primary Table1 M2 transcription 211.8 N', passed=abs(h['mean_N'][-1] - 211.8) <= pr['metrics']['source_table_tolerance'], injected_wrong_value_rejected=abs(21.18 - 211.8) > pr['metrics']['source_table_tolerance']))
    checks.append(dict(control='Upper relative share cannot emit absolute force without arch resultant; jaw mismatch rejected', passed=ft['total_force_N'] is None and lookup(16, pr)['unit'] != 'N' and (lookup(36, pr)['unit'] == 'N'), injected_wrong_value_rejected=lookup(16, pr)['unit'] != 'N'))
    supports = json.loads(SOURCE.read_text())
    support_rows = []
    for s in supports:
        admissible = bool(s.get('verified')) and s.get('unit') == 'N/mm' and ('axial' in s.get('direction', '').lower()) and any((q in s.get('tooth', '').lower() for q in ['molar', 'premolar'])) and (s.get('min') is not None) and (s.get('max') is not None)
        support_rows.append(dict(id=s['id'], accepted=admissible, reason='posterior axial stiffness matches port' if admissible else 'wrong tooth/direction/quantity or UNKNOWN/unverified numerical bounds'))
    checks.append(dict(control='Support port specificity (lateral incisors are not axial posterior support)', passed=not any((s['accepted'] for s in support_rows)), injected_wrong_value_rejected=not next((s['accepted'] for s in support_rows if s['id'] == 'NAT-INC-U-LAT-K'))))
    for r in pred['stress_predictions']:
        z = np.load(r['basis']['path'])
        b = z['stress_tensors_MPa']
        area = z['area_shares']
        F = r['tooth_force_scenario']['mean']
        base = peak(np.einsum('j,jea->ea', area, b))
        old = max((peak(t) for t in b))
        new = max((peak(t * (F / 100)) for t in b))
        newbase = peak(np.einsum('j,jea->ea', area, b * (F / 100)))
        ratio = new / newbase
        err = abs(ratio - r['old_shape_ratio']) / r['old_shape_ratio']
        alpha = np.random.default_rng(r['case'] * 100 + r['fdi']).dirichlet(np.ones(len(b)))
        forces = F * alpha
        A = np.array([p['area_mm2'] for p in next((rr for rr in json.loads((X18 / 'raw/RESULTS_R1_ROWS.json').read_text()) if rr.get('status') == 'SCORED' and rr['case'] == r['case'] and (rr['fdi'] == r['fdi'])))['arms']['0.1']['original']['patches']])
        pressure = forces / A
        integrated = float(pressure @ A)
        wrong_integrated = float(2 * pressure @ A)
        stress.append(dict(**r, verified_new_shape_ratio=ratio, narrowing_factor=r['old_shape_ratio'] / ratio, scalar_information_shrink_percent=100 * (1 - ratio / r['old_shape_ratio']), parity_relative=err, narrowing_gate=r['old_shape_ratio'] / ratio >= pr['metrics']['narrowing_factor_min'], force_conservation_error_N=abs(integrated - F), injected2x_pressure_rejected=abs(wrong_integrated - F) > pr['metrics']['force_sum_tolerance_N']))
        checks.append(dict(control=f"Simplex force integral and stress scaling {r['case']}/{r['fdi']}", passed=err <= pr['metrics']['stress_ratio_parity_relative_max'] and abs(integrated - F) <= pr['metrics']['force_sum_tolerance_N'], injected_wrong_value_rejected=abs(wrong_integrated - F) > pr['metrics']['force_sum_tolerance_N']))
    out = dict(round='R2', claim_type='information_link', external_referent=pr['external_referent'], decision='NARROWING_PASS' if all((r['narrowing_gate'] for r in stress)) else 'NEGATIVE_PRESSURE_IDENTIFIABILITY', source_tables=pred['tables'], pressure_predictions=pred['pressure_predictions'], stress_results=stress, controls=checks, support_port=dict(required='Measured posterior axial stiffness N/mm at specified ramp/time', attempted=len(support_rows), retained=sum((s['accepted'] for s in support_rows)), rejected=sum((not s['accepted'] for s in support_rows)), rejection_fraction=sum((not s['accepted'] for s in support_rows)) / len(support_rows), rows=support_rows, physical_support_interval='UNKNOWN'), cost=dict(wall_s=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), physical_pressure_accuracy='UNKNOWN', physical_stress_accuracy='UNKNOWN', negative_result='Published tooth resultants do not constrain intra-tooth pressure fractions. Shape uncertainty is unchanged; absolute load scenarios may be wider than100N.', next_construction='Calibrating pressure force/moment ports while preserving unresolved null-space; no acquisition available locally')
    dump(H / 'rounds/R2.json', out)
    (H / 'HANDOFF_R2.md').write_text('R2 NEGATIVE_PRESSURE_IDENTIFIABILITY. Same nonnegative contact simplex after tooth-force information: original 1.896–14.057 shape ratio is unchanged (scalar scaling invariance). Absolute mandibular scenarios saved; upper pressure remains per N of unknown arch force. All posterior axial support candidates rejected. Next: explicitly design minimal registered spatial-force measurements and demonstrate which ambiguity moments still leave (new PREREG required).\n')
    state('R2_DECIDED', out['decision'], 'Freeze minimal registered force/moment calibration construction')
    print(json.dumps(dict(decision=out['decision'], shape_ratios=[r['verified_new_shape_ratio'] for r in stress], support_retained=out['support_port']['retained']), indent=2))
if __name__ == '__main__':
    {'predict': predict, 'evaluate': evaluate}[sys.argv[1]]()
