from dental_release.paths import expand as _release_expand
import os, sys, json, math, hashlib, datetime, csv, copy, time, subprocess
from pathlib import Path
from statistics import NormalDist
from analyze_r1 import R, read, dump, sha, state
from compare_lab import W, H, Z, PROTOCOLS, REQUIRED, calibration, compare, validate, integrity

def source_lock():
    for (p, d) in read('SOURCE_MANIFEST.json')['files'].items():
        if sha(p) != d['sha256']:
            raise ValueError('source manifest drift: ' + p)
    if (R / 'GEOMETRY_REFERENCES.json').exists():
        for f in read('GEOMETRY_REFERENCES.json')['files']:
            if hashlib.sha256(Path(f['path']).read_bytes()).hexdigest() != f['sha256']:
                raise ValueError('geometry reference drift: ' + f['path'])
    for p in R.glob('PREREG*.json'):
        if sha(str(p.relative_to(R))) != p.with_suffix('.sha256').read_text().strip():
            raise ValueError('prereg drift ' + str(p))

def geometry_refs():
    b = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X1B'))
    refs = []
    for d in ['D1', 'M1', 'M2']:
        for name in ['crown.stl', 'die.stl']:
            p = b / 'exports' / d / name
            if p.exists():
                refs.append(dict(design=d, path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size, license='STS-3D CC BY4.0; inherited X1b export, unchanged'))
    dump('GEOMETRY_REFERENCES.json', dict(geometry_status='READ_ONLY_REFERENCE; no new FE or modified export', files=refs))

def freeze():
    if (R / 'FROZEN_PREDICTIONS.json').exists():
        integrity()
        return
    geometry_refs()
    locked = {str(p.relative_to(R)): sha(str(p.relative_to(R))) for p in sorted([*R.glob('PREREG*.json'), *R.glob('code/*.py'), R / 'LAB_PROTOCOL.md', R / 'GEOMETRY_REFERENCES.json', R / 'SOURCE_MANIFEST.json', R / 'RUNTIME_MEASUREMENTS.json', R / 'run_all.sh'])}
    (lo, hi) = (math.exp(-H), math.exp(H))
    powers = []
    zpower = NormalDist().inv_cdf(0.8)
    effect = math.log(1.6)
    s = 1 + W ** 2 + (1 - W) ** 2
    for cv in [0.2, 0.35, 0.4]:
        sem = cv * math.sqrt(s / 12)
        minn = math.ceil((Z + zpower) ** 2 * cv ** 2 * s / (effect - H) ** 2)
        pow = 1 - NormalDist().cdf(Z - (effect - H) / sem)
        powers.append(dict(planning_CV=cv, n_per_group=12, probability_reject_upper_window_if_true_force_factor1_6=float(pow), n_per_group_for80percent_power_rejecting_window=minn, assumptions='normal approximation, independent group means, equalCV, no censoring; NOT measured laboratory scatter', resolution='PHENOMENOLOGICAL', replacement_measurement='actual endpoint SDs plus independent validation repeat'))
    designs = []
    for (prot, p) in PROTOCOLS.items():
        for (d, t) in [('D1', 0.8), ('M1', 1.0), ('M2', 1.5)]:
            designs.append(dict(design=d, protocol=prot, t_mm=t, role='BLIND_VALIDATION' if d == 'M1' else 'CALIBRATION_MEASUREMENT', prediction_formula='mu_D1^0.645016... * mu_M2^0.354983...' if d == 'M1' else None, actual_weights={'D1': 1 - W, 'M2': W} if d == 'M1' else None, operational_window_multipliers=[lo, hi] if d == 'M1' else None, absolute_force_N=None, validity='PROSPECTIVE_OPERATIONAL_HYPOTHESIS' if d == 'M1' else 'UNMEASURED_CALIBRATION_REQUIRED', resolution='PER_TOOTH', source_protocol=p))
    for d in designs:
        if d['role'] == 'BLIND_VALIDATION':
            d['prediction_formula'] = f'mu_D1^{1 - W:.16g} * mu_M2^{W:.16g}'
    f = dict(schema='x1c-frozen-predictions-v1', claim_type='information_link', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), measurement_status='NO_PHYSICAL_MEASUREMENTS', locked_files=locked, designs=designs, w=W, operational_halfwidth_log=H, operational_width_ratio=hi / lo, interval_kind='FIXED_OPERATIONAL_ERROR_WINDOW; NOT statistical population interval', new_information='two measured endpoint group means in each exact lab protocol;48 calibration crowns', blind_validation_crowns=24, power=powers, source_direct_ratio_hypotheses=read('raw/R1_RESULTS.json')['predictions'], old_X1b_Q={'value': 0.474375594239966, 'status': 'SUSPENDED; support, contact, speed and cement changed, so no isolated-angle claim', 'source_sha256': sha('inputs/x1b_frozen.json')}, no_lifetime_prediction=True, external_referent=dict(kind='independent_measurement', locator='doi:10.4047/jap.2021.13.5.269 Table1; doi:10.3390/ma17020365 Table2', compared_quantity='published interior-force residual, some Weibull-derived; future physical M1 outcome is absent', refutes_us=True), scientific_limit='No exact3Y protocol replication. R3 local-scope statistics are exploratory reanalysis of exposed R2 corpus.')
    dump('FROZEN_PREDICTIONS.json', f)
    (R / 'FROZEN_PREDICTIONS.sha256').write_text(sha('FROZEN_PREDICTIONS.json') + '\n')
    rows = []
    for (prot, p) in PROTOCOLS.items():
        for d in ['D1', 'M1', 'M2']:
            for i in range(1, 13):
                rows.append({**{k: '' for k in REQUIRED}, 'specimen_id': f'{prot}_{d}_{i:02}', 'protocol_id': prot, 'design_id': d, 'force_unit': 'N', 'load_angle_deg': p['angle'], 'die_E_MPa': p['E'], 'indenter_diameter_mm': p['diameter'], 'crosshead_mm_min': p['speed'], 'crown_product': p['crown'], 'cement_product': p['cement'], 'die_product': p['die'], 'surface_protocol': p['surface'], 'storage_days': 7, 'status': 'UNMEASURED'})
    with (R / 'LAB_MEASUREMENTS_TEMPLATE.csv').open('w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=REQUIRED)
        wr.writeheader()
        wr.writerows(rows)

def faults():
    vals = {}
    rr = []
    now = datetime.datetime.now(datetime.timezone.utc)
    for (prot, p) in PROTOCOLS.items():
        for d in ['D1', 'M2', 'M1']:
            base = {'D1': 1000.0, 'M2': 3000.0, 'M1': 1000 ** (1 - W) * 3000 ** W}[d]
            for i in range(12):
                r = {k: '' for k in REQUIRED}
                r.update(specimen_id=f'SYNTHETIC_{prot}_{d}_{i}', protocol_id=prot, design_id=d, force_unit='N', fracture_force_N=str(base * (0.9 + 0.2 * i / 11)), die_E_MPa=str(p['E']), indenter_diameter_mm=str(p['diameter']), crosshead_mm_min=str(p['speed']), load_angle_deg=str(p['angle']), material_batch='SYNTHETIC', crown_product=p['crown'], cement_product=p['cement'], cement_batch='SYNTHETIC', die_product=p['die'], die_batch='SYNTHETIC', surface_protocol=p['surface'], storage_days='7', interlayer_spec='SYNTHETIC_fixed_spec', manufacturing_spec_sha256='1' * 64, failure_mode='crown_tensile_fracture', fracture_origin='intaglio_tensile_zone', force_trace_sha256='2' * 64, measured_utc=(now + datetime.timedelta(seconds=60 if d == 'M1' else -60)).isoformat(), status='measured')
                rr.append(r)
    anchors = [r for r in rr if r['design_id'] != 'M1']
    cal = calibration(anchors)
    base = compare(rr, cal)
    vals['baseline_operational_test'] = all((p['force_verdict'] == 'SUPPORTED_OPERATIONAL_TEST' for p in base['protocols']))
    bad = copy.deepcopy(rr)
    for r in bad:
        if r['design_id'] == 'M1':
            r['fracture_force_N'] = str(float(r['fracture_force_N']) * 1.6)
    vals['force_x1_6_rejected'] = all((p['force_verdict'] == 'REJECTED' for p in compare(bad, cal)['protocols']))
    for (key, field, value) in [('wrong_unit', 'force_unit', 'kN'), ('wrong_angle', 'load_angle_deg', '17'), ('wrong_die_E', 'die_E_MPa', '18'), ('wrong_speed', 'crosshead_mm_min', '8'), ('wrong_indenter', 'indenter_diameter_mm', '5'), ('wrong_cement', 'cement_product', 'unknown'), ('wrong_crown_product', 'crown_product', '5Y unknown'), ('wrong_trace_hash', 'force_trace_sha256', '')]:
        bad = copy.deepcopy(rr)
        bad[0][field] = value
        try:
            compare(bad, cal)
            vals[key] = False
        except ValueError:
            vals[key] = True
    bad = copy.deepcopy(rr)
    bad[0]['fracture_force_N'] = '999'
    try:
        compare(bad, cal)
        vals['anchor_changed_after_freeze'] = False
    except ValueError:
        vals['anchor_changed_after_freeze'] = True
    bad = copy.deepcopy(rr)
    for r in bad:
        if r['design_id'] == 'M1':
            r['measured_utc'] = (now - datetime.timedelta(seconds=60)).isoformat()
    try:
        compare(bad, cal)
        vals['measurement_before_freeze'] = False
    except ValueError:
        vals['measurement_before_freeze'] = True
    bad = copy.deepcopy(rr)
    bad[0]['specimen_id'] = bad[1]['specimen_id']
    try:
        compare(bad, cal)
        vals['duplicate_specimen'] = False
    except ValueError:
        vals['duplicate_specimen'] = True
    bad = [r for r in rr if r['specimen_id'] != 'SYNTHETIC_C0_PROTT_M1_0']
    vals['missing_specimen_unknown'] = compare(bad, cal)['outcome'] == 'UNKNOWN'
    bad = copy.deepcopy(rr)
    for r in bad:
        if r['protocol_id'] == 'C0_PROTT' and r['design_id'] == 'M1':
            r['failure_mode'] = 'die_fracture'
            r['fracture_origin'] = 'die'
    vals['dominant_die_failure_rejected'] = any((q['status'] == 'FAIL_MECHANISM' for p in compare(bad, cal)['protocols'] for q in p['mechanism_gates']))
    bad = copy.deepcopy(rr)
    bad[-1]['status'] = 'censored'
    vals['censoring_unknown'] = compare(bad, cal)['outcome'] == 'UNKNOWN'
    import analyze_r1
    orig = analyze_r1.read
    primary = copy.deepcopy(read('inputs/primary_tables.json'))
    for r in primary:
        if r['pmcid'] == 'PMC10004144':
            for t in r['tables']:
                if t['label'] == 'Table 2':
                    t['text'] = t['text'].replace('14.65 ± 2.25', '146.5 ± 2.25')
    analyze_r1.read = lambda p: primary if p == 'inputs/primary_tables.json' else orig(p)
    try:
        analyze_r1.curated()
        vals['source_mean_x10_rejected'] = False
    except AssertionError:
        vals['source_mean_x10_rejected'] = True
    finally:
        analyze_r1.read = orig
    vals['wrong_n1_external_gate'] = read('raw/R1_RESULTS.json')['controls']['wrong_n1_rejected']
    vals['wrong_equal_information_control_x1_6_rejected'] = math.log(1.6) > H
    dump('raw/FAULT_INJECTIONS.json', dict(external_referent=dict(kind='our_own_fixture', locator='code/run_demo.py faults()', compared_quantity='software gates only, not physical prediction accuracy', refutes_us=False), checks=vals, all_pass=all(vals.values()), force_fault_output=compare([dict(r, fracture_force_N=str(float(r['fracture_force_N']) * 1.6)) if r['design_id'] == 'M1' else r for r in rr], cal)))
    if not all(vals.values()):
        raise ValueError('fault injection failed ' + str(vals))
    return vals

def figure():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = read('raw/CURATED_ROWS.json')
    (fig, axs) = plt.subplots(1, 2, figsize=(10.5, 4.3))
    cols = ['#156e8c', '#a7591e', '#824596']
    free = [r for r in rows if r['protocol'] == 'free_disc_3ball_0.5']
    for (prod, c) in zip(['ESS', 'EMX', 'LP'], cols):
        rr = sorted([r for r in free if r['product'] == prod], key=lambda r: r['t'])
        ref = rr[0]
        axs[0].plot([r['t'] for r in rr], [r['mean'] / ref['mean'] for r in rr], 'o', color=c, label=prod + ' measured')
    ts = [0.4 + i * 0.012 for i in range(101)]
    axs[0].plot(ts, [(t / 0.4) ** 2 for t in ts], '-', color='black', label='n=2 free-disc law')
    axs[0].plot(ts, [(t / 0.4) ** 1.398 for t in ts], '--', color='gray', label='inherited n=1.398')
    axs[0].set(xlabel='Disc thickness [mm]', ylabel='Mean force / 0.4mm reference', title='External Table 2: free discs only')
    axs[0].legend(fontsize=8)
    prott = sorted([r for r in rows if r['protocol'] == 'crown_0_Prott_resin' and 0.8 <= r['t'] <= 1.5], key=lambda r: r['t'])
    (a, b) = (prott[0], prott[-1])
    xx = [0.8 + i * 0.007 for i in range(101)]
    ff = [a['mean'] ** (1 - math.log(t / 0.8) / math.log(1.5 / 0.8)) * b['mean'] ** (math.log(t / 0.8) / math.log(1.5 / 0.8)) for t in xx]
    axs[1].plot(xx, ff, color=cols[0], label='endpoint-calibrated shape')
    axs[1].fill_between(xx, [v * math.exp(-H) for v in ff], [v * math.exp(H) for v in ff], alpha=0.2, color=cols[0], label='operational +/-0.20 log')
    axs[1].plot([r['t'] for r in prott], [r['mean'] for r in prott], 'ko', label='published Weibull-derived means')
    axs[1].set(xlabel='Crown nominal thickness [mm]', ylabel='Arithmetic mean force [N]', title='3Y axial: one study; exploratory')
    axs[1].legend(fontsize=8)
    fig.suptitle('X1c: separate load protocols, then calibrate locally', fontsize=13)
    fig.tight_layout()
    fig.savefig(R / 'crown_literature.png', dpi=160)
    fig.savefig(R / 'crown_literature.pdf', metadata={'CreationDate': datetime.datetime(2026, 10, 3, tzinfo=datetime.timezone.utc)})
    plt.close(fig)

def results(vals):
    r1 = read('raw/R1_RESULTS.json')
    r2 = read('raw/R2_RESULTS.json')
    r3 = read('raw/R3_RESULTS.json')
    freeze = read('FROZEN_PREDICTIONS.json')
    out = dict(schema='x1c-results-v1', round_tag='X1c-crown-literature', claim_type='information_link', status='PENDING_INDEPENDENT_REVIEW', capability='Run protocol-separated source checks and explicit leave-study-out abstention, then freeze two local endpoint-conditioned M1 force windows before physical testing', outcome='NEW_FALSIFIABLE_LOCAL_CALIBRATION; ABSOLUTE_FORCE_AND_EXACT3Y_TRANSFER_UNKNOWN', external_referent=dict(kind='independent_measurement', locator='doi:10.3390/ma16051997 Table2; doi:10.3390/ma16176006 Table3; doi:10.4047/jap.2021.13.5.269 Table1; doi:10.3390/ma17020365 Table2', compared_quantity='published mean-force thickness ratios / held-out interior means conditional on two endpoint measurements; Prott arithmetic means derived from published Weibull parameters', refutes_us=True), quantity_resolution={'external_force_means': 'PER_TOOTH', 'summary_metrics': 'POPULATION', 'free_disc_force_ratio': 'PHENOMENOLOGICAL', 'operational_tolerance': 'PHENOMENOLOGICAL', 'planning_power': 'PHENOMENOLOGICAL'}, phenomenological_debts=freeze['power'], rounds={'R1': dict(outcome=r1['outcome'], gates=r1['gates'], free_disc_logRMSE=r1['free_disc_logRMSE'], LOSO_abstention=r1['LOSO_abstention']), 'R2': dict(gates=r2['gates'], stats=r2['stats'], outcome='FAILED_FULL_RANGE_FIXED_WINDOW'), 'R3': dict(gates=r3['gates'], stats=r3['stats'], outcome='LOCAL_OPERATIONAL_FIT_ONLY; exact3Y replication UNKNOWN', scope='reanalysis of exposed R2 corpus; not independent external confirmation'), 'R4': dict(outcome='PROSPECTIVE_FROZEN_PROTOCOL_AND_SCORER; NO PHYSICAL MEASUREMENTS', operational_width_ratio=freeze['operational_width_ratio'], calibration_crowns=48, validation_crowns=24, total_crowns=72)}, matched_material_port=read('raw/MATCHED_MATERIAL_PORT.json'), same_information_control=dict(R1_GLS_max_slope_difference=r1['controls']['block_GLS_max_slope_difference'], R3_affine_interpolation_max_force_difference_N=r3['largest_control_difference_N'], outcome='Equivalent standard control; this is an information link, no algorithmic superiority claimed'), source_audit=read('SOURCE_AUDIT.json'), frozen_predictions_sha256=sha('FROZEN_PREDICTIONS.json'), lab_protocol_sha256=sha('LAB_PROTOCOL.md'), fault_injections=dict(checks=len(vals), all_rejecting_gates_pass=all(vals.values()), source='raw/FAULT_INJECTIONS.json', kind='our_own_fixture; software only'), geometry_references=read('GEOMETRY_REFERENCES.json'), edges=[dict(source='published primary force cell: study/product/protocol/thickness', consumer='DENT-DESIGN-CROWN-PROBLEM conditional force response', quantity='group mean force N and thickness mm', resolution_level='PER_TOOTH', timescale='SIMULTANEOUS', status='PENDING_INDEPENDENT_REVIEW', numeric_fusion=False), dict(source='physical endpoint force summaries, not yet measured', consumer='M1 premeasurement freeze', quantity='protocol-specific endpoint means and local secant slope', resolution_level='PER_TOOTH', timescale='HANDOVER', status='PROPOSED_UNMEASURED', numeric_fusion=False)], cost={'preparation_bytes': sum((p.stat().st_size for p in (R / 'inputs').iterdir())), 'fit_seconds_R1': read('RUNTIME_MEASUREMENTS.json')['R1_seconds'], 'validation_seconds_R2': read('RUNTIME_MEASUREMENTS.json')['R2_seconds'], 'validation_seconds_R3': read('RUNTIME_MEASUREMENTS.json')['R3_seconds'], 'CPU_threads': 1, 'GPU': False, 'fresh_FE_solves': 0, 'lab_calibration_crowns': 48, 'lab_validation_crowns': 24, 'actual_lab_cost': 'UNKNOWN: not performed', 'human_questions': 0, 'research_human_and_token_cost': 'UNKNOWN; not silently zero', 'fallback': 'No numeric new N intervals until endpoint measurements; no overwrite or target refit'}, raw_files={str(p.relative_to(R)): dict(sha256=sha(str(p.relative_to(R))), bytes=p.stat().st_size) for p in sorted((R / 'raw').glob('*.json')) if p.name != 'GRAPH_PACKET.json'}, next_construction='Independent exact-product3Y same-geometry 0.8/1.0/1.5 crown measurements under each fixed protocol; if residual fails, couple registered contact and cement failure modes before fitting strength')
    dump('results.json', out)
    dump('ATTEMPTS.json', [dict(round='R1', changed_operation='exact protocol keys + source-corrected rows', outcome=r1['outcome'], next='two endpoint local calibration'), dict(round='R2', changed_operation='full source-span endpoint calibration', outcome='FAIL RMSE/coverage', next='actual0.8..1.5 lab span'), dict(round='R3', changed_operation='actual lab span, same tolerance', outcome='numeric gates pass; exact3Y independent replication missing; exploratory', next='two-stage physical prediction freeze'), dict(round='R4', changed_operation='immutable physical endpoints -> blind target with matched protocols', outcome='software rejection gates pass; physical outcome UNKNOWN', next=out['next_construction'])])
    return out

def material_port():
    import xml.etree.ElementTree as ET, re
    root = ET.parse(R / 'inputs/PMC10817558.xml').getroot()
    cells = []
    mat = thick = None
    for wrap in root.iter('table-wrap'):
        if ''.join(wrap.find('label').itertext()) != 'Table 2':
            continue
        for row in wrap.iter('tr'):
            cs = [''.join(c.itertext()).strip() for c in row]
            if cs[0] == 'Zirconia':
                continue
            if len(cs) == 5:
                (mat, thick, abrasion, cement, force) = cs
            elif len(cs) == 4:
                (thick, abrasion, cement, force) = cs
            elif len(cs) == 3:
                (abrasion, cement, force) = cs
            else:
                raise ValueError('unexpected source row span')
            numbers = re.findall('\\d+(?:\\.\\d+)?', force)
            cells.append(dict(material=mat, t=float(thick), abrasion=abrasion, cement=cement, mean=float(numbers[0]), sd=float(numbers[1])))
    a = next((c for c in cells if (c['material'], c['t'], c['abrasion'], c['cement']) == ('3Y-Z', 1.0, 'Yes', 'RMGI')))
    b = next((c for c in cells if (c['material'], c['t'], c['abrasion'], c['cement']) == ('5Y-Z', 1.0, 'Yes', 'RMGI')))
    assert a['mean'] == 1378.1 and a['sd'] == 143.22 and (b['mean'] == 663.78) and (b['sd'] == 106.8)
    table1 = ' '.join((t for w in root.iter('table-wrap') if ''.join(w.find('label').itertext()) == 'Table 1' for t in w.itertext()))
    assert 'Katana HT' in table1 and 'Katana UTML' in table1
    e = math.log(a['mean'] / b['mean'])
    se = math.hypot(a['sd'] / a['mean'], b['sd'] / b['mean']) / math.sqrt(10)
    ratio = math.exp(e)
    ci = [math.exp(e - 1.96 * se), math.exp(e + 1.96 * se)]
    expected = {'Katana HT': a['mean'], 'Katana UTML': b['mean']}

    def bind(product, value):
        if expected[product] != value:
            raise ValueError('material/force source association mismatch')
    try:
        bind('Katana HT', b['mean'])
        rejected = False
    except ValueError:
        rejected = True
    out = dict(claim_type='information_link', external_referent=dict(kind='independent_measurement', locator='doi:10.3390/ma17020365 Tables1–2', compared_quantity='matched1mm mean fracture-force ratio Katana HT3Y/UTML5Y, sameRMGIC/abrasion,30deg,NextDent support', refutes_us=False), design='M1 nominal1.0mm', protocol='C30_CHEN', products=['Katana HT3Y', 'Katana UTML5Y'], ratio=ratio, measurement_95CI=ci, n_per_product=10, source_study_rank='HT3Y > UTML5Y' if ci[0] > 1 else 'UNKNOWN', resolution='PER_TOOTH', material_swap_rejected=rejected, transport_to_X1c_STS_molar='UNKNOWN: primary is premolar uniform coping; no5Ylabgroup in current72', other_design_protocol_ranks='UNKNOWN: no same-product matched-thickness primarycontrast for D1/M2/C0')
    dump('raw/MATCHED_MATERIAL_PORT.json', out)
    return out

def main():
    t0 = time.perf_counter()
    source_lock()
    if (R / 'FROZEN_PREDICTIONS.json').exists():
        integrity()
    for command in [['code/analyze_r1.py'], ['code/analyze_r2.py'], ['code/analyze_r2.py', '--local']]:
        subprocess.run([sys.executable, '-s', *command], cwd=R, check=True, stdout=subprocess.DEVNULL)
    freeze()
    vals = faults()
    material = material_port()
    vals['material_swap_rejected'] = material['material_swap_rejected']
    figure()
    out = results(vals)
    integrity()
    seconds = time.perf_counter() - t0
    dump('VERIFICATION.json', dict(status='PASS', frozen_sha256=sha('FROZEN_PREDICTIONS.json'), source_manifest_verified=True, prereg_verified=True, fault_checks=len(vals), all_fault_checks_pass=all(vals.values()), elapsed_seconds=seconds, physical_measurements=0, scientific_status='PENDING_INDEPENDENT_REVIEW', runtime=dict(python=sys.version, executable=sys.executable)))
    state('DELIVERABLE_VERIFIED_PENDING_PHYSICAL_TEST', dict(software=True, R1_crown_transfer=False, R2_full_range=False, R3_exact3Y_replication=False, physical_test='UNKNOWN'), 'Run 48 calibration crowns, freeze numeric M1 windows, then24 blind validation crowns; independent protocol replication and contact/cement branch are next')
    print(json.dumps(dict(status='PASS_SOFTWARE_ONLY', physical_measurements=0, fault_checks=len(vals), frozen_sha256=sha('FROZEN_PREDICTIONS.json'), elapsed_seconds=seconds), indent=2))
if __name__ == '__main__':
    main()
