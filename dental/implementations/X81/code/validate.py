"""Independent data/decision controls and actual targeted value mutations."""
import copy, hashlib, json, math, re
import numpy as np
from lxml import etree
from scipy.stats import t
from common import ROOT, read, save, verify_freeze, scalar_decision
from run_r3 import ready_edges, close_prerequisites, added
from run_r2 import welch_hw, interp_control

def weighted_median(values, counts):
    order = np.argsort(values)
    v = np.asarray(values)[order]
    c = np.cumsum(np.asarray(counts)[order], dtype=np.int64)
    n = int(c[-1])
    at = lambda k: float(v[np.searchsorted(c, k + 1, side='left')])
    return (at((n - 1) // 2) + at(n // 2)) / 2

def compare_map(a, b, tol=1e-10):
    return set(a) == set(b) and all((abs(a[k] - b[k]) <= tol for k in a))

def replica_source_map():
    text = (ROOT / 'raw/input_replica_text.txt').read_text()
    section = text.split('Table 1. Means and standard deviation')[1].split('In each column')[0]
    source = []
    for line in section.splitlines():
        pairs = re.findall('(\\d+\\.\\d+)±(\\d+\\.\\d+)', line)
        if pairs:
            source.append([(float(a), float(b)) for (a, b) in pairs])
    assert len(source) == 4 and all((len(x) == 4 for x in source))
    regions = ['marginal', 'axial', 'transition', 'occlusal']
    systems = ['Amann', 'Dentsply_Sirona']
    out = {}
    for (sy, idx) in zip(systems, [0, 2]):
        for (j, reg) in enumerate(regions):
            for (method, i) in [('replica', idx), ('ct', idx + 1)]:
                out[sy, reg, method, 'mean'] = source[i][j][0]
                out[sy, reg, method, 'sd'] = source[i][j][1]
    return out

def replica_snapshot_map(rows):
    return {(r['system'], r['region'], m, kind): float(r[f'{m}_{kind}_um']) for r in rows for m in ['replica', 'ct'] for kind in ['mean', 'sd']}

def paired_current_preload(p):
    d = p['first']['mean_N'] - p['tenth']['mean_N']
    hw = t.ppf(0.975, p['n'] - 1) * (p['first']['sd_N'] + p['tenth']['sd_N']) / math.sqrt(p['n'])
    return [float(d - hw), float(d + hw)]

def coverage_check(out):
    n = read('net')
    numeric = {e['id'] for e in n['edges'] if e.get('gap', {}).get('value') is not None}
    c = read('contracts')
    phy = {e for (e, v) in c['edges'].items() if v.get('status') == 'BLOCKED_ON_PHYSICAL_MEASUREMENT'}
    groups = {g for e in phy for g in c['edges'][e].get('measurement_ids', [])}
    return numeric == {r['edge_id'] for r in out['numeric_gap_rows']} and phy == set(out['physical_edge_obligations']) and (groups == {r['id'] for r in out['physical_group_rows']})

def raw_baseline_check():
    audits = []
    sources = []
    for rec in read('cbct_baseline') + read('cbct_held'):
        meta = rec['raw_samples']
        p = __import__('pathlib').Path(meta['path'])
        b = p.read_bytes()
        assert hashlib.sha256(b).hexdigest() == meta['sha256']
        z = np.load(p)
        numeric = {}
        for row in rec['inserts']:
            k = row['material']
            m = weighted_median(z[k + '_gray'], z[k + '_legacy_multiplicity'])
            numeric[k] = m
            assert abs(m - row['gray_median']) <= 1e-10
        audits.append({'device': rec['device'], 'scan': rec.get('zip_member', rec['scan']), 'weighted_medians': numeric, 'max_error': max((abs(numeric[x['material']] - x['gray_median']) for x in rec['inserts'])), 'n_ROIs': len(numeric), 'resolution': 'PER_SURFACE_REGION', 'median_is_not_SE': True})
        sources.append({'path': str(p), 'sha256': meta['sha256'], 'bytes': len(b), 'role': 'Raw independently published image ROI values, inherited verified archive-member mapping'})
    return (audits, sources)

def source_preload_check(p):
    root = etree.parse(str(ROOT / 'raw/input_preload_xml.xml'))
    text = ' '.join(root.xpath('//p//text()'))
    matches = re.search('329\\.9\\s*N\\s*±\\s*33\\.3.*?253\\.7\\s*N\\s*±\\s*36\\.8', text)
    expected = [329.9, 33.3, 253.7, 36.8]
    actual = [p['first']['mean_N'], p['first']['sd_N'], p['tenth']['mean_N'], p['tenth']['sd_N']]
    table = root.xpath('//table-wrap[@id="Tab1"]//tbody//tr')
    rows = []
    for tr in table:
        cells = [''.join(x.itertext()).strip() for x in tr.findall('td')]
        if len(cells) == 6:
            rows.append([float(x) for x in cells[:5]])
    return bool(matches) and actual == expected and (rows == [[r['from'], r['to'], r['difference_N'], r['lo_N'], r['hi_N']] for r in p['rows']])

def source_force_map():
    root = etree.parse(str(added('chen_xml')))
    text = ' '.join(root.xpath('//table-wrap[@id="materials-17-00365-t002"]//text()'))
    pairs = re.findall('(\\d+\\.\\d+)\\s*±\\s*(\\d+\\.\\d+)', text)
    return {float(a): float(b) for (a, b) in pairs}

def main():
    for f in ['INPUT_LOCK.json', 'INPUT_ADDITION_R3.json', 'PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json', 'DECOMPOSITION_R1.json', 'DECOMPOSITION_R2.json', 'DECOMPOSITION_R3.json', 'FROZEN_PREDICTIONS_R1.json', 'FROZEN_PREDICTIONS_R2.json', 'FROZEN_PREDICTIONS.json']:
        verify_freeze(f)
    r1 = json.loads((ROOT / 'rounds/R1_results.json').read_text())
    r2 = json.loads((ROOT / 'rounds/R2_results.json').read_text())
    r3 = json.loads((ROOT / 'rounds/R3_results.json').read_text())
    (raw, sources) = raw_baseline_check()
    faults = []

    def fault(name, rejected, details=''):
        faults.append({'fault': name, 'rejected': bool(rejected), 'details': details})
        assert rejected, name
    smap = replica_source_map()
    rep = read('replica')
    assert compare_map(smap, replica_snapshot_map(rep))
    wrong = copy.deepcopy(rep)
    wrong[3]['replica_mean_um'] = '197.35'
    fault('replica_primary_table_plus10um', not compare_map(smap, replica_snapshot_map(wrong)))
    p = read('preload')
    assert source_preload_check(p)
    wrong = copy.deepcopy(p)
    wrong['tenth']['mean_N'] += 10
    fault('preload_primary_response_plus10N', not source_preload_check(wrong))
    fmap = source_force_map()
    emp = [q for q in r3['force_intervals'] if q['source_origin'] == 'MEASURED_GROUP_MEAN']
    assert all((abs(fmap[q['observed_N']] - q['sample_SD_N']) <= 1e-10 for q in emp))
    q = copy.deepcopy(emp[0])
    q['sample_SD_N'] += 10
    fault('fracture_source_SD_plus10N', abs(fmap[q['observed_N']] - q['sample_SD_N']) > 1e-10)
    bad_baseline = read('cbct_baseline')
    bad_baseline[0]['inserts'][0]['gray_median'] += 1
    fault('ROI_median_plus1gray', abs(raw[0]['weighted_medians']['AIR'] - bad_baseline[0]['inserts'][0]['gray_median']) > 1e-10)
    original = (ROOT / 'raw/input_cbct_baseline.json').read_bytes()
    expected = json.loads((ROOT / 'INPUT_LOCK.json').read_text())['cbct_baseline']['sha256']
    fault('input_bytes_appended', hashlib.sha256(original + b' ').hexdigest() != expected)
    original = (ROOT / 'PREREG_R2.json').read_bytes()
    wrong = original.replace(b'"cbct_gate_HU_ref": 40', b'"cbct_gate_HU_ref": 1000')
    assert wrong != original
    fault('frozen_gate_40_to1000', hashlib.sha256(wrong).hexdigest() != (ROOT / 'PREREG_R2.json.sha256').read_text().strip())
    assert coverage_check(r3)
    wrong = copy.deepcopy(r3)
    wrong['numeric_gap_rows'].pop()
    fault('omit_numeric_edge', not coverage_check(wrong))
    wrong = copy.deepcopy(r3)
    wrong['physical_group_rows'].pop()
    fault('omit_measurement_group', not coverage_check(wrong))
    groups = read('contracts')['groups']
    phy = r3['physical_edge_obligations']
    pilot = r3['lab_frontier']['pilot']
    assert ready_edges(pilot['ids'], phy) == pilot['type_ready_edges']
    wrong_ids = [q for q in pilot['ids'] if q != 'M01']
    fault('package_drops_asbuilt_prerequisite', ready_edges(wrong_ids, phy) != pilot['type_ready_edges'])
    force = next((q for q in emp if q['source_mean_gate'] == 'PASS'))
    (lo, hi) = force['simultaneous95_response_mean_interval_N']
    (al, ah) = force['fixed_prediction_acceptance_N']
    fault('accepted_force_response_shift_plus1000N', not (al <= lo + 1000 and hi + 1000 <= ah))
    laser = r2['thermal']['rows']
    bad = next((q for q in laser if q['minimum_drift_bound_0p4'].get('reason') == 'BIAS_EXCEEDS_CONTRAST'))
    wrong_n = 10 ** 12
    wrong_projected_lower = bad['dry_minus_cooled_degC'] - float(welch_hw(bad['sample_SD_dry_degC'], bad['sample_SD_cooled_degC'], wrong_n, wrong_n, 9)) - 0.4 - 0.01
    fault('infinite_repeats_claim_for_drift_limited_contrast', wrong_projected_lower <= 0, 'Actual injected n=10^12 still fails; differential bias not averaged')
    plan = next((q for q in laser if q['minimum_no_drift']['total_specimens'] > 20))
    sd1 = plan['sample_SD_dry_degC']
    sd2 = plan['sample_SD_cooled_degC']
    fault('replace_thermal_allocation_by2plus2', float(welch_hw(sd1, sd2, 2, 2, 9)) + 0.01 >= plan['dry_minus_cooled_degC'])
    minpair = r2['paired']['preload_minimum']
    assert minpair['n'] == 6
    hw = float(t.ppf(0.975, 4) * (p['first']['sd_N'] + p['tenth']['sd_N']) / math.sqrt(5) + 0.1)
    fault('preload_minimum6_replaced_by5', p['first']['mean_N'] - p['tenth']['mean_N'] <= hw)
    best = r2['calibration']['selected']
    base = read('cbct_baseline')
    ref = next((q for q in base if q['device'] == 'Varian'))
    h = np.array([z['gray_median'] for z in ref['inserts']])
    sel = best['indices']
    maxdiff = 0
    for q in base:
        if q['device'] == 'Varian':
            continue
        g = np.array([z['gray_median'] for z in q['inserts']])
        fit = best['fits'][q['device']]
        x = (g - fit['center']) / fit['scale']
        ctrl = interp_control(x[sel], h[sel], x)
        ours = np.polynomial.polynomial.polyval(x, fit['coefficients'])
        maxdiff = max(maxdiff, float(abs(ctrl - ours).max()))
        wrong = fit['coefficients'].copy()
        wrong[0] += 10
        fault('calibration_' + q['device'] + '_intercept_plus10HU', float(abs(ctrl - np.polynomial.polynomial.polyval(x, wrong)).max()) > 1e-08)
    assert maxdiff <= 1e-08

    def physical_claim_ok(x):
        return x['full_physical_decisions_certified'] == 0 and all((q['full_physical_gate'] == 'UNKNOWN' for q in x['force_intervals']))
    assert physical_claim_ok(r3)
    wrong = copy.deepcopy(r3)
    wrong['full_physical_decisions_certified'] = 1
    wrong['force_intervals'][0]['full_physical_gate'] = 'PASS'
    fault('promote_missing_bound_to_physical_pass', not physical_claim_ok(wrong))
    for s in [r1['sufficiency'], r3['paired_sufficiency']]:
        assert s['identity_error'] == 0
    assert r1['sufficiency']['downstream_decision_difference_count'] == 2
    assert r3['paired_sufficiency']['decisions'][0] != r3['paired_sufficiency']['decisions'][1]
    assert scalar_decision(-1, 1, 0) == 'UNRESOLVED'
    fault('report_straddling_interval_as_PASS', scalar_decision(-1, 1, 0) != 'PASS')
    out = {'status': 'PASS', 'all_faults_rejected': all((q['rejected'] for q in faults)), 'faults': faults, 'raw_cbct_ROI_audit': raw, 'raw_array_sources': sources, 'raw_CBCT_scans_checked': len(raw), 'raw_ROI_medians_checked': sum((q['n_ROIs'] for q in raw)), 'replica_numeric_source_values': len(smap), 'preload_adjusted_source_contrasts': len(p['rows']), 'fracture_empirical_means_checked': len(emp), 'calibration_control_max_difference': maxdiff, 'parity_control': 'Independent primary table parsing, raw weighted medians, explicit dependency loops and Lagrange interpolation', 'correctness_does_not_validate_physical_models': True, 'rigorous_physical_enclosure': 'MISSING', 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    save('VALIDATION.json', out)
    print('validation', out['status'], len(faults), 'injected faults rejected;', out['raw_ROI_medians_checked'], 'ROI medians verified')
if __name__ == '__main__':
    main()
