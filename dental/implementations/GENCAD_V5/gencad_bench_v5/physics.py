"""Fail-closed physical ports. Reuse primary rows and corrected v4 observation gates."""
from common import *
import csv, gzip, math, importlib.util, collections
spec = importlib.util.spec_from_file_location('reviewed_v4_controls', ROOT.parent / 'LANE_XREVIEW_BATCH10/patches/PROOF_LANE_GENCAD_V4/review_controls.py')
ctrl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctrl)
PROTOCOL = ['geometry_sha256', 'material', 'process_batch', 'support', 'contact', 'cement', 'aging', 'quantity', 'units', 'resolution']

def calibrate_candidate(candidate, calibration):
    missing = [k for k in PROTOCOL if candidate.get(k) in [None, 'UNKNOWN', '']]
    if not calibration:
        return dict(status='UNKNOWN_NO_MATCHED_MEASUREMENT', force05_N=None, force05_interval_N=None, cementfilm_um=None, cementfilm_interval_um=None, missing=missing + ['matched external observation'], resolution='PER_TOOTH')
    mismatch = [k for k in PROTOCOL if candidate.get(k) != calibration.get(k)]
    if missing or mismatch:
        return dict(status='UNKNOWN_PROTOCOL_MISMATCH', value=None, missing=missing, mismatch=mismatch)
    if candidate['quantity'] == 'fracture_force05':
        (F, m) = (calibration.get('characteristic_load_N'), calibration.get('Weibull_m'))
        if not F or not m or (not np.isfinite([F, m]).all()) or (min(F, m) <= 0):
            return dict(status='UNKNOWN_PARAMETER')
        return dict(status='MATCHED_MODEL_QUANTILE', value=F * (-math.log(0.95)) ** (1 / m), units='N', uncertainty='Parameter uncertainty must be supplied separately; no empirical geometry transfer')
    if candidate['quantity'] == 'seated_film':
        if calibration.get('measurement_state') != 'SEATED_CEMENTED':
            return dict(status='UNKNOWN_WRONG_OBSERVABLE', value=None)
        value = calibration.get('observed_um')
        if value is None or not np.isfinite(value) or value < 0:
            return dict(status='UNKNOWN_MEASUREMENT')
        return dict(status='MATCHED_OBSERVATION', value=value, units='um', uncertainty=calibration.get('uncertainty', 'UNKNOWN'))
    return dict(status='UNKNOWN_WRONG_OBSERVABLE', value=None)

def quantile(F, m):
    return F * (-math.log(0.95)) ** (1 / m)

def run():
    start = time.perf_counter()
    p = read(V4 / 'FROZEN_LITERATURE_PREDICTIONS.json')
    observations = []
    tests = []
    for r in p['rows']:
        observations.append(dict(r, quantity='dry_marginal_gap' if r['kind'] == 'cement' else 'crown_group_mean_force', resolution='POPULATION', source_spatial_support='PER_SURFACE_REGION' if r['kind'] == 'cement' else 'PER_TOOTH', gate=ctrl.observed_gate(r)))
        tests.append(dict(id=r.get('row_id', r['study'] + '_' + str(r['x'])), fault_rejected=ctrl.fault_probe(r)))
    cement = list(csv.DictReader((V4 / 'payload/literature/CEMENT.csv').open()))
    cm = {r['row_id']: r for r in cement}
    curated = read(V4 / 'payload/literature/CROWN_CURATED.json')
    direct = []
    for r in p['rows']:
        if r['kind'] == 'cement':
            anchors = [cm[i] for i in r['calibration_rows']]
            X = np.array([[1, 1 / float(a['internal_spacer_um'])] for a in anchors])
            y = np.array([float(a['measured_mean_um']) for a in anchors])
            value = float(np.array([1, 1 / r['x']]) @ np.linalg.solve(X, y))
        else:
            anchors = [a for a in curated if a['study'] == r['study'] and a['protocol'] == r['protocol'] and (a['product'] == r['product']) and (a['t'] in r['calibration_thicknesses_mm'])]
            anchors = sorted(anchors, key=lambda a: a['t'])
            X = np.array([[1, math.log(a['t'])] for a in anchors])
            y = np.log([a['mean'] for a in anchors])
            value = math.exp(float(np.array([1, math.log(r['x'])]) @ np.linalg.solve(X, y)))
        direct.append(dict(kind=r['kind'], x=r['x'], predicted=value, stored_candidate=r['predicted'], absolute_difference=abs(value - r['predicted']), same_information=True))
    groups = read(ROOT.parent / 'LANE_X1B_CROWN_LOOP/raw/MATCHED_R2.json')['groups']
    quant = []
    for r in groups:
        if r['study'] != 'PMC8558575':
            continue
        (F, m) = (r['characteristic_load_N'], r['Weibull_m'])
        q = quantile(F, m)
        corners = [quantile(f, mm) for f in r['F0_90CI_N'] for mm in r['m_90CI']]
        quant.append(dict(thickness_mm=r['thickness_mm'], force05_N=q, parameter_rectangle_range_N=[min(corners), max(corners)], resolution='POPULATION', interval_scope='Analytic extrema of monotone positive-parameter rectangle; no directed-rounding enclosure. Marginal 90CI rectangle is NOT a joint 90CI', locator=r['locator'], external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.4047/jap.2021.13.5.269', compared_quantity='published fitted characteristic force and Weibull m', refutes_us=True)))
    x60 = read(ROOT.parent / 'LANE_X60_CROWN_OPTIMIZER/raw/SCORED_ROWS.json')
    extras = [r for r in x60 if r['participant'] in ['volume_1.0', 'volume_1.1']]
    dump(ROOT / 'raw/X60_PARTICIPANT_ROWS.json', extras)
    rows = []
    import item_analysis
    for r in item_analysis.load_rows():
        if r['status'] != 'DESIGN':
            continue
        port = calibrate_candidate({'material': 'KATANA_HTML_PLUS', 'geometry_sha256': None, 'quantity': 'fracture_force05', 'units': 'N', 'resolution': 'PER_TOOTH'}, None)
        rows.append(dict(task_id=r['task_id'], participant=r['participant'], track='v4_roof', force05_N=port['force05_N'], force05_interval_N=None, cementfilm_um=None, cementfilm_interval_um=None, physical_status=port['status'], nominal_CAD_gap='geometry only; dry fit and film transfer not established', resolution='PER_TOOTH'))
    for r in extras:
        if r['status'] != 'DESIGN':
            continue
        rows.append(dict(task_id=r['task_id'], participant='X60_' + r['participant'], track='X60_roof', force05_N=None, force05_interval_N=None, cementfilm_um=None, cementfilm_interval_um=None, physical_status='UNKNOWN_BATCH_SUPPORT_CONTACT_AND_SPATIAL_FILM', conditional_force_ratio=r.get('optimization', {}).get('model_force_quantile_ratio'), conditional_force05_N=r.get('conditional_model_force05_N'), conditional_force_scope=r.get('conditional_model_force05_scope'), conditional_uncertainty='No rigorous design-box or empirical transfer enclosure', resolution='PER_TOOTH'))
    if (ROOT / 'raw/FUNCTIONAL_ROWS.json').exists():
        for r in read(ROOT / 'raw/FUNCTIONAL_ROWS.json') + (read(ROOT / 'raw/R6_EXTERNAL_SUPPORT.json')['rows'] if (ROOT / 'raw/R6_EXTERNAL_SUPPORT.json').exists() else []):
            if r['status'] != 'SCORED':
                continue
            rows.append(dict(task_id=r['key'], participant=r['participant'], track=r['track'], force05_N=None, force05_interval_N=None, cementfilm_um=None, cementfilm_interval_um=None, physical_status='UNKNOWN_UNMATCHED_SPECIMEN_AND_SPATIAL_FILM', digital_spacing=r.get('digital_spacing'), resolution='PER_TOOTH'))
    blob = json.dumps(clean(rows), separators=(',', ':')).encode()
    (ROOT / 'raw/PHYSICAL_ROWS.json.gz').write_bytes(gzip.compress(blob, mtime=0))
    out = dict(claim_type='capability', candidate_records=len(rows), calibrated_generated_force_count=0, calibrated_generated_film_count=0, published_force05=quant, calibration_rows=observations, direct_same_information_control=direct, corrected_fault_tests=tests, dropout=dict(literature_input_rows=148, primary_region_cells=94, other_or_incompatible_region_cells=54, eligible_force_calibration_targets=sum((r['kind'] == 'crown' for r in observations)), rejected_calibration_blocks=p['rejected_blocks'], generated_calibration_rejections=len(rows), reason='No matched generated-specimen external observations; source groups not independent crown predictions'), X1B_failed_transfer=read(ROOT.parent / 'LANE_X1B_CROWN_LOOP/raw/CALIBRATION_R1.json')['gates'], X13_physical_debt='Measured spatial seated-film field plus matched manufacture/seat process; regional dry-fit means are insufficient', X23_debt='Stress+surface flaw intensity+load/support/cement jointly needed; process-specific bar strength cannot certify an arbitrary crown', same_information_control='Corrected v4 direct reciprocal interpolation and log-force interpolation recomputed at source; no method superiority claim', seconds=time.perf_counter() - start)
    dump(ROOT / 'raw/R3.json', out)
    (ROOT / 'history/HANDOFF_R3.md').write_text('R3 records every generated design physical port and rechecks published calibration with corrected observation tests. Generated absolute fracture quantiles and seated-film remain unidentifiable. Next operation: use matched lab observations, not stronger metadata or extrapolated means; frozen blinded images and measurement templates prepared in R4.\n')
    return out
if __name__ == '__main__':
    run()
