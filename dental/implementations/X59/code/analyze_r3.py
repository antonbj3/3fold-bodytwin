"""Independent source stress test and a detector contrast at unchanged film state."""
from decimal import Decimal
from fractions import Fraction as F
import hashlib, itertools, json, math, time
import numpy as np
from extract_sources import ROOT, save, extract_cunali, extract_pasha, extract_smartphone, sha
from analyze_r2 import contrast

def add(a, b):
    return (a[0] + b[0], a[1] + b[1])

def neg(a):
    return (-a[1], -a[0])

def sub(a, b):
    return add(a, neg(b))

def mul(a, b):
    vals = [x * y for x in a for y in b]
    return (min(vals), max(vals))

def div(a, b):
    if b[0] <= 0 <= b[1]:
        raise ValueError('Denominator enclosure includes zero')
    return mul(a, (1 / b[1], 1 / b[0]))

def sq(a):
    if a[0] <= 0 <= a[1]:
        return (F(0), max((x * x for x in a)))
    return (min((x * x for x in a)), max((x * x for x in a)))

def isum(xs):
    ans = (F(0), F(0))
    for x in xs:
        ans = add(ans, x)
    return ans

def box(s):
    return (F(s) - F(1, 200), F(s) + F(1, 200))

def floats(a):
    return [math.nextafter(float(a[0]), -math.inf), math.nextafter(float(a[1]), math.inf)]

def ols_enclosure(rows):
    out = []
    for (sys, published) in [('Amann', ('0.4736', '36.206')), ('Dentsply_Sirona', ('0.4695', '44.553'))]:
        rr = [r for r in rows if r['system'] == sys]
        xs = [box(r['replica_mean_um']) for r in rr]
        ys = [box(r['ct_mean_um']) for r in rr]
        n = (F(4), F(4))
        sx = isum(xs)
        sy = isum(ys)
        sxy = isum([mul(x, y) for (x, y) in zip(xs, ys)])
        sxx = isum([sq(x) for x in xs])
        num = sub(sxy, div(mul(sx, sy), n))
        den = sub(sxx, div(sq(sx), n))
        slope = div(num, den)
        intercept = sub(div(sy, n), mul(slope, div(sx, n)))
        sp = F(published[0])
        ip = F(published[1])
        sfig = (sp - F(1, 20000), sp + F(1, 20000))
        ifig = (ip - F(1, 2000), ip + F(1, 2000))
        intersect = lambda a, b: max(a[0], b[0]) <= min(a[1], b[1])
        corners = []
        for bits in itertools.product([0, 1], repeat=8):
            x = [xs[j][bits[j]] for j in range(4)]
            y = [ys[j][bits[j + 4]] for j in range(4)]
            scalar_a = (sum((a * b for (a, b) in zip(x, y))) - sum(x) * sum(y) / 4) / (sum((a * a for a in x)) - sum(x) ** 2 / 4)
            scalar_b = sum(y) / 4 - scalar_a * sum(x) / 4
            assert slope[0] <= scalar_a <= slope[1] and intercept[0] <= scalar_b <= intercept[1]
            corners.append((scalar_a, scalar_b))
        out.append(dict(system=sys, resolution='PER_SURFACE_REGION', slope_rounding_enclosure=floats(slope), intercept_rounding_enclosure_um=floats(intercept), exact_rational_slope_enclosure=[str(v) for v in slope], exact_rational_intercept_enclosure=[str(v) for v in intercept], published=dict(slope=published[0], intercept_um=published[1]), published_rounding_intervals=dict(slope=floats(sfig), intercept_um=floats(ifig)), coefficient_rounding_gate='NOT_EXCLUDED_BY_RIGOROUS_ROUNDING_ENCLOSURE' if intersect(slope, sfig) and intersect(intercept, ifig) else 'EXCLUDED', independent_corner_checks=len(corners), proof='Exact rational inclusion arithmetic over continuous eight-dimensional source-rounding box; corner enumeration is an additional executed check.', limitations='Separate enclosure membership need not imply joint attainability. No specimen or target prediction error is enclosed. R2 failed point-coefficient gate remains failed.'))
    return out

def freeze_predictions(donors, targets):
    pred = []
    for r in targets:
        reg = 'axial' if r['region'] == 'Axial' else 'occlusal' if r['region'] == 'Occlusal cavosurface area' else None
        if reg is None:
            continue
        for donor in [d for d in donors if d['region'] == reg]:
            p = Decimal(r['replica_mean_um']) - Decimal(donor['bias_replica_minus_ct_um'])
            pred.append(dict(source=r['source'], system=r['system'], region=r['region'], donor_system=donor['system'], predicted_ct_mean_um=str(p), resolution='PER_SURFACE_REGION', alignment='OUT_OF_SOURCE_SCOPE_RESTORATION_GEOMETRY_AND_STATE' if reg == 'axial' else 'OUT_OF_SOURCE_SCOPE_QUANTITY_GEOMETRY_AND_STATE', source_locator=r['source_locator']))
    freeze = dict(round='R3', created_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), claim_type='information_link', prospective=False, disclosure='Retrospective: source tables were already inspected during selection. Saved before scorer, not before original publication or blind physical measurement.', predictions=pred, source_hashes={'cunali': sha(ROOT / 'raw/cunali2017.pdf'), 'pasha': sha(targets[0]['source_path'])}, prereg_sha256=sha(ROOT / 'PREREG_R3.json'))
    path = ROOT / 'FROZEN_PREDICTIONS_R3.json'
    if path.exists():
        old = json.loads(path.read_text())
        if old['predictions'] != pred or old['source_hashes'] != freeze['source_hashes'] or old['prereg_sha256'] != freeze['prereg_sha256']:
            raise ValueError('R3 freeze drift')
    else:
        save('FROZEN_PREDICTIONS_R3.json', freeze)
        (ROOT / 'FROZEN_PREDICTIONS_R3.sha256').write_text(sha(path) + '\n')
    return pred

def run():
    start = time.perf_counter()
    cunali = extract_cunali()
    pasha = extract_pasha()
    smart = extract_smartphone()
    donors = [contrast(r) for r in cunali]
    pred = freeze_predictions(donors, pasha)
    scored = []
    for p in pred:
        r = next((r for r in pasha if r['region'] == p['region'] and r['system'] == p['system']))
        error = abs(Decimal(p['predicted_ct_mean_um']) - Decimal(r['ct_mean_um']))
        scored.append(dict(**p, observed_ct_mean_um=r['ct_mean_um'], absolute_error_um=str(error), raw_no_correction_error_um=str(abs(Decimal(r['replica_mean_um']) - Decimal(r['ct_mean_um']))), numeric_20um_gate='PASS' if error <= 20 else 'FAIL', physical_transport_gate='FAIL_DIFFERENT_OBSERVATION_DOMAIN', external_referent=dict(kind='independent_measurement', locator=r['source_locator'], compared_quantity='reported same-onlay cohort method mean', refutes_us=True)))
    device = []
    for r in smart:
        diff = Decimal(r['phone_mean_um']) - Decimal(r['microscope_mean_um'])
        device.append(dict(**r, phone_minus_microscope_mean_um=str(diff), rounding_interval_um=[str(diff - Decimal('.010')), str(diff + Decimal('.010'))], identified_effect='Detector/localization/pixel-scale group contrast conditional on same sectioned film', actual_loa=None, covariance_status='UNKNOWN_UNSPECIFIED_ICC_IS_NOT_PEARSON_CORRELATION', external_referent=dict(kind='independent_measurement', locator=r['source_locator'], compared_quantity='device-specific cohort means on same film', refutes_us=True)))
    pasha_contrasts = [contrast(r) for r in pasha]
    save('INDEPENDENT_METHOD_CONTRASTS.json', pasha_contrasts)
    save('SAME_FILM_DEVICE_CONTRASTS.json', device)
    result = dict(round='R3', claim_type='information_link', outcome='OBSERVATION_STATE_SEPARATION_DELIVERED_UNIVERSAL_TRANSPORT_NOT_SUPPORTED', external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.7759/cureus.40020 TAB1/TAB2; https://doi.org/10.1038/s41598-024-55711-4 Tab2', compared_quantity='same-sample cross-method region means and same-film device contrasts', refutes_us=True), transfer_stress_test=scored, accepted_in_domain_validation_strata=0, specimen_pairs_available=0, table_pair_rows=len(pasha), pasha_missing_table_rows=2, pasha_table_row_exclusion_fraction=2 / 12, device_region_pairs=len(device), device_mean_contrast_min_max_um=[min((float(r['phone_minus_microscope_mean_um']) for r in device)), max((float(r['phone_minus_microscope_mean_um']) for r in device))], device_contrast_resolution='PER_SURFACE_REGION', rounding_ols=ols_enclosure(cunali), sufficiency_test_ref='RESULTS_R2.json /sufficiency_test; exact identity error 0.0 for marginal summaries, classification error changes by 4/4', equally_informed_control='Direct conventional paired-cohort arithmetic yields the same contrast; exact rational OLS enclosure is checked on all 256 rounding-box vertices.', next_operation='Registered same-coping dry/PVS-film CT plus optical measurements; record actual pair rows and force/pose. No calibration is inferred from high ecological correlation or unspecified ICC.', seconds=time.perf_counter() - start)
    save('RESULTS_R3.json', result)
    print(json.dumps(dict(round='R3', stress_comparisons=len(scored), numeric_gates=[r['numeric_20um_gate'] for r in scored], device_region_pairs=len(device), device_mean_contrast_range_um=result['device_mean_contrast_min_max_um'], rounding_gates=[r['coefficient_rounding_gate'] for r in result['rounding_ols']])))
    return result
if __name__ == '__main__':
    run()
