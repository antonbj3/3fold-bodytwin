"""Same-part independent signed metrology port with supplied pointwise bounds.

The bound must cover the signed observable AND its coordinate/correspondence
error; an instrument accuracy specification is not that certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from metrology import write_json, REGIONS

def site_hash(data):
    h = hashlib.sha256()
    for key in ['xyz_mm', 'normals', 'design_face', 'regions']:
        a = np.ascontiguousarray(data[key])
        h.update(str(a.dtype).encode())
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()

def compare(field_path, reference_path, specimen_id, out):
    data = np.load(field_path)
    ref = json.loads(Path(reference_path).read_text())
    if ref.get('kind') not in ('independent_measurement', 'our_own_fixture'):
        raise ValueError('Explicit reference provenance kind required; unknown source cannot calibrate')
    if 'specimen_id' not in data.files or str(data['specimen_id']) != specimen_id or specimen_id == 'UNKNOWN':
        raise ValueError('Scanner field has no matching physical specimen identity')
    if ref.get('site_sha256') != site_hash(data):
        raise ValueError('Reference sites/frame/normals do not match')
    if ref.get('specimen_id') != specimen_id:
        raise ValueError('Independent reference is for a different physical part')
    if ref.get('units') != 'um' or ref.get('observable') != 'signed_normal_deviation_at_frozen_design_sites':
        raise ValueError('Reference units or observable mismatch')
    if not ref.get('locator') or not ref.get('bound_locator'):
        raise ValueError('Reference and pointwise bound need locators')
    r = np.asarray(ref['reference_deviation_um'], float)
    b = np.asarray(ref['pointwise_error_bound_um'], float)
    d = data['normal_deviation_mm'].mean(0) * 1000
    if r.shape != d.shape or b.shape != d.shape or (not np.isfinite(r).all()) or (not np.isfinite(b).all()) or np.any(b < 0):
        raise ValueError('Finite nonnegative bound required independently for every site')
    lower = np.nextafter(r - b, -np.inf)
    upper = np.nextafter(r + b, np.inf)
    diff_lower = np.nextafter(d - upper, -np.inf)
    diff_upper = np.nextafter(d - lower, np.inf)
    regions = {}
    for group in REGIONS:
        take = data['regions'] == group
        w = data['area_weights_mm2'][take]
        gamma = len(w) * np.finfo(float).eps / (1 - len(w) * np.finfo(float).eps)

        def weighted_interval(lo, hi):
            a = lo[take] * w
            c = hi[take] * w
            numerator = [a.sum() - gamma * abs(a).sum(), c.sum() + gamma * abs(c).sum()]
            denominator = [w.sum() * (1 - gamma), w.sum() * (1 + gamma)]
            values = [x / y for x in numerator for y in denominator]
            return [float(np.nextafter(min(values), -np.inf)), float(np.nextafter(max(values), np.inf))]
        regions[group] = dict(resolution='PER_SURFACE_REGION', manufacture_signed_mean_interval_um=weighted_interval(lower, upper), mean_scan_minus_reference_interval_um=weighted_interval(diff_lower, diff_upper))
    result = dict(claim_type='capability', specimen_id=specimen_id, site_sha256=site_hash(data), reference_locator=ref['locator'], bound_locator=ref['bound_locator'], resolution='PER_POINT', timescale='HANDOVER', regions=regions, manufacture_lower_um=lower.tolist(), manufacture_upper_um=upper.tolist(), mean_scan_minus_true_reference_lower_um=diff_lower.tolist(), mean_scan_minus_true_reference_upper_um=diff_upper.tolist(), reference_kind=ref.get('kind', 'UNKNOWN'), empirical_calibration_status='UNKNOWN_SIMULATED_REFERENCE' if ref.get('kind') == 'our_own_fixture' else 'CONDITIONAL_ON_SUPPLIED_REFERENCE_CERTIFICATE', rigorous_enclosure='Conditional on supplied independent signed pointwise bound; outward float rounding', limitations=['Finite sites, no full-surface enclosure', 'Finite repeat mean, not asymptotic scanner bias', 'Reference bound validity is supplied, not established by this comparison'])
    write_json(out, result)
    return result

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fields', required=True)
    p.add_argument('--reference', required=True)
    p.add_argument('--specimen-id', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    try:
        compare(a.fields, a.reference, a.specimen_id, a.output)
    except (ValueError, KeyError) as e:
        p.exit(2, 'REFERENCE_REJECT: ' + str(e) + '\n')
if __name__ == '__main__':
    main()
