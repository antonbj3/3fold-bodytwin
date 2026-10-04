"""Explicit measured-error/tissue contract. Empty or other-patient data never substitute."""
from pathlib import Path
import hashlib, math

def rejects(fn):
    try:
        fn()
        return False
    except (ValueError, TypeError, KeyError):
        return True

def registration_bound(record):
    if record.get('kind') != 'independent_directed_uniform_anchor_error':
        raise ValueError('Surface residual quantile is not a uniform directed registration bound')
    if record.get('independently_measured') is not True or not record.get('measurement_locator'):
        raise ValueError('No independent calibration locator')
    eps = record.get('epsilon_mm')
    if eps is None or not math.isfinite(eps) or eps < 0:
        raise ValueError('Unknown/invalid epsilon')
    return eps

def admit_tissue(record, tissue, patient_id):
    if record is None:
        raise ValueError('Same-patient tissue geometry missing')
    if record.get('patient_id') != patient_id:
        raise ValueError('Foreign patient tissue cannot be borrowed')
    if record.get('tissue') != tissue:
        raise ValueError('Wrong anatomical quantity')
    if record.get('frame') != 'CBCT_native':
        raise ValueError('CBCT-native transform missing')
    p = record.get('path')
    if not p or not Path(p).is_file():
        raise ValueError('Geometry path missing')
    if not record.get('locator') or record.get('units') != 'mm':
        raise ValueError('Independent locator or unit missing')
    if hashlib.sha256(Path(p).read_bytes()).hexdigest() != record.get('sha256'):
        raise ValueError('Geometry hash mismatch')
    eta = record.get('independent_Hausdorff_error_mm')
    if eta is None or not math.isfinite(eta) or eta < 0:
        raise ValueError('Uniform tissue geometry error missing')
    return record

def conditional_distance_interval(d, L, epsilon, eta_tissue, eta_prep=0.0):
    if any((x is None for x in [d, L, epsilon, eta_tissue, eta_prep])):
        return None
    if any((not math.isfinite(x) or x < 0 for x in [d, L, epsilon, eta_tissue, eta_prep])):
        raise ValueError('Invalid conditional distance inputs')
    e = math.nextafter(math.nextafter(math.nextafter(L * epsilon, math.inf) + eta_tissue, math.inf) + eta_prep, math.inf)
    return [max(0.0, math.nextafter(d - e, -math.inf)), math.nextafter(d + e, math.inf)]
