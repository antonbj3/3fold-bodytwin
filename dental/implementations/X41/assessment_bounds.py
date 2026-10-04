"""Conditional reference intervals. Submitted declarations are not anatomical proof."""
import math

def reference_interval(gap_mm, threshold_mm, manifest=None, numerical_mm=1e-09):
    if not all((math.isfinite(x) and x >= 0 for x in (gap_mm, threshold_mm, numerical_mm))):
        raise ValueError('Finite nonnegative gap, threshold and numerical allowance required')
    required = ('reference_error_bound_mm', 'registration_error_bound_mm', 'body_error_bound_mm')
    if manifest is None:
        return {'decision': 'UNKNOWN_MISSING_BOUND_INFORMATION', 'missing': list(required), 'physical_status': 'UNKNOWN', 'resolution_level': 'PER_SURFACE_REGION'}
    missing = [k for k in required if manifest.get(k) is None]
    if missing:
        return {'decision': 'UNKNOWN_MISSING_BOUND_INFORMATION', 'missing': missing, 'physical_status': 'UNKNOWN', 'resolution_level': 'PER_SURFACE_REGION'}
    values = [manifest[k] for k in required]
    if any((isinstance(v, bool) or not isinstance(v, (int, float)) or (not math.isfinite(v)) or (v < 0) for v in values)):
        raise ValueError('Hard bounds must be finite, nonnegative mm numbers; SD is not a bound')
    if manifest.get('coverage') != 'COMPLETE_SHAPES':
        return {'decision': 'UNKNOWN_INCOMPLETE_BOUND_COVERAGE', 'physical_status': 'UNKNOWN', 'reason': 'A reference subset cannot lower-bound the complete hazard minimum', 'resolution_level': 'PER_SURFACE_REGION'}
    total = sum(values) + numerical_mm
    lo = max(0.0, gap_mm - total)
    hi = gap_mm + total
    decision = 'ABOVE' if lo >= threshold_mm else 'BELOW' if hi < threshold_mm else 'ABSTAIN'
    return {'decision': decision, 'clearance_interval_mm': [lo, hi], 'total_declared_bound_plus_roundoff_mm': total, 'numerical_allowance_mm': numerical_mm, 'physical_status': 'UNKNOWN', 'bound_status': 'CALLER_DECLARED_CONDITIONAL_GEOMETRY_ONLY', 'resolution_level': 'PER_SURFACE_REGION', 'time_scale': 'SIMULTANEOUS', 'source_locators': manifest.get('source_locators', []), 'requirements': 'Independent complete-shape anatomical and registration bounds need review; no injury or clinical guarantee'}
