"""Decision-relative resolution budgets. No empirical parameter is fitted here.

Grid pitch is a deterministic numerical bound, never a standard deviation.
The caller must keep physical calibration separate from a conditional error law.
"""
from dataclasses import dataclass, asdict
from math import isfinite, sqrt, erf

@dataclass(frozen=True)
class ScalarPort:
    margin: float
    unit: str
    bias_bound: float = 0.0
    standard_deviations: tuple = ()
    coefficient: float = sqrt(3)
    exponent: float = 1.0
    k: float = 3.0
    covariance: str = 'UNKNOWN_WORST_COVARIANCE'
    local_gaussian_law_validated: bool = False
    physical_model_validated: bool = False

    def evaluate(self, pitch):
        vals = (self.margin, self.bias_bound, self.coefficient, self.exponent, self.k, pitch, *self.standard_deviations)
        if not all((isfinite(x) for x in vals)):
            raise ValueError('finite quantities required; missing values are UNKNOWN')
        if any((x < 0 for x in (self.bias_bound, self.coefficient, pitch, *self.standard_deviations))):
            raise ValueError('sigma/bounds/pitch must be nonnegative; signed bias needs an explicit map')
        if self.exponent <= 0 or self.k <= 0 or (not self.unit):
            raise ValueError('positive convergence exponent, k and declared unit required')
        s = sum(self.standard_deviations)
        floor = self.bias_bound + self.k * s
        slack = abs(self.margin) - floor
        klass = 1 if slack > 0 else 2
        limit = (slack / self.coefficient) ** (1 / self.exponent) if slack > 0 and self.coefficient else None
        sufficient = limit * 0.5 if limit is not None else pitch if slack > 0 else None
        envelope = floor + self.coefficient * pitch ** self.exponent
        certified = abs(self.margin) > envelope
        physical = certified and self.local_gaussian_law_validated and self.physical_model_validated
        return {'class': klass, 'margin': self.margin, 'unit': self.unit, 'bias_bound': self.bias_bound, 'sigma_envelope': s, 'floor': floor, 'numerical_error_bound': self.coefficient * pitch ** self.exponent, 'envelope': envelope, 'pitch': pitch, 'pitch_limit_strict': limit, 'suggested_sufficient_pitch': sufficient, 'conditional_status': ('CERTIFIED_POSITIVE' if self.margin > 0 else 'CERTIFIED_NEGATIVE') if certified else 'ABSTAIN', 'conditional_gaussian_point_coverage': erf(self.k / sqrt(2)), 'confidence_scope': 'conditional single scalar; NOT simultaneous surface coverage', 'physical_status': 'CERTIFIED' if physical else 'UNKNOWN', 'floor_is_conditional': True, 'port': asdict(self)}

def interval_decision(lower, upper, threshold, *, adverse_error=0.0):
    """All admissible values, with a deterministic error. No probability attached."""
    if not all((isfinite(x) for x in (lower, upper, threshold, adverse_error))):
        raise ValueError('finite interval required')
    if lower > upper or adverse_error < 0:
        raise ValueError('ordered interval and nonnegative error required')
    (lo, hi) = (lower - adverse_error, upper + adverse_error)
    if hi < threshold:
        return 'BELOW'
    if lo > threshold:
        return 'ABOVE'
    return 'ABSTAIN'

def refinement_sequence(values, pitches):
    """Diagnostic only: finite samples cannot prove divergence or a floor."""
    if len(values) != len(pitches) or len(values) < 3:
        return {'class': None, 'status': 'UNKNOWN_NO_REFINEMENT_SEQUENCE'}
    if not all((isfinite(x) for x in values + pitches)):
        return {'class': None, 'status': 'UNKNOWN_NONFINITE'}
    if not all((pitches[i] > pitches[i + 1] > 0 for i in range(len(pitches) - 1))):
        raise ValueError('strictly decreasing positive pitches required')
    return {'class': None, 'status': 'DIAGNOSTIC_ONLY_REQUIRES_ANALYTIC_OR_ERROR_ENVELOPE', 'successive_changes': [values[i + 1] - values[i] for i in range(len(values) - 1)]}
