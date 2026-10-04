"""Domain-guarded source fixture port. An endpoint is not a validated force law."""
from fractions import Fraction
import math
FIXTURE = 'TPBT_L10_E4H_PEEK_unbonded'

def query(direction, displacement_mm, fixture, unit='mm'):
    if fixture != FIXTURE:
        raise ValueError('unmeasured fixture/material/span transfer')
    if unit != 'mm':
        raise ValueError('displacement unit must be mm')
    if direction not in ['AP', 'vertical']:
        raise ValueError('unmeasured direction')
    d = float(displacement_mm)
    if not math.isfinite(d):
        raise ValueError('nonfinite displacement')
    delta = Fraction(str(d))
    source_delta = Fraction('.06' if direction == 'AP' else '.02')
    force = Fraction('1.27' if direction == 'AP' else '.95')
    if delta < 0 or delta > source_delta:
        raise ValueError('outside measured source excursion')
    lo = (force - Fraction('.005')) / source_delta * delta
    hi = (force + Fraction('.005')) / source_delta * delta
    return dict(force_linear_scenario_N=float(force / source_delta * delta), printed_rounding_linear_scenario_interval_N=[float(lo), float(hi)], exact_rounding_rationals=[str(lo), str(hi)], endpoint_replay=delta == source_delta, status='PRINTED_ENDPOINT_REPLAY' if delta == source_delta else 'CONDITIONAL_LINEAR_SCENARIO', physical_force_guarantee=False, constitutive_remainder='UNKNOWN', resolution='PHENOMENOLOGICAL', timescale='SIMULTANEOUS')
