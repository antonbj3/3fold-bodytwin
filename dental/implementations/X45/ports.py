"""Dimensioned observation ports. No patient/material promotion by interpolation."""
import math
from fractions import Fraction as F
import mpmath as mp

class PortError(ValueError):
    pass

def fraction(value):
    return F(str(value))

def rational_bounds(lo, hi):
    """Convert exact rational bounds outward to IEEE binary64."""
    (a, b) = (float(lo), float(hi))
    if F.from_float(a) > lo:
        a = math.nextafter(a, -math.inf)
    if F.from_float(b) < hi:
        b = math.nextafter(b, math.inf)
    return [a, b]

def interpolate_enclosed(times, means, query, half_rounding='0.05'):
    t = [fraction(v) for v in times]
    y = [fraction(v) for v in means]
    q = fraction(query)
    if len(t) != len(y) or len(t) < 2 or any((a >= b for (a, b) in zip(t, t[1:]))):
        raise PortError('non-increasing time state')
    if q < t[0] or q > t[-1]:
        raise PortError('unmeasured time extrapolation')
    j = next((i for i in range(len(t) - 1) if t[i] <= q <= t[i + 1]))
    w = (q - t[j]) / (t[j + 1] - t[j])
    value = (1 - w) * y[j] + w * y[j + 1]
    e = fraction(half_rounding)
    return {'mean': float(value), 'rounding_enclosure': rational_bounds(value - e, value + e), 'weights': [str(1 - w), str(w)], 'resolution': 'POPULATION', 'timescale': 'HANDOVER', 'physical_interpolation_remainder': 'UNKNOWN', 'individual_prediction': 'UNKNOWN'}

def promote_isq(quantity):
    if quantity not in ['group_mean_ISQ', 'tabular_ISQ']:
        raise PortError('ISQ does not identify ' + quantity)

def weibull_quantile(m, scale, probability):
    if not (m > 0 and scale > 0 and (0 < probability < 1)):
        raise PortError('invalid Weibull parameter/probability')
    return scale * (-math.log1p(-probability)) ** (1 / m)

def weibull_rounding_enclosure(m, scale, probability):
    """mpmath directed interval arithmetic, decimal inputs, enclosing printed rounding only."""
    mp.iv.dps = 40
    mf = fraction(m)
    sf = fraction(scale)
    e = F(5, 1000)
    if mf - e <= 0 or sf - e <= 0:
        raise PortError('rounding interval includes invalid material')

    def exact_str(a):
        return mp.nstr(mp.mpf(a.numerator) / a.denominator, 50)

    def dec(a):
        from decimal import Decimal, localcontext
        with localcontext() as ctx:
            ctx.prec = 80
            return str(Decimal(a.numerator) / Decimal(a.denominator))
    mi = mp.iv.mpf([dec(mf - e), dec(mf + e)])
    si = mp.iv.mpf([dec(sf - e), dec(sf + e)])
    pi = mp.iv.mpf(str(probability))
    v = si * (-mp.iv.ln(1 - pi)) ** (1 / mi)
    return [math.nextafter(float(v.a), -math.inf), math.nextafter(float(v.b), math.inf)]
COUPON_PROTOCOL = {'temperature_C': 122, 'pressure_bar': 2, 'hours': 8, 'cyclic_load_N': 50, 'cycles': 240000, 'frequency_Hz': 1.1, 'specimen_diameter_mm': 15, 'specimen_thickness_mm': 1}

def material_query(record, probability, context='source_coupon', requested_protocol=None):
    if context != 'source_coupon':
        raise PortError('unvalidated stressed-volume/geometry transfer: ' + context)
    if requested_protocol is not None and requested_protocol != COUPON_PROTOCOL:
        raise PortError('exposure mismatch')
    q = weibull_quantile(record['m'], record['sigma0_MPa'], probability)
    return {'quantile_MPa': q, 'printed_rounding_enclosure_MPa': weibull_rounding_enclosure(record['m'], record['sigma0_MPa'], probability), 'sampling_CI': 'UNKNOWN', 'resolution': 'POPULATION', 'timescale': 'SIMULTANEOUS', 'clinical_years': 'UNKNOWN', 'geometry_transfer': 'UNKNOWN'}

def fit_rate(cycles, masses, intercept=False):
    """Exact rational normal equations and active-boundary enumeration."""
    x = [fraction(v) for v in cycles]
    y = [fraction(v) for v in masses]
    n = F(len(x))
    xx = sum((v * v for v in x))
    xy = sum((a * b for (a, b) in zip(x, y)))
    sx = sum(x)
    sy = sum(y)
    k = max(F(0), xy / xx)
    candidates = [(F(0), k)]
    if intercept:
        candidates.append((max(F(0), sy / n), F(0)))
        slope = (n * xy - sx * sy) / (n * xx - sx * sx)
        a = (sy - slope * sx) / n
        if a >= 0 and slope >= 0:
            candidates.append((a, slope))
    (a, k) = min(candidates, key=lambda ak: sum(((u - ak[0] - ak[1] * v) ** 2 for (v, u) in zip(x, y))))
    return (a, k)

def rate_query(coefficients, cycles):
    (a, k) = coefficients
    pred = a + k * fraction(cycles)
    return {'mean_mg': float(pred), 'arithmetic_enclosure_mg': rational_bounds(pred, pred), 'physical_extrapolation_remainder': 'UNKNOWN', 'resolution': 'POPULATION'}

def wear_envelope(cycles, masses, query):
    x = [fraction(v) for v in cycles]
    y = [fraction(v) for v in masses]
    q = fraction(query)
    if len(x) != len(y) or any((a >= b for (a, b) in zip(x, x[1:]))):
        raise PortError('unordered wear state')
    if any((a > b for (a, b) in zip(y, y[1:]))):
        raise PortError('nonmonotone mass-loss anchors')
    if not x[0] <= q <= x[-1]:
        raise PortError('unmeasured cycle extrapolation')
    j = next((i for i in range(len(x) - 1) if x[i] <= q <= x[i + 1]))
    e = F(5, 1000)
    if q == x[j]:
        (lo, hi) = (y[j] - e, y[j] + e)
    elif q == x[j + 1]:
        (lo, hi) = (y[j + 1] - e, y[j + 1] + e)
    else:
        (lo, hi) = (y[j] - e, y[j + 1] + e)
    return {'mass_interval_mg': rational_bounds(lo, hi), 'width_mg': float(hi - lo), 'assumption': 'cumulative mean mass loss is monotone; UNKNOWN hydration/deposition bias', 'SEM': 'not included in hard envelope', 'resolution': 'POPULATION', 'timescale': 'HANDOVER'}

def mass_to_volume(mass_mg, density_g_cm3=None):
    if density_g_cm3 is None or density_g_cm3 <= 0:
        raise PortError('same-material density missing')
    return mass_mg / density_g_cm3
