"""Nonlinear fatigue stress port and rigorous brackets for its declared closure.

For compression c>=0, t=g*F>c, q=(1-Rmachine)/2, B=min(q*t/(t-c),1):
q*t=C*B**alpha. C>0,0<alpha<1. The uncapped equation is
q**(1-alpha)*t**(1-alpha)*(t-c)**alpha=C, strictly increasing in t.
Thus t increases with C,c and decreases with alpha. Force=t/g.
Natural arbitrary-precision interval evaluation certifies each bisection sign.
Only this scalar empirical closure is enclosed, not titanium validity or FE error.
"""
import math, re
from mpmath import iv
iv.dps = 50

def bracket_stress(C, c, alpha, R):
    if not (C > 0 and 0 <= c and (0 < alpha < 1) and (0 <= R < 1)):
        raise ValueError('invalid tensile fatigue port')
    q = (1 - R) / 2
    if not iv.mpf(c) < iv.mpf(C) / ((1 - iv.mpf(R)) / 2):
        raise ValueError('compression offset reaches resistance ceiling; tensile branch not certified')
    low = c
    high = math.nextafter(float((iv.mpf(C) / ((1 - iv.mpf(R)) / 2)).b), math.inf)
    steps = 0
    while steps < 100:
        mid = (low + high) / 2
        if mid == low or mid == high:
            break
        ti = iv.mpf(mid)
        ci = iv.mpf(c)
        qi = (1 - iv.mpf(R)) / 2
        B = qi * ti / (ti - ci)
        if B.a >= 1:
            B = iv.mpf(1)
        elif B.b > 1:
            break
        residual = qi * ti - iv.mpf(C) * B ** iv.mpf(alpha)
        if residual.b < 0:
            low = mid
        elif residual.a > 0:
            high = mid
        else:
            break
        steps += 1
    return (low, high, steps)

def force_bounds(g, C, c, alpha, R):
    if not (0 < g[0] <= g[1] and 0 < C[0] <= C[1] and (0 <= c[0] <= c[1]) and (0 < alpha[0] <= alpha[1] < 1)):
        raise ValueError('invalid measured interval')
    if not iv.mpf(c[1]) < iv.mpf(C[0]) / ((1 - iv.mpf(R)) / 2):
        raise ValueError('entire box not in tensile branch')
    (low, _, n1) = bracket_stress(C[0], c[0], alpha[1], R)
    (_, high, n2) = bracket_stress(C[1], c[1], alpha[0], R)
    return {'force_interval_N': [math.nextafter(low / g[1], -math.inf), math.nextafter(high / g[0], math.inf)], 'interval_arithmetic_precision_decimal_digits': 50, 'bisection_sign_steps': [n1, n2], 'enclosure_kind': 'nonlinear monotonic box enclosure with interval-certified signs', 'scope': 'declared scalar closure only; same critical region, elastic proportional loading, compressive fixed offset,0<alpha<1'}

def force_point(g, C, c, alpha, R):
    q = (1 - R) / 2
    low = c
    high = C / q
    for _ in range(100):
        t = (low + high) / 2
        if t <= c:
            low = t
            continue
        B = min(q * t / (t - c), 1.0)
        if q * t > C * B ** alpha:
            high = t
        else:
            low = t
    return (low + high) / 2 / g

def measured_query(port, R):
    required = ['system_id', 'critical_region', 'g_MPa_per_N', 'C_MPa', 'compression_offset_MPa', 'alpha', 'measurement_locators', 'provenance_kind', 'horizon_cycles', 'frozen_before_measurement_sha256', 'units']
    if any((k not in port for k in required)):
        return {'status': 'UNKNOWN', 'reason': 'missing same-configuration physical port', 'required_fields': required}
    if port['provenance_kind'] not in ['independent_measurement', 'published_dataset']:
        return {'status': 'UNKNOWN', 'reason': 'a fixture or simulation is not independent physical calibration'}
    if port['units'] != {'g': 'MPa/N', 'C': 'MPa', 'compression_offset': 'MPa', 'alpha': '1', 'horizon': 'cycles'}:
        return {'status': 'UNKNOWN', 'reason': 'stress/force/unit contract mismatch'}
    for k in ['g_MPa_per_N', 'C_MPa', 'compression_offset_MPa', 'alpha']:
        if not isinstance(port[k], list) or len(port[k]) != 2 or any((not isinstance(x, (int, float)) or not math.isfinite(x) for x in port[k])):
            return {'status': 'UNKNOWN', 'reason': 'missing/invalid physical interval: ' + k}
    if not re.fullmatch('[0-9a-fA-F]{64}', str(port['frozen_before_measurement_sha256'])):
        return {'status': 'UNKNOWN', 'reason': 'missing valid frozen-plan hash'}
    if len(port['measurement_locators']) < 3 or any((not s or 'UNKNOWN' in s for s in port['measurement_locators'])):
        return {'status': 'UNKNOWN', 'reason': 'requires stress slope, absolute static offset and reference fatigue force locators'}
    try:
        bounds = force_bounds(port['g_MPa_per_N'], port['C_MPa'], port['compression_offset_MPa'], port['alpha'], R)
    except ValueError as e:
        return {'status': 'UNKNOWN', 'reason': str(e)}
    return {'status': 'MODEL_CONDITIONAL_PREDICTION', 'machine_R': R, 'horizon_cycles': port['horizon_cycles'], 'system_id': port['system_id'], 'critical_region': port['critical_region'], 'resolution': 'PER_TOOTH', 'risk_certificate': False, **bounds}
