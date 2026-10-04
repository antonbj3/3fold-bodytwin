"""Licensed specialization of the committed two-contact circular grasp.

Exact rational model only. The result does not certify an anatomical hand,
pressure, torque, finite patch deformation, or empirical friction bounds.
"""
from fractions import Fraction as Q
from itertools import product
import hashlib
import json


def rational(x):
    if type(x) not in (int, str, Q):
        raise ValueError('exact integer or rational string required')
    try:
        return Q(x)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError('finite rational required') from exc


def canonical(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'), ensure_ascii=False)


def digest(x):
    return hashlib.sha256(canonical(x).encode()).hexdigest()


def interval(lo, hi, unit, kind='EXACT_MODEL_BOUND'):
    lo, hi = rational(lo), rational(hi)
    if lo > hi:
        raise ValueError('unordered interval')
    return {'lower': str(lo), 'upper': str(hi), 'unit': unit,
            'uncertainty_kind': kind}


def normalize(request):
    if set(request) != {'J', 'inverse_mass', 'compliance', 'step_s',
                        'gaps_m', 'mu', 'free_velocity_m_s'}:
        raise ValueError('unknown or missing model fields')
    if request['J'] != [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                        [-1, 0, 0], [0, 1, 0], [0, 0, 1]]:
        raise ValueError('fixed opposed pure translation frames required')
    if list(map(rational, request['inverse_mass'])) != [Q(1)] * 3:
        raise ValueError('M=I specialization required')
    if list(map(rational, request['compliance'])) != [Q(1)] * 2:
        raise ValueError('normal compliance=1 specialization required')
    h = rational(request['step_s'])
    v = list(map(rational, request['free_velocity_m_s']))
    if h <= 0 or len(v) != 3 or v[0] != 0:
        raise ValueError('positive step and zero free normal velocity required')
    g = [tuple(map(rational, b)) for b in request['gaps_m']]
    mu = [tuple(map(rational, b)) for b in request['mu']]
    if len(g) != 2 or len(mu) != 2:
        raise ValueError('two contacts required')
    if any(len(b) != 2 or not b[0] <= b[1] < 0 for b in g):
        raise ValueError('strictly closed ordered gap boxes required')
    if any(len(b) != 2 or not 0 <= b[0] <= b[1] for b in mu):
        raise ValueError('nonnegative ordered friction boxes required')
    return h, v, g, mu


def certify(request):
    h, v, g, mu = normalize(request)
    # Normal block [[2,-1],[-1,2]] has principal minors 2 and 3.
    # Its inverse has nonnegative coefficients. Every pn is positive and
    # affine decreasing in both gaps. Capacity increases in each mu and pn.
    rows = []
    for g1, g2, m1, m2 in product(*g, *mu):
        pn = [-(2*g1+g2)/(3*h), -(g1+2*g2)/(3*h)]
        capacity = m1*pn[0]+m2*pn[1]
        need2 = v[1]**2+v[2]**2
        hold = capacity**2 >= need2
        tangents = None
        if hold and capacity:
            tangents = [[-v[j]*mi*pi/capacity for j in [1, 2]]
                        for mi, pi in zip([m1, m2], pn)]
            assert all(sum(t*t for t in ts) <= (mi*pi)**2
                       for ts, mi, pi in zip(tangents, [m1, m2], pn))
        assert 2*pn[0]-pn[1]+g1/h == 0
        assert -pn[0]+2*pn[1]+g2/h == 0
        rows.append({'gaps_m': [str(g1), str(g2)], 'mu': [str(m1), str(m2)],
                     'normal_impulses_N_s': list(map(str, pn)),
                     'capacity_N_s': str(capacity), 'vx_m_s': str(pn[0]-pn[1]),
                     'tangent_witness_N_s': None if tangents is None else
                     [list(map(str, ts)) for ts in tangents]})
    minimum = min(Q(r['capacity_N_s']) for r in rows)
    maximum = max(Q(r['capacity_N_s']) for r in rows)
    need2 = v[1]**2+v[2]**2
    all_hold = minimum**2 >= need2
    any_slip = minimum**2 < need2
    every_slip = maximum**2 < need2
    # For v_t != 0 maximum dissipation implies impulse=-C v_t/|v_t|,
    # hence |vf_t|=|v_t|+C. If |vf_t|<=C, v_t must be zero.
    bounds = {
        'normal_force_per_finger': [interval(min(Q(r['normal_impulses_N_s'][i])/h for r in rows),
                                           max(Q(r['normal_impulses_N_s'][i])/h for r in rows), 'N')
                                    for i in range(2)],
        'tangential_capacity': interval(minimum/h, maximum/h, 'N'),
        'required_tangent_force_squared': interval(need2/h**2, need2/h**2, 'N^2'),
        'margin_squared': interval(minimum**2-need2, maximum**2-need2, '(N*s)^2'),
        'vx': interval(min(Q(r['vx_m_s']) for r in rows), max(Q(r['vx_m_s']) for r in rows), 'm/s')}
    if all_hold:
        bounds.update(vy=interval(0, 0, 'm/s'), vz=interval(0, 0, 'm/s'))
    return {'schema': 'hum-circular-grasp-certificate/v1', 'request': request,
            'input_sha256': digest(request),
            'license': 'SPD_NORMAL_AFFINE_AND_MONOTONE_DISK_CAPACITY',
            'decision': 'HOLDS_ENTIRE_BOX' if all_hold else
                        'SLIPS_ENTIRE_BOX' if every_slip else 'REGIME_DEPENDENT',
            'model_status': 'EXACT_CONDITIONAL_MODEL_BOUND',
            'physical_status': 'UNKNOWN', 'any_slip_witness': any_slip,
            'bounds': bounds, 'corner_witnesses': rows}


def verify(certificate):
    try:
        return certificate == certify(certificate['request'])
    except (KeyError, TypeError, ValueError, ArithmeticError, AssertionError):
        return False


def source_request():
    return {'J': [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, 1, 0], [0, 0, 1]],
            'inverse_mass': ['1']*3, 'compliance': ['1']*2, 'step_s': '1/100',
            'gaps_m': [['-1/500', '-3/2500']]*2,
            'mu': [['1/2', '1/2']]*2,
            'free_velocity_m_s': ['0', '-981/10000', '-3/100']}
