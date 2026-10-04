"""Small exact planar-wall contact port. Units: um and dimensionless projections.

This port requires physical signed clearance bounds. RMS, color-map region maxima
and group medians are refused. A bridge shares a single axial insertion state.
"""
import itertools, math

def seat(clearance, projections):
    if len(clearance) != len(projections) or not clearance:
        raise ValueError('Empty or incompatible contact arrays')
    if any((not math.isfinite(x) for x in clearance + projections)) or min(projections) <= 0:
        raise ValueError('Need finite positive axis projections; different insertion branch UNKNOWN')
    lift = max(0.0, max((-b / a for (b, a) in zip(clearance, projections))))
    gap = [b + a * lift for (b, a) in zip(clearance, projections)]
    return (lift, gap)

def box_bounds(spacers, error_intervals, projections):
    lo = [s + e[0] for (s, e) in zip(spacers, error_intervals)]
    hi = [s + e[1] for (s, e) in zip(spacers, error_intervals)]
    if len(spacers) != len(error_intervals) or len(spacers) != len(projections) or (not spacers):
        raise ValueError('Port dimensions differ')
    if any((not math.isfinite(x) for x in lo + hi + projections)) or min(projections) <= 0 or any((l > h for (l, h) in zip(lo, hi))):
        raise ValueError('Invalid finite interval/projection')
    gl = []
    gu = []
    for (j, aj) in enumerate(projections):
        gl.append(max([0.0, lo[j]] + [lo[j] - aj * hi[i] / ai for (i, ai) in enumerate(projections) if i != j]))
        gu.append(max([0.0, hi[j]] + [hi[j] - aj * lo[i] / ai for (i, ai) in enumerate(projections) if i != j]))
    hl = max(0.0, max((-b / a for (b, a) in zip(hi, projections))))
    hu = max(0.0, max((-b / a for (b, a) in zip(lo, projections))))
    return {'gap_lower_um': gl, 'gap_upper_um': gu, 'lift_interval_um': [hl, hu]}

def vertex_control(spacers, intervals, projections):
    gaps = []
    lifts = []
    residual = 0.0
    for eps in itertools.product(*intervals):
        b = [s + e for (s, e) in zip(spacers, eps)]
        (h, g) = seat(b, projections)
        gaps.append(g)
        lifts.append(h)
        residual = max(residual, max(0.0, -min(g)))
    return {'gap_lower_um': [min((g[j] for g in gaps)) for j in range(len(spacers))], 'gap_upper_um': [max((g[j] for g in gaps)) for j in range(len(spacers))], 'lift_interval_um': [min(lifts), max(lifts)], 'penetration_um': residual}

def certify_port(port):
    if port.get('observation_kind') != 'whole_contact_signed_clearance_bound':
        return {'status': 'ABSTAIN_UNKNOWN', 'reason': 'Need full contact-region signed normal clearance bound; summary statistic cannot substitute'}
    if not port.get('datum_id') or port.get('datum_id') != port.get('consumer_datum_id'):
        return {'status': 'ABSTAIN_UNKNOWN', 'reason': 'Physical datum missing or mismatched'}
    if port.get('process_bound_known') is not True:
        return {'status': 'ABSTAIN_UNKNOWN', 'reason': 'Manufacturing bound UNKNOWN'}
    bounds = box_bounds(port['spacers_um'], port['error_intervals_um'], port['axis_projections'])
    if bounds['lift_interval_um'][1] > port['valid_lift_um']:
        return {'status': 'ABSTAIN_UNKNOWN', 'reason': 'Seating leaves validated geometry regime', 'bounds': bounds}
    j = port['margin_index']
    limit = port['gap_limit_um']
    status = 'ACCEPT_CONDITIONAL_BOX' if bounds['gap_upper_um'][j] < limit else 'REJECT_CONDITIONAL_BOX' if bounds['gap_lower_um'][j] >= limit else 'ABSTAIN_STRADDLES_BOUNDARY'
    return {'status': status, 'bounds': bounds, 'scope': 'DECLARED_SIGNED_BOUNDS_AND_PLANAR_CONTACT_ONLY'}
