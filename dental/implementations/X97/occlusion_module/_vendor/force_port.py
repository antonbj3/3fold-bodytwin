from fractions import Fraction as Q

def project(observation):
    try:
        meta = observation['context']
        if meta['time_relation'] != 'SIMULTANEOUS' or meta['unit'] != 'force_share_same_total':
            return dict(status='UNKNOWN', reason='time_or_unit_mismatch')
        if not meta.get('same_pose_id') or not meta.get('common_denominator_id') or (not meta.get('incidence_sha256')):
            return dict(status='UNKNOWN', reason='unbound_pose_denominator_or_incidence')
        if meta['evidence_kind'] == 'MEASURED':
            if not meta.get('calibration_locator') or meta.get('gain_uncertainty_status') != 'BOUNDED':
                return dict(status='UNKNOWN', reason='sensor_to_force_operator_unknown')
        elif meta['evidence_kind'] != 'SIMULATED':
            return dict(status='UNKNOWN', reason='unknown_evidence_kind')
        (ql, qh) = map(Q, observation['upper_region_interval_share'])
        (dl, dh) = map(Q, observation['signed_net_crossing_interval_share'])
        if not 0 <= ql <= qh <= 1 or dl > dh or dl < -1 or (dh > 1):
            return dict(status='REJECTED', reason='invalid_measurement_interval')
        lo = max(Q(0), ql + dl)
        hi = min(Q(1), qh + dh)
        if lo > hi:
            return dict(status='REJECTED', reason='no_conserved_nonnegative_lower_force')
        return dict(status='CONDITIONAL_BOUND', lower_region_interval_share=list(map(str, [lo, hi])), relation='lower=upper+signed_net_crossing', rigor='Exact rational interval sum intersect[0,1]; conservative when input errors correlated; no constitutive or floating sensitivity approximation', evidence_kind=meta['evidence_kind'], level='PER_ARCH', input_resolution='PER_TOOTH')
    except (KeyError, ValueError, TypeError, ZeroDivisionError) as e:
        return dict(status='REJECTED', reason='malformed_observation:' + type(e).__name__)
