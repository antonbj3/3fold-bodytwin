"""Exact measurement-box transport. Fraction calculations bound all vertices.

The declared +/-1 um error is a prospective instrumentation requirement,
not an uncertainty estimated from published sensor resolution.
"""
from fractions import Fraction as Q
from pathlib import Path
import json, math, itertools, hashlib, datetime
from motion_port import local_threshold
LANE = Path(__file__).resolve().parents[1]

def inv2(a):
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    if det == 0:
        raise ValueError('rank deficient')
    return [[a[1][1] / det, -a[0][1] / det], [-a[1][0] / det, a[0][0] / det]]

def mm(a, b):
    return [[sum((x * y for (x, y) in zip(row, col))) for col in zip(*b)] for row in a]

def directed(x, direction):
    v = float(x)
    if direction == 'lo' and Q(v) > x:
        v = math.nextafter(v, -math.inf)
    if direction == 'hi' and Q(v) < x:
        v = math.nextafter(v, math.inf)
    return v

def response_bounds(z_markers, W, w, z_query, implant_readings, bone_readings, error_mm):
    if bone_readings is None:
        return {'status': 'UNKNOWN_MISSING_LOCAL_BONE_REFERENCE'}
    H = [[Q(1), Q(z)] for z in z_markers]
    Hq = [[Q(1), Q(z)] for z in z_query]
    try:
        left = mm(Hq, inv2(H))
        right = mm(inv2(W), [[x] for x in w])
    except ValueError:
        return {'status': 'UNKNOWN_RANK_DEFICIENT_CALIBRATION'}
    vals = []
    coeffs = []
    for lrow in left:
        coef = [lrow[i] * right[j][0] for i in range(2) for j in range(2)]
        center = sum((coef[i * 2 + j] * (implant_readings[i][j] - bone_readings[i][j]) for i in range(2) for j in range(2)))
        radius = sum((abs(c) * 2 * error_mm for c in coef))
        vals.append([center - radius, center + radius])
        coeffs.append(coef)
    return {'status': 'BOUNDED_UNDER_DECLARED_LINEAR_BRANCH', 'rational_intervals': vals, 'coefficients': coeffs}

def main():
    z = [Q(-4), Q(4)]
    W = [[Q(1), Q(0)], [Q(0), Q(4)]]
    w = [Q(1), Q(4)]
    zq = [Q(-4), Q(0), Q(4)]
    err = Q(1, 1000)
    freeze = {'round': 'R3', 'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'reading_error_mm': str(err), 'marker_z_mm': [str(x) for x in z], 'predictions': 'interval maxima separate <50 and >150 in the two R2 fixtures', 'kind': 'virtual measurement only; published source data predates this freeze'}
    f = LANE / 'FROZEN_PREDICTIONS_R3.json'
    if not f.exists():
        f.write_text(json.dumps(freeze, indent=2) + '\n')
    (LANE / 'FROZEN_PREDICTIONS_R3.sha256').write_text(hashlib.sha256(f.read_bytes()).hexdigest() + '  ' + f.name + '\n')
    outcases = {}
    for (name, K2) in [('spread_support', Q(1024)), ('central_support', Q(64))]:
        C = [[Q(1, 64), Q(0)], [Q(0), 1 / K2]]
        H = [[Q(1), a] for a in z]
        bone = [[Q(1, 8), Q(1, 16)], [Q(1, 8), Q(1, 16)]]
        local = mm(mm(H, C), W)
        implant = [[local[i][j] + bone[i][j] for j in range(2)] for i in range(2)]
        r = response_bounds(z, W, w, zq, implant, bone, err)
        intervals = r['rational_intervals']
        failures = 0
        max_rounding_outside = 0
        for sign in itertools.product([-1, 1], repeat=8):
            yi = [[implant[i][j] + sign[i * 2 + j] * err for j in range(2)] for i in range(2)]
            yb = [[bone[i][j] + sign[4 + i * 2 + j] * err for j in range(2)] for i in range(2)]
            exact = response_bounds(z, W, w, zq, yi, yb, Q(0))['rational_intervals']
            for (actual, (lo, hi)) in zip(exact, intervals):
                if not lo <= actual[0] <= hi:
                    failures += 1
                (flo, fhi) = (directed(lo, 'lo'), directed(hi, 'hi'))
                if not Q(flo) <= actual[0] <= Q(fhi):
                    max_rounding_outside += 1
        absbounds = []
        for (lo, hi) in intervals:
            abslo = Q(0) if lo <= 0 <= hi else min(abs(lo), abs(hi))
            absbounds.append([abslo, max(abs(lo), abs(hi))])
        peak = [max((x[0] for x in absbounds)) * 1000, max((x[1] for x in absbounds)) * 1000]
        peakf = [directed(peak[0], 'lo'), directed(peak[1], 'hi')]
        held = mm(mm([[Q(1), a] for a in zq], C), [[x] for x in w])
        held_pass = all((lo - 2 * err <= actual[0] <= hi + 2 * err for (actual, (lo, hi)) in zip(held, intervals)))
        bad = implant.copy()
        bad = [row.copy() for row in bad]
        bad[1][1] += Q(10, 1000)
        badr = response_bounds(z, W, w, zq, bad, bone, err)['rational_intervals']
        bad_pass = all((lo - 2 * err <= actual[0] <= hi + 2 * err for (actual, (lo, hi)) in zip(held, badr)))
        outcases[name] = {'rigorous_sampled_peak_interval_um': peakf, 'rational_sampled_peak_interval_um': [str(v) for v in peak], 'corner_count': 256, 'rational_corner_failures': failures, 'binary64_corner_failures': max_rounding_outside, 'declared_error_um_each_reading': 1, 'error_status': 'PHENOMENOLOGICAL_PROSPECTIVE_SPECIFICATION', 'held_third_wrench_pass': held_pass, 'injected_10um_error_pass': bad_pass, 'injected_error_rejected': not bad_pass, 'mechanical_sampled_region_decision': local_threshold(peakf, True), 'biological_decision': 'UNKNOWN: sampled linear mechanical model has no healing closure', 'spatial_whole_interface_peak': 'UNKNOWN: no bound on local nonrigid residual beyond sampled regions', 'resolution': 'PER_SURFACE_REGION', 'time_scale': 'SIMULTANEOUS', 'missing_bone_reference_query': response_bounds(z, W, w, zq, implant, None, err), 'absolute_interval_width_um': peakf[1] - peakf[0]}
        assert failures == 0 and max_rounding_outside == 0 and held_pass and (not bad_pass)
    out = {'round': 'R3', 'claim_type': 'capability', 'verdict': 'EXACT_LINEAR_MEASUREMENT_BOX_QUERY_DELIVERED; EMPIRICAL_INTERFACE_CALIBRATION_UNKNOWN', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'R2 declared spring equilibrium, exact Fraction corners in code/round3.py', 'compared_quantity': 'affine image of declared implant and bone reading-error boxes', 'refutes_us': False}, 'independent_measurement_constraint': {'kind': 'independent_measurement', 'locator': 'doi:10.1186/s12903-024-03854-1 Table 1', 'compared_quantity': '43/43 um published axial mean with 57/112 um total means; need transverse/rotation information', 'refutes_us': True}, 'cases': outcases, 'rigorous_enclosure_scope': 'all declared independent reading-error combinations; calibrated positions and wrenches exact; same linear branch; no dental constitutive uncertainty enclosed', 'minimum_real_measurement': '2 separated implant markers and colocated bone reference markers under 2 independent controlled in-plane wrenches; third wrench validation; full 3D needs 6 DOF and sufficient wrenches', 'remaining_obstacle': 'spatially unresolved interface slip and evolving nonlinear contact; no source supplies this calibration per implant type and human bone class', 'controls': 'exact rational full-corner enumeration; 10 um injected reading rejected; no empirical accuracy claim from fixture', 'full_cost': {'calibration_readings': 8, 'independent_wrenches': 2, 'validation_wrenches': 1, 'candidate_worlds': 2, 'corner_evaluations': 512, 'query_fit_count': 0, 'large_arrays': 0}}
    (LANE / 'rounds/R3').mkdir(parents=True, exist_ok=True)
    (LANE / 'rounds/R3/results.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(outcases, indent=2))
if __name__ == '__main__':
    main()
