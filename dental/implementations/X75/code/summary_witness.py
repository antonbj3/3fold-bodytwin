"""Exact-summary witnesses and enclosed, hypothetical contrast transport.

No CBCT gray-to-temperature or gray-to-density conversion is made.
The reference is the published homogeneous whole-line heat kernel.
"""
import hashlib, json, math, pathlib, time
from fractions import Fraction
import mpmath as mp
import numpy as np
from scipy.special import erf
P = pathlib.Path(__file__).resolve().parents[1]

def put(n, x):
    (P / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def bounds_for_bins(n, sub=256):
    mp.iv.dps = 35
    I = mp.iv
    norm = I.sqrt(4 * I.pi)
    step = I.mpf(4) / (n * sub)
    out = []
    for j in range(n):
        lower = I.mpf(0)
        upper = I.mpf(0)
        for k in range(sub):
            a = I.mpf(4 * (j * sub + k)) / (n * sub)
            b = I.mpf(4 * (j * sub + k + 1)) / (n * sub)
            lower += I.exp(-b * b / 4) / norm * step
            upper += I.exp(-a * a / 4) / norm * step
        out.append(I.mpf([lower.a, upper.b]))
    return out

def enclosed(vals, weights):
    I = mp.iv
    s = I.mpf(0)
    for (v, w) in zip(vals, weights):
        s += I.mpf(v.numerator) / v.denominator * w
    return [float(np.nextafter(float(s.a), -np.inf)), float(np.nextafter(float(s.b), np.inf))]

def run():
    start = time.perf_counter()
    pr = json.loads((P / 'PREREG_R2.json').read_text())
    assert hashlib.sha256((P / 'PREREG_R2.json').read_bytes()).hexdigest() == (P / 'PREREG_R2.sha256').read_text().strip()
    rows = [json.loads(l) for l in (P / 'raw/PER_SITE_R2.jsonl').open()]
    geometric = []
    for r in rows:
        if r['status'] != 'MEASURED':
            continue
        q = r['rays'][1]
        hi = q['buccal_extent_mm']
        side = q['lingual_extent_mm']
        if hi is None:
            continue
        widthA = hi + side
        widthB = side + hi
        A = min(hi - 1, side + 1) - r['radius_mm']
        B = min(side - 1, hi + 1) - r['radius_mm']
        geometric.append({'case': r['case'], 'site': r['site'], 'width_A_mm': widthA, 'width_B_mm': widthB, 'identity_error_mm': abs(widthA - widthB), 'summary_bitwise_identical': np.float64(widthA).tobytes() == np.float64(widthB).tobytes(), 'fixed_shift_mm': 1, 'regional_min_gap_A_mm': A, 'regional_min_gap_B_mm': B, 'downstream_difference_mm': abs(A - B), 'minimal_extension_for_this_pair': 'signed side extents (or width plus signed center offset) in process frame', 'referent_kind': 'our_own_fixture', 'interpretation': 'reflection of measured source geometry, not a second measured patient or a clinical allowable shift'})
    heat = None
    for r in sorted(rows, key=lambda r: (r['case'], r['site'])):
        if r['status'] != 'MEASURED' or not r['kind'].startswith('X8'):
            continue
        with np.load(r['raw_file']['path']) as raw:
            g = raw['core_gray'][:32].copy()
        if len(g) < 32 or max(g) <= min(g):
            continue
        lo = Fraction.from_float(float(min(g)))
        span = Fraction.from_float(float(max(g))) - lo
        exact = sorted(((Fraction.from_float(float(v)) - lo) / span for v in g))
        reverse = exact[::-1]
        n = len(exact)
        edges = np.linspace(0, 4, n + 1)
        w = (erf(edges[1:] / 2) - erf(edges[:-1] / 2)) / 2
        A = np.array([float(v) for v in exact])
        B = np.array([float(v) for v in reverse])
        candA = float(A @ w)
        candB = float(B @ w)
        mp.mp.dps = 60
        wm = [(mp.erf(mp.mpf(4 * (i + 1)) / n / 2) - mp.erf(mp.mpf(4 * i) / n / 2)) / 2 for i in range(n)]
        refA = float(sum((mp.mpf(v.numerator) / v.denominator * z for (v, z) in zip(exact, wm))))
        refB = float(sum((mp.mpf(v.numerator) / v.denominator * z for (v, z) in zip(reverse, wm))))
        wb = bounds_for_bins(n, 256)
        bA = enclosed(exact, wb)
        bB = enclosed(reverse, wb)
        meanA = float(sum(exact, Fraction(0)) / n)
        meanB = float(sum(reverse, Fraction(0)) / n)
        histA = np.sort(A).tobytes()
        histB = np.sort(B).tobytes()
        gates = {'identity_exact': histA == histB and np.float64(meanA).tobytes() == np.float64(meanB).tobytes(), 'downstream_separation': bB[0] - bA[1] >= 0.01, 'closed_form_control': max(abs(candA - refA), abs(candB - refB)) <= 1e-10, 'model_enclosure': bA[0] <= refA <= bA[1] and bB[0] <= refB <= bB[1] and (max(bA[1] - bA[0], bB[1] - bB[0]) <= 0.002)}
        np.savez_compressed(P / 'raw/EXACT_SUMMARY_WITNESS.npz', gray_values=g, initial_A=A, initial_B=B, bin_edges=edges, weights=w)
        heat = {'source_case': r['case'], 'source_site': r['site'], 'source_raw_file': r['raw_file'], 'state_definition': 'constructed ascending/reversed unit contrast in32 bins over[0,4]; unbounded physical gray-to-heat closure', 'identity_error_histogram': 0 if histA == histB else 1, 'identity_error_mean': abs(meanA - meanB), 'mean_A': meanA, 'mean_B': meanB, 'functional_A_unit_contrast': candA, 'functional_B_unit_contrast': candB, 'difference_unit_contrast': candB - candA, 'rigorous_model_intervals_A_B': [bA, bB], 'difference_lower_bound': bB[0] - bA[1], 'external_closed_form_reference': [refA, refB], 'max_control_error': max(abs(candA - refA), abs(candB - refB)), 'minimal_extension_for_fixed_query': 'one kernel-weighted spatial projection; histogram/mean alone fails', 'extension_for_new_location_time_or_transport': 'retain spatial state plus independently calibrated transport/observation operator', 'numerical_guarantee': 'outward interval Gaussian Riemann sums; monotone kernel on y>=0, nonnegative exact rational initial values', 'physical_remainder_enclosure': 'MISSING; no absolute physical sensitivity or temperature prediction', 'resolution': 'PHENOMENOLOGICAL', 'external_referent': pr['external_referents'][1], 'constructed_states_referent': {'kind': 'our_own_fixture', 'locator': str(P / 'raw/EXACT_SUMMARY_WITNESS.npz'), 'compared_quantity': 'exact histogram twins', 'refutes_us': not all(gates.values())}, 'gates': gates}
        break
    grayrows = json.loads((P / 'raw/CONTROLS_R2.json').read_text())
    validgray = next((q for q in grayrows if q['kind'] == 'gray'))
    actualgrayerror = validgray['error']
    geom = next((q for q in grayrows if q['kind'] == 'ray_box' and q['candidate'] is not None))
    mutatedray = list(geom['candidate'])
    mutatedray[0] += 0.1
    injectedrayerror = max((abs(x - y) for (x, y) in zip(mutatedray, geom['control'])))
    fault = {'gray_guard_rejects_plus1': actualgrayerror == 0 and abs(actualgrayerror + 1) > 0, 'geometry_guard_rejects_plus0p1mm': injectedrayerror > 1e-08, 'heat_control_rejects_plus0p1': heat is not None and abs(heat['functional_A_unit_contrast'] + 0.1 - heat['external_closed_form_reference'][0]) > 1e-10, 'interval_rejects_shifted_prediction': heat is not None and (not heat['rigorous_model_intervals_A_B'][0][0] <= heat['functional_A_unit_contrast'] + 0.1 <= heat['rigorous_model_intervals_A_B'][0][1])}
    put('raw/GEOMETRY_SUMMARY_WITNESSES.json', geometric)
    put('SUMMARY_SUFFICIENCY.json', {'geometry': {'pairs': len(geometric), 'identity_max_error_mm': max([g['identity_error_mm'] for g in geometric], default=0), 'max_downstream_difference_mm': max([g['downstream_difference_mm'] for g in geometric], default=0), 'nonzero_downstream_pairs': sum((g['downstream_difference_mm'] > 0 for g in geometric)), 'scope': 'EXPLORATORY_DESCRIPTIVE_NOT_PREREGISTERED: source-derived mathematical reflection; not independent patient truth; side-offset is minimum separating input for tested gap'}, 'heat': heat, 'fault_injection': fault, 'cost_wall_s': time.perf_counter() - start, 'physical_measurements': 0})
    print(json.dumps({'heat': heat, 'fault_injection': fault}, indent=2))
if __name__ == '__main__':
    run()
