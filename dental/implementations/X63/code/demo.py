from dental_release.paths import expand as _release_expand
import argparse, datetime, hashlib, json, os, time, resource, zipfile, struct
from pathlib import Path
import numpy as np
from lxml import etree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from assembly import beam_matrix, state, independent_state, exact_first_event
ROOT = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
OFJ = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Open-Full-Jaw'))
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X63-implant'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, d):
    Path(p).write_text(json.dumps(d, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def source_read():
    p = CORPUS / 'PMC10033796.xml'
    r = etree.parse(str(p))
    rows = []
    for (i, tr) in enumerate(r.findall('.//table-wrap[@id="Tab1"]//tr')):
        s = [' '.join(x.itertext()).strip() for x in tr]
        if len(s) == 6 and s[0].isdigit():
            rows.append({'from': int(s[0]), 'to': int(s[1]), 'difference_N': float(s[2]), 'lo_N': float(s[3]), 'hi_N': float(s[4]), 'p': float(s[5]), 'locator': str(p) + '#table-wrap[@id="Tab1"]//tr[' + str(i + 1) + ']'})
    text = ' '.join(r.getroot().itertext())
    for literal in ['329.9', '253.7', '255.7', '383.9', '200.1', '332.5', '25']:
        assert literal in text, literal
    return {'doi': '10.1186/s40729-023-00473-3', 'path': str(p), 'sha256': sha(p), 'rows': rows, 'first': {'mean_N': 329.9, 'sd_N': 33.3, 'observed_range_N': [255.7, 383.9]}, 'tenth': {'mean_N': 253.7, 'sd_N': 36.8, 'observed_range_N': [200.1, 332.5]}, 'n': 25, 'torque_Ncm': 25, 'resolution_level': 'POPULATION', 'regime': 'Dry Nobel Replace carbon-coated screw; intentionally separated 0.1mm force-measuring fixture. Tightening reuse events, not chewing cycles or elapsed time.', 'limitations': 'Table intervals are adjusted contrasts, not simultaneous support bounds or patient prediction intervals.'}

def anatomy():
    DATA.mkdir(parents=True, exist_ok=True)
    assert sum((p.stat().st_size for p in DATA.glob('*') if p.is_file())) < 3000000000
    p = OFJ / _release_expand('dataset/@DENTAL_CASE_ID@.zip')
    with zipfile.ZipFile(p) as z:
        member = _release_expand('@DENTAL_CASE_ID@/input/mandible/teeth_axes_mandible.json')
        raw = z.read(member)
        axes = json.loads(raw)
        points = np.array([axes[str(t)]['c'] for t in [20, 22, 27, 29]])
        bones = z.read(_release_expand('@DENTAL_CASE_ID@/input/mandible/bone/mandible.stl'))
    n = struct.unpack_from('<I', bones, 80)[0]
    assert len(bones) == 84 + 50 * n
    dtype = np.dtype([('normal', '<f4', (3,)), ('verts', '<f4', (3, 3)), ('attr', '<u2')])
    tris = np.frombuffer(bones, offset=84, count=n, dtype=dtype)['verts']
    view = tris[::max(1, n // 15000)].reshape(-1, 3).copy()
    np.savez_compressed(DATA / 'patient_view.npz', points=view)
    lengths = np.round(np.linalg.norm(np.diff(points, axis=0), axis=1), 2)
    return {'case': _release_expand('@DENTAL_CASE_ID@'), 'proxy_tooth_ids': [20, 22, 27, 29], 'points_mm': points.tolist(), 'segment_lengths_raw_mm': np.linalg.norm(np.diff(points, axis=0), axis=1).tolist(), 'segment_lengths_quantized_mm': lengths.tolist(), 'member': member, 'member_sha256': hashlib.sha256(raw).hexdigest(), 'zip_sha256': sha(p), 'bone_member_sha256': hashlib.sha256(bones).hexdigest(), 'length_quantization_max_error_mm': float(np.max(abs(lengths - np.linalg.norm(np.diff(points, axis=0), axis=1)))), 'resolution_level': 'PER_POINT', 'support_roles': 'Virtual support proxies, tooth centers on dentate jaw. Not implanted or approved surgical sites.'}

def freeze(round_id):
    src = source_read()
    geo = anatomy()
    write(ROOT / 'raw/SOURCE_PRELOAD.json', src)
    write(ROOT / 'raw/PATIENT_MODEL.json', geo)
    pred = {'frozen_at': now(), 'round': round_id, 'PREREG_sha256': sha(ROOT / f'PREREG_{round_id}.json'), 'code_sha256': {p.name: sha(p) for p in (ROOT / 'code').glob('*.py')}, 'source_sha256': src['sha256'], 'patient_member_sha256': geo['member_sha256'], 'static_preload_prediction_tenth_N': 329.9, 'model': 'Unfolded vertical rail; four tensile bolts and compression seats', 'fixed_parameters': {'E_MPa': 110000, 'rail_width_mm': 4, 'rail_height_mm': 4, 'kb_N_per_mm': 20000, 'kc_N_per_mm': 200000}, 'patterns_mm': {'zero': [0, 0, 0, 0], 'alternating': [0.02, -0.02, 0.02, -0.02], 'end': [0.02, 0, 0, 0], 'rigid_shift': [0.02, 0.02, 0.02, 0.02]}, 'load_direction': [0, 0, 0, 1], 'predictions_before_evaluation': ['All-gap summary cannot identify internal contact forces; require identical represented closed gaps.', 'Rigid shift has identical reaction/frontier to zero; local distortion changes it.', 'A constant first-use force fails source tenth-use mean within5%.', 'Patient absolute contact reserve remains unvalidated.']}
    write(ROOT / f'FROZEN_PREDICTIONS_{round_id}.json', pred)
    (ROOT / f'FROZEN_PREDICTIONS_{round_id}.sha256').write_text(sha(ROOT / f'FROZEN_PREDICTIONS_{round_id}.json') + '\n')
    if round_id == 'R1':
        write(ROOT / 'FROZEN_PREDICTIONS.json', pred)
    return (pred, src, geo)

def serialize(s):
    return {k: v.tolist() if isinstance(v, np.ndarray) else float(v) if isinstance(v, np.floating) else v for (k, v) in s.items()}

def beam_verification():
    from fractions import Fraction as Q
    from assembly import rational_solve
    K = beam_matrix([10], exact=True)
    A = [r[2:] for r in K[2:]]
    x = rational_solve(A, [Q(1), Q(0)])
    closed = Q(10) ** 3 / (3 * Q(110000) * Q(64, 3))
    return {'computed_tip_compliance_fraction_mm_per_N': str(x[0]), 'closed_form_fraction': str(closed), 'relative_error': float(abs(x[0] / closed - 1)), 'gate_pass': x[0] == closed}

def run_r1():
    t = time.perf_counter()
    (pred, src, geo) = freeze('R1')
    lengths = geo['segment_lengths_quantized_mm']
    K = beam_matrix(lengths)
    rows = []
    for (name, z) in pred['patterns_mm'].items():
        a = state(K, np.array(z), 329.9, np.zeros(4))
        b = independent_state(K, np.array(z), 329.9, np.zeros(4))
        cert = exact_first_event(lengths, z, 329.9, pred['load_direction'])
        rows.append({'pattern': name, 'source_preload_scenario_N': 329.9, 'assembly': serialize(a), 'independent_control': serialize(b), 'control_max_force_error_N': float(max(np.max(abs(a['bolt'] - b['bolt'])), np.max(abs(a['contact'] - b['contact'])))), 'frontier': cert, 'resolution_level': 'PHENOMENOLOGICAL'})
    (zero, alt) = rows[:2]
    ident = float(np.max(abs(np.array(zero['assembly']['gap']) - np.array(alt['assembly']['gap']))))
    delta = abs(zero['frontier']['first_event_N'] - alt['frontier']['first_event_N'])
    err = abs(329.9 / 253.7 - 1)
    out = {'round': 'R1', 'claim_type': ['capability', 'information_link'], 'source': src, 'patient_model': geo, 'rows': rows, 'constant_preload': {'predicted_N': 329.9, 'measured_tenth_mean_N': 253.7, 'relative_error': err, 'frozen_tolerance': 0.05, 'gate': 'FAIL' if err > 0.05 else 'PASS', 'resolution_level': 'POPULATION'}, 'sufficiency': {'summary': 'four closed assembly gaps', 'identity_error_mm': ident, 'downstream_first_event_difference_N': delta, 'gate': 'SUMMARY_REFUTED' if ident == 0 and delta > 1 else 'NOT_REFUTED', 'minimal_extension_fixed_query': 'Per-support signed contact reaction plus the load-response direction and branch guards. A preload scalar alone is insufficient.'}, 'beam_verification': beam_verification(), 'control_gate': all((r['control_max_force_error_N'] < 1e-05 for r in rows)), 'physical_patient_gate': 'UNKNOWN_NO_MATCHED_PRELOAD_COMPLIANCE_MISFIT_MEASUREMENT', 'cost': {'wall_s': time.perf_counter() - t, 'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'fit_s': 0, 'threads': 1, 'GPU': False}}
    write(ROOT / 'raw/RESULTS_R1.json', out)
    return out

def run_r2():
    t = time.perf_counter()
    (pred, src, geo) = freeze('R2')
    lengths = geo['segment_lengths_quantized_mm']
    rows = []
    for h in [3, 4, 5]:
        K = beam_matrix(lengths, h=h)
        for g in [0, 0.02, 0.05, 0.1, 0.15]:
            z = [g, -g, g, -g]
            for (use, p) in [('first_mean', 329.9), ('tenth_mean', 253.7), ('tenth_observed_min', 200.1), ('tenth_observed_max', 332.5)]:
                a = state(K, np.array(z), p, np.zeros(4))
                cert = exact_first_event(lengths, z, p, [0, 0, 0, 1], h=h)
                rows.append({'height_mm': h, 'misfit_pattern_mm': z, 'use_condition': use, 'preload_scenario_N': p, 'closed_assembly_gaps': a['gap'].tolist(), 'contact_reactions_N': a['contact'].tolist(), 'bolt_tension_N': a['bolt'].tolist(), 'frontier': cert, 'screen_150um': g <= 0.15, 'research_test_100N_status': 'NO_INITIAL_FULL_CONTACT' if cert['status'] == 'ASSEMBLY_BRANCH_INVALID' else 'MODEL_RESERVE_GE100N' if cert['first_event_N'] >= 100 else 'MODEL_RESERVE_LT100N', 'resolution_level': 'PHENOMENOLOGICAL'})
    from scipy.optimize import linprog
    A = []
    b = []
    for d in src['rows']:
        row = np.zeros(10)
        row[d['from'] - 1] = 1
        row[d['to'] - 1] = -1
        A.extend([row, -row])
        b.extend([d['hi_N'], -d['lo_N']])
    for i in range(9):
        row = np.zeros(10)
        row[i + 1] = 1
        row[i] = -1
        A.append(row)
        b.append(0)
    bounds = [(0, None)] * 10
    bounds[0] = (329.9, 329.9)
    bounds[9] = (253.7, 253.7)
    envelope = []
    for j in range(10):
        c = np.zeros(10)
        c[j] = 1
        lo = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
        hi = linprog(-c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
        assert lo.success and hi.success
        envelope.append({'tightening_event': j + 1, 'lo_mean_scenario_N': float(lo.fun), 'hi_mean_scenario_N': float(-hi.fun)})
    out = {'round': 'R2', 'claim_type': ['capability', 'information_link'], 'rows': rows, 'source_condition_envelope': envelope, 'envelope_scope': 'Feasible population-mean scenarios under quoted marginal contrast intervals plus monotone closure. NOT a simultaneous confidence region or a future screw bound.', 'source_preload_mean_drop_N': 76.2, 'source_preload_mean_drop_fraction': 1 - 253.7 / 329.9, 'sample_range_transfer': 'Observed minmax is a sensitivity scenario, not a physical enclosure.', 'physical_patient_gate': 'UNKNOWN_NO_MATCHED_MEASUREMENTS', 'cost': {'wall_s': time.perf_counter() - t, 'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'GPU': False}}
    write(ROOT / 'raw/RESULTS_R2.json', out)
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--round', choices=['R1', 'R2'], required=True)
    args = p.parse_args()
    out = run_r1() if args.round == 'R1' else run_r2()
    print(json.dumps({'round': out['round'], 'wall_s': out['cost']['wall_s'], 'physical_patient_gate': out['physical_patient_gate']}))
if __name__ == '__main__':
    main()
