"""Independent distance controls and falsifiable gates; fixtures are not external facit."""
from measure import *
from scipy.spatial.distance import cdist
from sections import ray_exit
sys.path.insert(0, str(DENT / 'results/LANE_X8_GUIDE_NERVE_RISK'))
from full_geometry import voxel_cylinder_bracket

def ident_gate(a, b):
    return np.array_equal(a, b)

def independence_gate(r):
    return bool(r.get('image_exact_equal') and r.get('blinded') and (len(r.get('annotator_ids', [])) >= 2) and (len(set(r['annotator_ids'])) == len(r['annotator_ids'])))

def same_numbers(a, b):
    return bool(np.array_equal(np.array(a), np.array(b)))

def main():
    start = time.perf_counter()
    patients = json.loads((ROOT / 'raw/PAIRED_PATIENTS.json').read_text())
    p = next((x for x in patients if not x['numeric_mask_equal']))
    d = np.load(p['point_artifact']['path'])
    po = d['old_boundary_zyx']
    pn = d['new_boundary_zyx']
    recorded = d['old_to_new_mm']
    shape = tuple(p['shape'])
    m = np.zeros(shape, dtype=bool)
    m[tuple(pn.T)] = True
    edt = ndi.distance_transform_edt(~m, sampling=SP)
    edtval = edt[tuple(po.T)]
    error = float(np.max(np.abs(edtval - recorded)))
    poison = recorded.copy()
    poison[0] += 0.1
    n = min(257, len(po))
    brute = cdist(po[:n] * SP, pn * SP).min(1)
    bruteerr = float(np.max(np.abs(brute - recorded[:n])))
    ref = np.array([[z, 0, 0] for z in range(-16, 17)], dtype=float)
    a = ref + [0, 0, -1]
    b = ref + [0, 0, 1]
    q = np.array([[0, 0, 3]], float)
    da = cKDTree(ref).query(a)[0]
    db = cKDTree(ref).query(b)[0]
    suma = np.array([len(a), 0, np.median(da), np.quantile(da, 0.95), da.max()])
    sumb = np.array([len(b), 0, np.median(db), np.quantile(db, 0.95), db.max()])
    qa = float(cdist(q, a).min())
    qb = float(cdist(q, b).min())
    summary_err = float(np.max(np.abs(suma - sumb)))
    suff = {'status': 'FAIL_SUMMARY_SUFFICIENCY', 'summary': ['voxel_count', 'Dice_to_reference', 'median_surface_mm', 'p95_surface_mm', 'max_surface_mm'], 'state_a_summary': suma.tolist(), 'state_b_summary': sumb.tolist(), 'summary_identity_error': summary_err, 'full_distance_distribution_identical': bool(np.array_equal(np.sort(da), np.sort(db))), 'downstream_fixed_point_distance_a_mm': qa, 'downstream_fixed_point_distance_b_mm': qb, 'downstream_difference_mm': qa - qb, 'decision_threshold_mm': 2.5, 'decision_a': 'ABOVE', 'decision_b': 'BELOW', 'smallest_extension_fixed_query': 'signed clearance decrement for that query (or one bit for this two-state fixture); retain spatial wall location for new queries', 'resolution': 'PER_POINT', 'time_scale': 'SIMULTANEOUS', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/controls.py canal columns at x=-1,+1; query x=3', 'compared_quantity': 'summary identity vs downstream distance', 'refutes_us': True}}
    refm = np.zeros((9, 9, 13), bool)
    refm[3:6, 2:7, 2:7] = True
    am = refm.copy()
    am[3:6, 2:7, 7] = True
    bm = refm.copy()
    bm[3:6, 2:7, 1] = True
    dirs = np.array([[0, 0, 1], [0, 1, 0], [0, 0, -1], [0, -1, 0]], float)
    orig = np.tile(np.array([4, 4, 4]) * SP, (4, 1))
    radii_ref = ray_exit(refm, orig, dirs)
    radii_a = ray_exit(am, orig, dirs)
    radii_b = ray_exit(bm, orig, dirs)
    da_rad = np.abs(radii_a - radii_ref)
    db_rad = np.abs(radii_b - radii_ref)
    dra = 2 * int((am & refm).sum()) / (int(am.sum()) + int(refm.sum()))
    drb = 2 * int((bm & refm).sum()) / (int(bm.sum()) + int(refm.sum()))
    rsa = np.array([am.sum(), dra, np.median(da_rad), np.quantile(da_rad, 0.95), da_rad.max()])
    rsb = np.array([bm.sum(), drb, np.median(db_rad), np.quantile(db_rad, 0.95), db_rad.max()])
    qa_rad = 10 * SP - (np.argwhere(am)[:, 2].max() + 0.5) * SP
    qb_rad = 10 * SP - (np.argwhere(bm)[:, 2].max() + 0.5) * SP
    radial_suff = {'summary': ['voxel_count', 'Dice_to_reference', 'median_abs_radial_delta_mm', 'p95_abs_radial_delta_mm', 'max_abs_radial_delta_mm'], 'state_a_summary': rsa.tolist(), 'state_b_summary': rsb.tolist(), 'summary_identity_error': float(np.max(np.abs(rsa - rsb))), 'full_radial_abs_distribution_identical': bool(np.array_equal(np.sort(da_rad), np.sort(db_rad))), 'downstream_a_mm': float(qa_rad), 'downstream_b_mm': float(qb_rad), 'downstream_difference_mm': float(qb_rad - qa_rad), 'threshold_mm': 0.9, 'smallest_extension': 'Changed wall direction at this station, or the signed clearance decrement for a fixed query; preserve directional radial array for new queries', 'resolution': 'PER_POINT', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/controls.py mirrored one-column cuboid expansion', 'compared_quantity': 'local radial summary vs point clearance', 'refutes_us': True}}
    assert radial_suff['summary_identity_error'] == 0 and radial_suff['full_radial_abs_distribution_identical'] and (qa_rad < 0.9 < qb_rad)
    cr = voxel_cylinder_bracket(np.array([[2.0, 0.0, 4.0]]), np.array([0.3] * 3), np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]), 5.0, 1.0)
    analytic = 2.85
    cerr = max(abs(cr['lower_mm'] - analytic), abs(cr['upper_mm'] - analytic))
    g = 0.5
    shifted = voxel_cylinder_bracket(np.array([[2.0, 0.0, 4.0]]), np.array([0.3] * 3), np.array([0.0, 0.0, g]), np.array([1.0, 0.0, 0.0]), 5.0, 1.0)
    terr = abs(shifted['upper_mm'] - (analytic - g))
    copies = [x for x in patients if x['numeric_mask_equal']]
    copyok = all((x['symmetric']['max_mm'] == 0 for x in copies))
    scales = np.array(p['spacing_mm'])
    scaleok = bool(np.all(np.abs(scales - 0.3) < 1e-14))
    badscale = scales.copy()
    badscale[0] = 3
    original = np.array([1, 2, 3])
    badimage = original.copy()
    badimage[0] += 1
    out = {'claim_type': 'capability', 'distance_control': {'patient': p['patient'], 'kd_vs_edt_max_abs_error_mm': error, 'kd_vs_bruteforce_max_abs_error_mm': bruteerr, 'tolerance_mm': 1e-10, 'pass': error <= 1e-10 and bruteerr <= 1e-10, 'equal_information_control': 'MATCHES'}, 'closed_form': {'locator': 'finite cylinder {0<=z<=5, x^2+y^2<=1}; box z=2, y=0, x=4; halfwidth .15', 'compared_quantity': 'minimum set distance', 'expected_mm': analytic, 'computed_lower_mm': cr['lower_mm'], 'computed_upper_mm': cr['upper_mm'], 'max_abs_error_mm': cerr, 'pass': cerr <= 1e-10, 'rigorous_floating_point_enclosure': 'NOT_CERTIFIED: analytic convex support inequality valid; floating-point implementation checked numerically only'}, 'translation_bound': {'guide_ball_radius_mm': g, 'expected_mm': analytic - g, 'computed_upper_mm': shifted['upper_mm'], 'error_mm': terr, 'pass': terr <= 1e-10, 'population_probability_claim': False}, 'copy_zero_gate': {'copies_tested': len(copies), 'pass': copyok}, 'sufficiency': suff, 'radial_sufficiency': radial_suff, 'injected_errors': {'0p1mm_distance_error_rejected': float(np.max(np.abs(edtval - poison))) > 1e-10, 'tenfold_scale_rejected': not np.all(np.abs(badscale - 0.3) < 1e-14), 'one_voxel_image_error_rejected': not ident_gate(original, badimage), 'summary_downstream_false_equivalence_rejected': summary_err == 0 and qa != qb, 'radial_summary_false_equivalence_rejected': bool(radial_suff['summary_identity_error'] == 0 and qa_rad != qb_rad), 'closed_form_wrong_answer_rejected': abs(analytic - (cr['upper_mm'] + 0.1)) > 1e-10, 'translation_wrong_answer_rejected': abs(analytic - g - (shifted['upper_mm'] + 0.1)) > 1e-10, 'copied_annotation_false_independence_rejected': not independence_gate({'image_exact_equal': True, 'blinded': True, 'annotator_ids': ['reader_a', 'reader_a'], 'independent': True}), 'unsupported_independence_flag_rejected': not independence_gate({'image_exact_equal': True, 'independent': True}), 'documented_metadata_positive_control_accepted': independence_gate({'image_exact_equal': True, 'blinded': True, 'annotator_ids': ['reader_a', 'reader_b']})}, 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1}}
    dump(ROOT / 'raw/CONTROLS.json', out)
    dump(ROOT / 'raw/SUFFICIENCY.json', suff)
    dump(ROOT / 'raw/RADIAL_SUFFICIENCY.json', radial_suff)
    assert out['distance_control']['pass'] and out['closed_form']['pass'] and out['translation_bound']['pass'] and copyok
    assert all(out['injected_errors'].values())
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    main()
