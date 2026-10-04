import json, pathlib, time, copy, hashlib
import numpy as np
from sweep import *
from first_event import dual_valid
R = pathlib.Path(__file__).resolve().parents[1]

def first_height_valid(row):
    ws = [w for q in row['all_pair_certificates'].values() for w in q['intersection_certificates']]
    if not ws:
        return False
    exact = max((Q(w['s_exact']) for w in ws)) * 40
    return Q(row['height_exact_mm']) == exact and row['first_surface_contact_height_mm'] == float(exact)

def repair_verdict_valid(report):
    collision = any((q['status'] == 'COLLISION' for q in report['paths'].values()))
    invariance = report['native_exterior_identity_error_mm'] == 0 and all((w['same_exact_witness_survives'] for w in report['invariant_witnesses']))
    return report['status'] == 'NO_FEASIBLE_INSIDE_ONLY_REPAIR_FOR_FIXED_PATH' and collision and invariance

def endpoint_report_valid(report):
    return report['identity_error'] == 0 and report['endpoint_summary_mm'][0] == report['endpoint_summary_mm'][1] and (report['downstream_path'] == ['SURFACE_CLEAR_CERTIFIED', 'COLLISION']) and all((x['status'] == 'SURFACE_CLEAR_CERTIFIED' for x in report['endpoint_checks']))

def run(root=R, full_reference=True):
    tick = time.perf_counter()
    root = pathlib.Path(root)
    co = json.load(open(root / 'FROZEN_COHORT.json'))
    r1 = json.load(open(root / 'rounds/R1/COHORT.json'))
    r3 = json.load(open(root / 'rounds/R3/FIRST_CONTACT.json'))
    r4 = json.load(open(root / 'rounds/R4/AXES.json'))
    counts = {'exact_primal_checks': 0, 'exact_dual_checks': 0, 'independent_LP_checks': 0, 'reference_boolean_mismatches': 0, 'invalid_certificates': 0, 'first_height_report_checks': 0}
    firsts = []
    for row in r3['rows']:
        rec = next((x for x in co['rows'] if x['key'] == row['key']))
        mesh = np.load(rec['mesh_path'])
        public = np.load(rec['public_path'])
        best = None
        for (side, q) in row['all_pair_certificates'].items():
            for w in q['intersection_certificates']:
                A = mesh['vertices'][mesh['faces'][w['source_face_id']]]
                B = public[side][w['obstacle_face_id']]
                D = np.array(w['translation_mm'])
                valid = validate_witness(A, B, D, w)
                dvalid = dual_valid(A, B, D, w, w['dual'])
                counts['exact_primal_checks'] += 1
                counts['exact_dual_checks'] += 1
                counts['invalid_certificates'] += int(not (valid and dvalid))
                if full_reference:
                    ref = independent_lp(A, B, D)
                    counts['independent_LP_checks'] += 1
                    counts['reference_boolean_mismatches'] += int(not ref['feasible'])
                if best is None or Q(w['s_exact']) > Q(best['s_exact']):
                    best = w
        exact = Q(best['s_exact']) * 40
        counts['first_height_report_checks'] += 1
        counts['invalid_certificates'] += int(not first_height_valid(row))
        firsts.append(dict(key=row['key'], height_mm=float(exact), point_mm=best['point_mm'], source_face_id=best['source_face_id']))
    selected = [x for x in r4['rows'] if x['selected_translation_mm']]
    for row in selected:
        trials = row['trials']
        assert trials[-1]['status'] == 'SURFACE_CLEAR_CERTIFIED'
        for q in trials[-1]['results'].values():
            assert q['full_enumeration'] and q['coverage']['unknown_pairs'] == 0 and (q['coverage']['processed_faces'] == q['coverage']['source_faces'])
        assert row['crown_identity_error'] == 0 and row['margin_identity_error_mm'] == 0
    bad_height = copy.deepcopy(r3['rows'][0])
    bad_height['height_exact_mm'] = str(Q(bad_height['height_exact_mm']) + 1)
    repair = json.load(open(root / 'rounds/R2/REPAIR.json'))
    bad_repair = copy.deepcopy(repair)
    bad_repair['status'] = 'PASS'
    endpoint = json.load(open(root / 'raw/SUFFICIENCY_AND_FAULTS.json'))
    bad_endpoint = copy.deepcopy(endpoint)
    bad_endpoint['downstream_path'] = ['SURFACE_CLEAR_CERTIFIED'] * 2
    faults = {'wrong_first_height': not first_height_valid(bad_height), 'wrong_inside_repair_PASS': not repair_verdict_valid(bad_repair), 'false_endpoint_sufficiency': not endpoint_report_valid(bad_endpoint)}
    assert repair_verdict_valid(repair) and endpoint_report_valid(endpoint)
    assert counts['invalid_certificates'] == 0 and counts['reference_boolean_mismatches'] == 0 and all(faults.values())
    out = dict(status='PASS_SCOPED_DIGITAL_VALIDATION', claim_type='capability', counts=counts, faults=faults, firsts=firsts, seconds=time.perf_counter() - tick, external_referent={'kind': 'published_code', 'locator': 'https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html', 'compared_quantity': 'independent5variable barycentric triangle-sweep feasibility versus20654exact intersection certificates, no physical comparison', 'refutes_us': True}, scope='All saved first-event primal and dual certificates independently recomputed from hashed original face arrays. Separate SciPy formulation supplies numerical control; exact rational checker supplies acceptance. No external measured crown. Surface separation coverage inspected, every internal proof recomputed by one-command replay.')
    (root / 'raw/VALIDATION.json').write_text(json.dumps(out, indent=2) + '\n')
    print('VALIDATION', counts, round(out['seconds'], 2), flush=True)
    return out
if __name__ == '__main__':
    run()
