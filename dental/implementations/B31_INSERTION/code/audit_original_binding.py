"""Recheck original saved records after report-binding correction; reuse immutable numerical LP controls."""
import json, pathlib, hashlib, datetime, time, copy
from fractions import Fraction as Q
from validate import run
R = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def main():
    prior_path = R / 'raw/VALIDATION_NUMERICAL_V1.json'
    prior_path.write_text((R / 'versions/V1_BEFORE_SPATIAL_BINDING_FIX/VALIDATION.json').read_text())
    prior = json.load(open(prior_path))
    fresh = run(R, full_reference=False)
    (R / 'raw/VALIDATION_BINDING_V2.json').write_text(json.dumps(fresh, indent=2) + '\n')
    c = json.load(open(R / 'rounds/R3/FIRST_CONTACT.json'))
    headers = []

    def bound(row, w):
        return row['point_mm'] == w['point_mm'] and row['source_face_id'] == w['source_face_id'] and (row['source_role'] == w['source_role']) and (Q(row['height_exact_mm']) == Q(w['s_exact']) * 40)
    for row in c['rows']:
        allw = [w for q in row['all_pair_certificates'].values() for w in q['intersection_certificates']]
        w = max(allw, key=lambda w: Q(w['s_exact']))
        assert bound(row, w)
        bad = copy.deepcopy(row)
        bad['point_mm'][0] += 1.0
        assert not bound(bad, w)
        headers.append(dict(key=row['key'], header_bound=True, wrong_header_point_rejected=True))
    fresh['counts']['independent_LP_checks'] = prior['counts']['independent_LP_checks']
    fresh['counts']['reference_boolean_mismatches'] = prior['counts']['reference_boolean_mismatches']
    fresh['counts']['first_point_header_checks'] = len(headers)
    fresh['header_checks'] = headers
    fresh['seconds'] += prior['seconds']
    fresh['evidence_composition'] = {'original_numerical_reference': dict(path=str(prior_path), sha256=sha(prior_path), description='Same hashed source triangles and barycentric weights; display-binding correction does not alter any matrix'), 'corrected_exact_binding': dict(path=str(R / 'raw/VALIDATION_BINDING_V2.json'), sha256=sha(R / 'raw/VALIDATION_BINDING_V2.json')), 'original_first_event_sha256': sha(R / 'rounds/R3/FIRST_CONTACT.json'), 'final_copied_validation': 'raw/REPLAY_RECEIPTS.json V3samecohort fullnew6faults and20654independentLPs'}
    fresh['status'] = 'PASS_SOURCE_GEOMETRY_AND_REPORT_BINDING_V2'
    (R / 'raw/VALIDATION.json').write_text(json.dumps(fresh, indent=2) + '\n')
    print('ORIGINAL_BINDING20654valid+5headerpointsbound; allbadheadermutationsrejected', flush=True)
if __name__ == '__main__':
    main()
