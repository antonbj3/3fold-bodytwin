"""Shared real observation model; optional external predicate annotations."""
import argparse
import json
import time
from pathlib import Path
from common import read, write
from measurement import answer

def demo(annotations=None):
    start = time.perf_counter()
    real = read('raw/RESULTS_R4.json')['real_shared_case']
    model = next((x for x in read('raw/PATIENT_MODELS_R2.json') if x['snapshot_row'] == real['snapshot_row']))
    obs = [o for o in read('raw/OBSERVATIONS_R2.json') if o['snapshot_row'] == real['snapshot_row']]
    a = annotations or {}
    if a:
        for key in ['patient_id', 'snapshot_row', 'FDI']:
            if a.get(key) != real[key]:
                raise ValueError('Annotation does not match frozen case/row/tooth')
        if not a.get('observer') or not a.get('source_locator'):
            raise ValueError('Independent annotation needs observer and locator')
        if a.get('pulp_measurement_method') in ['cold_only', 'cold', 'cold(-)'] and a.get('BD08', {}).get('pulp_state_compatible') is True:
            raise ValueError('Cold-only response cannot certify pulpal biological compatibility')
    panswer = answer('BD10', a.get('BD10', {}))
    local = a.get('local_periodontal_origin_confirmed')
    if local is not None and type(local) is not bool:
        raise ValueError('Local origin predicate must be literal Boolean')
    if local is not True:
        panswer['status'] = 'UNKNOWN_LOCAL_ETIOLOGY' if local is None else 'UNKNOWN_STAGE_LOCAL_EVIDENCE_EXCLUDED'
        panswer['possible_outcomes'] = ['STAGE_UNRESOLVED']
        panswer['missing_local_attribution'] = 'Independent exclusion of endodontic drainage, root damage and non-periodontal local causes'
    output = {'patient_id': real['patient_id'], 'snapshot_row': real['snapshot_row'], 'FDI': real['FDI'], 'CBCT_header_and_locator': model['volume'], 'observations': [{k: o[k] for k in ['observation_id', 'quantity', 'teeth', 'surface', 'resolution', 'match_text', 'locator']} for o in obs], 'BD10': panswer, 'BD08': answer('BD08', a.get('BD08', {})), 'local_periodontal_origin_confirmed': local, 'annotation_used': bool(a), 'independent_annotation_quality': 'UNKNOWN' if a else 'NOT_PERFORMED', 'claims': ['information_link', 'capability'], 'timescale': 'SIMULTANEOUS', 'complete_clinical_validation': False, 'clinical_recommendation': False, 'query_wall_s': time.perf_counter() - start}
    return output
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--annotations', type=Path)
    p.add_argument('--out', type=Path)
    args = p.parse_args()
    out = demo(json.loads(args.annotations.read_text()) if args.annotations else None)
    if args.out:
        args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n')
    else:
        print(json.dumps(out, ensure_ascii=False, indent=2))
