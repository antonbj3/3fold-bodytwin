"""Prospective comparison: no fitting, no mutation of frozen predictions."""
import argparse, csv
from common import *

def compare(path):
    frozen = read(ROOT / 'FROZEN_PREDICTIONS.json')
    if digest(frozen['payload']) != frozen['payload_sha256']:
        raise ValueError('prediction freeze corrupt')
    predictions = {}
    for case in read(ROOT / 'LAB_CASES_20.json')['cases']:
        file = Path(case['prediction_file'])
        if sha(file) != frozen['payload']['point_files'][str(file)]:
            raise ValueError('prediction drift')
        predictions[case['case_id']] = read(file)
    rows = []
    with open(path, newline='') as f:
        for obs in csv.DictReader(f):
            out = dict(observation=obs)
            required = ['case_id', 'library', 'axes', 'point_index', 'mesh_sha256', 'tool_library_sha256', 'observed_reachable', 'measured_assembly_and_frame', 'cam_locator']
            if any((not obs.get(k) for k in required)):
                out.update(status='UNKNOWN_MISSING_FIELDS')
                rows.append(out)
                continue
            try:
                case = predictions[obs['case_id']]
                if obs['mesh_sha256'] != case['source_sha256'] or obs['tool_library_sha256'] != frozen['payload']['tool_library_sha256'] or obs['measured_assembly_and_frame'] != 'true':
                    out.update(status='UNKNOWN_UNMATCHED_SCENE')
                    rows.append(out)
                    continue
                v = next((x for x in case['points_and_tools'] if x['library'] == obs['library'] and x['axes'] == int(obs['axes'])))
                p = v['records'][int(obs['point_index'])]
                if obs['observed_reachable'] not in ['true', 'false']:
                    raise ValueError('boolean must be true or false')
                pred = p['smallest_found_green_mm'] is not None
                obs_ok = obs['observed_reachable'] == 'true'
                out.update(predicted_witness=pred, observed_reachable=obs_ok, resolution='PER_POINT', status='FAIL_POSITIVE_WITNESS' if pred and (not obs_ok) else 'GRID_COVERAGE_GAP' if not pred and obs_ok else 'CONSISTENT_DIGITAL_OBSERVATION')
            except (KeyError, ValueError, StopIteration, IndexError) as e:
                out.update(status='UNKNOWN_INVALID_OBSERVATION', reason=str(e))
            rows.append(out)
    counts = {s: sum((r['status'] == s for r in rows)) for s in sorted({r['status'] for r in rows})}
    return dict(claim_type='capability', prediction_sha256=sha(ROOT / 'FROZEN_PREDICTIONS.json'), rows=rows, counts=counts, physical_validation='UNKNOWN' if not rows or any((r['status'].startswith('UNKNOWN') for r in rows)) else 'FAIL' if counts.get('FAIL_POSITIVE_WITNESS', 0) else 'CONSISTENT_NOT_YET_INDEPENDENTLY_REVIEWED', fit_performed=False)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('observations')
    ap.add_argument('--output', default='CAM_COMPARISON.json')
    args = ap.parse_args()
    dump(args.output, compare(args.observations))
